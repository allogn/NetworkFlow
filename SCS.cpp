//
// Created by alvis on 20.03.16.
//

#include "SCS.h"

void SCS::startAugment(int max_length) {
    // Perform cost scaling phases
    IntVector path;
    BoolVector path_arc(_res_arc_num, false); // visited array
    int relabel_cnt = 0;
    int eps_phase_cnt = 0;
    for (; _epsilon >= 1; _epsilon = _epsilon < _alpha && _epsilon > 1 ?
                                     1 : _epsilon / _alpha) {
        ++eps_phase_cnt;

        // Price refinement heuristic
//        if (eps_phase_cnt >= PRICE_REFINEMENT_LIMIT) {
//          if (priceRefinement()) continue;
//        }

        // Initialize current phase
        initPhase();

        // Perform partial augment and relabel operations
        while (true) {

            // Select an active node (FIFO selection)
            while (_active_nodes.size() > 0 &&
                   _excess[_active_nodes.front()] <= 0) {
                _active_nodes.pop_front();
            }
            if (_active_nodes.size() == 0) break;
            int start = _active_nodes.front();

            // Find an augmenting path from the start node
            int tip = start;
            while (int(path.size()) < max_length && _excess[tip] >= 0) {
                int u;
                LargeCost rc, min_red_cost = std::numeric_limits<LargeCost>::max();
                LargeCost pi_tip = _pi[tip];
                int last_out = _first_out[tip + 1];
                for (int a = _next_out[tip]; a != last_out; ++a) {
                    if (_res_cap[a] > 0) {
                        u = _target[a];
                        rc = _cost[a] + pi_tip - _pi[u];
                        if (rc < 0) {
                            path.push_back(a);
                            _next_out[tip] = a;
                            if (path_arc[a]) {
                                goto augment;   // a cycle is found, stop path search
                            }
                            tip = u;
                            path_arc[a] = true;
                            goto next_step;
                        }
                        else if (rc < min_red_cost) {
                            min_red_cost = rc;
                        }
                    }
                }

                // Relabel tip node
                if (tip != start) {
                    int ra = _reverse[path.back()];
                    min_red_cost =
                            std::min(min_red_cost, _cost[ra] + pi_tip - _pi[_target[ra]]);
                }
                last_out = _next_out[tip];
                for (int a = _first_out[tip]; a != last_out; ++a) {
                    if (_res_cap[a] > 0) {
                        rc = _cost[a] + pi_tip - _pi[_target[a]];
                        if (rc < min_red_cost) {
                            min_red_cost = rc;
                        }
                    }
                }
                _pi[tip] -= min_red_cost + _epsilon;
                _next_out[tip] = _first_out[tip];
                ++relabel_cnt;

                // Step back
                if (tip != start) {
                    int pa = path.back();
                    path_arc[pa] = false;
                    tip = _source[pa];
                    path.pop_back();
                }

                next_step:;
            }

            // Augment along the found path (as much flow as possible)
            augment:
            Value delta;
            int pa, u, v = start;
            for (int i = 0; i != int(path.size()); ++i) {
                pa = path[i];
                u = v;
                v = _target[pa];
                path_arc[pa] = false;
                delta = std::min(_res_cap[pa], _excess[u]);
                _res_cap[pa] -= delta;
                _res_cap[_reverse[pa]] += delta;
                _excess[u] -= delta;
                _excess[v] += delta;
                if (_excess[v] > 0 && _excess[v] <= delta) {
                    _active_nodes.push_back(v);
                }
            }
            path.clear();
        }

    }

}

ProblemType SCS::init() {
    if (_res_node_num <= 1) return INFEASIBLE;

    // Check the sum of supply values
    _sum_supply = 0;
    for (int i = 0; i != _res_node_num; ++i) {
        _sum_supply += _supply[i];
    }
    if (_sum_supply > 0) return INFEASIBLE;

    // Check lower and upper bounds
    LEMON_DEBUG(checkBoundMaps(),
                "Upper bounds must be greater or equal to the lower bounds");


    // Initialize vectors
    for (int i = 0; i != _res_node_num; ++i) {
        _pi[i] = 0;
        _excess[i] = _supply[i];
    }

    // Remove infinite upper bounds and check negative arcs
    const Value MAX = std::numeric_limits<Value>::max();
    int last_out;
    if (_has_lower) {
        for (int i = 0; i != _res_node_num; ++i) {
            last_out = _first_out[i + 1];
            for (int j = _first_out[i]; j != last_out; ++j) {
                if (_forward[j]) {
                    Value c = _scost[j] < 0 ? _upper[j] : _lower[j];
                    if (c >= MAX) return UNBOUNDED;
                    _excess[i] -= c;
                    _excess[_target[j]] += c;
                }
            }
        }
    } else {
        for (int i = 0; i != _res_node_num; ++i) {
            last_out = _first_out[i + 1];
            for (int j = _first_out[i]; j != last_out; ++j) {
                if (_forward[j] && _scost[j] < 0) {
                    Value c = _upper[j];
                    if (c >= MAX) return UNBOUNDED;
                    _excess[i] -= c;
                    _excess[_target[j]] += c;
                }
            }
        }
    }
    Value ex, max_cap = 0;
    for (int i = 0; i != _res_node_num; ++i) {
        ex = _excess[i];
        _excess[i] = 0;
        if (ex < 0) max_cap -= ex;
    }
    for (int j = 0; j != _res_arc_num; ++j) {
        if (_upper[j] >= MAX) _upper[j] = max_cap;
    }

    // Initialize the large cost vector and the epsilon parameter
    _epsilon = 0;
    LargeCost lc;
    for (int i = 0; i != _res_node_num; ++i) {
        last_out = _first_out[i + 1];
        for (int j = _first_out[i]; j != last_out; ++j) {
            lc = _scost[j];//static_cast<LargeCost>(_scost[j]) * _res_node_num * _alpha;//todo
            _cost[j] = lc;
            if (lc > _epsilon) _epsilon = lc;
        }
    }
    _epsilon /= _alpha;

    // Initialize maps for Circulation and remove non-zero lower bounds
    IntVector cap(_graph.m);
    IntVector sup(_graph.m);

    for (uintT n = 0; n < _graph.n; n++) {
        sup[n] = _supply[_node_id[n]];
    }
    if (_has_lower) {
        for (uintT a = 0; a < _graph.m; a++) {
            int j = _arc_idf[a];
            Value c = _lower[j];
            cap[a] = _upper[j] - c;
            sup[_graph.E[a].fromid] -= c;
            sup[_graph.E[a].toid] += c;
        }
    } else {
        for (uintT a = 0; a < _graph.m; a++) {
            cap[a] = _upper[_arc_idf[a]];
        }
    }

    _sup_node_num = 0;
    for (uintT n = 0; n < _graph.n; n++) {
        if (sup[n] > 0) ++_sup_node_num;
    }

    // Find a feasible flow using Circulation
//        Circulation <Digraph, ConstMap<Arc, Value>, ValueArcMap, ValueNodeMap>
//                circ(_graph, low, cap, sup);
//        if (!circ.flowMap(flow).run()) return INFEASIBLE;

    // Set residual capacities and handle GEQ supply type
//    if (_sum_supply < 0) {
//        for (uintT a = 0; a < _graph.m; a++) {
//            Value fa = flow[a];
//            _res_cap[_arc_idf[a]] = cap[a] - fa;
//            _res_cap[_arc_idb[a]] = fa;
//            sup[_graph.E[a].fromid] -= fa;
//            sup[_graph.E[a].toid] += fa;
//        }
//        for (uintT n = 0; n < _graph.n; n++) {
//            _excess[_node_id[n]] = sup[n];
//        }
//        for (int a = _first_out[_root]; a != _res_arc_num; ++a) {
//            int u = _target[a];
//            int ra = _reverse[a];
//            _res_cap[a] = -_sum_supply + 1;
//            _res_cap[ra] = -_excess[u];
//            _cost[a] = 0;
//            _cost[ra] = 0;
//            _excess[u] = 0;
//        }
//    } else {
//        for (uintT a = 0; a < _graph.m; a++) {
//            Value fa = flow[a];
//            _res_cap[_arc_idf[a]] = cap[a] - fa;
//            _res_cap[_arc_idb[a]] = fa;
//        }
//        for (int a = _first_out[_root]; a != _res_arc_num; ++a) {
//            int ra = _reverse[a];
//            _res_cap[a] = 0;
//            _res_cap[ra] = 0;
//            _cost[a] = 0;
//            _cost[ra] = 0;
//        }
//    }

    // Initialize data structures for buckets
    _max_rank = _alpha * _res_node_num;
    _buckets.resize(_max_rank);
    _bucket_next.resize(_res_node_num + 1);
    _bucket_prev.resize(_res_node_num + 1);
    _rank.resize(_res_node_num + 1);

    return OPTIMAL;
}
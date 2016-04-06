//
// Created by alvis on 20.03.16.
//

#include "LSCS.h"

void LSCS::startAugment(long max_length) {


    IntVector path;
    BoolVector path_arc(_graph.m * 2, false); //reserve space for all possible edges
    long relabel_cnt = 0;
    long eps_phase_cnt = 0;
    for ( ; _epsilon >= 1; _epsilon = _epsilon < _alpha && _epsilon > 1 ?
                                      1 : _epsilon / _alpha )
    {
        ++eps_phase_cnt;

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
            long start = _active_nodes.front();

            // Find an augmenting path from the start node
            long tip = start;
            while (int(path.size()) < max_length && _excess[tip] >= 0) {
                long u;
                LargeCost rc, min_red_cost = std::numeric_limits<LargeCost>::max();
                LargeCost pi_tip = _pi[tip];

                for (long i = _next_out[tip]; i != _first_out[tip].size(); ++i) {
                    long a = _first_out[tip][i];
                    if (_res_cap[a] > 0) {
                        u = _target[a];
                        rc = _cost[a] + pi_tip - _pi[u];
                        if (rc < 0) {
                            path.push_back(a);
                            _next_out[tip] = i;
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


                //recalculate because some nodes before _next_out could have changed their values
                for (long i = 0; i < _next_out[tip]; ++i) {
                    long a = _first_out[tip][i];
                    if (_res_cap[a] > 0) {
                        rc = _cost[a] + pi_tip - _pi[_target[a]];
                        if (rc < min_red_cost) {
                            min_red_cost = rc;
                        }
                    }
                }

                // if next weight is greater than current minimum cost + target potential
                // cost + p(t) - p(s) > cost_min + p(t_min) - p(s)
                // cost + min_p(t) > cost; min_p(t) = 0 => cost - p(s) > min_red_cost
                while (QryCnt[tip] < _graph.fullE[tip].size() &&
                        _graph.E[_graph.fullE[tip][QryCnt[tip]]].weight + pi_tip <= min_red_cost) {
                    long i = _first_out[tip].size();
                    addArc(tip);
                    long a = _first_out[tip][i];
                    u = _target[a];
                    rc = _cost[a] + pi_tip - _pi[u];
                    if (rc < 0) {
                        path.push_back(a);
                        _next_out[tip] = i;
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

                _pi[tip] -= min_red_cost + _epsilon;
                _next_out[tip] = 0;
                ++relabel_cnt;

                // Step back
                if (tip != start) {
                    long pa = path.back();
                    path_arc[pa] = false;
                    tip = _source[pa];
                    path.pop_back();
                }

                next_step: ;
            }

            // Augment along the found path (as much flow as possible)
            augment:
            Value delta;
            long pa, u, v = start;
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

void LSCS::addArc(long fromid) {
    long eid = _graph.fullE[fromid][QryCnt[fromid]];
    QryCnt[fromid]++;

    long local_eid = _res_arc_num;
    //for spatial data this does not work
    _first_out[fromid].push_back(local_eid);
    long toid = _graph.get_pair(eid, fromid);
    _first_out[toid].push_back(local_eid + 1);//reverse edge
    _arc_idf.push_back(local_eid);
    _arc_idb.push_back(local_eid + 1);
    _forward.push_back(true);
    _forward.push_back(false);
    _source.push_back(fromid);
    _source.push_back(toid);
    _target.push_back(toid);
    _target.push_back(fromid);
    //skipped lower
    _upper.push_back(_graph.E[eid].capacity);
    _upper.push_back(_graph.E[eid].capacity);
    _scost.push_back(_graph.E[eid].weight);
    _scost.push_back(-_graph.E[eid].weight);
    LargeCost lc =
            static_cast<LargeCost>(_scost[local_eid]) * _res_node_num * _alpha; //COST MODIFICATION
    _cost.push_back(lc);
    _cost.push_back(-lc);

    _res_cap.push_back(_upper[local_eid]);
    _res_cap.push_back(0);
    _reverse.push_back(local_eid + 1);
    _reverse.push_back(local_eid);

    _res_arc_num += 2;
}

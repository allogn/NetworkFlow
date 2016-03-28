//
// Created by alvis on 20.03.16.
//

#include "SCS.h"

void SCS::startAugment(int max_length) {
    // Paramters for heuristics
//            const int PRICE_REFINEMENT_LIMIT = 2;
//            const double GLOBAL_UPDATE_FACTOR = 1.0;
//            const int global_update_skip = static_cast<int>(GLOBAL_UPDATE_FACTOR *
//                                                            (_res_node_num + _sup_node_num * _sup_node_num));
//            int next_global_update_limit = global_update_skip;

    // Perform cost scaling phases
//    IntVector path;
//    BoolVector path_arc(_res_arc_num, false);
//    int relabel_cnt = 0;
    int eps_phase_cnt = 0;
    for ( ; _epsilon >= 1; _epsilon = _epsilon < _alpha && _epsilon > 1 ?
                                      1 : _epsilon / _alpha )
    {
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
//            while (_active_nodes.size() > 0 &&
//                   _excess[_active_nodes.front()] <= 0) {
//                _active_nodes.pop_front();
//            }
//            if (_active_nodes.size() == 0) break;
//            int start = _active_nodes.front();


            // Dijkstra
            // put all excess nodes in first bucket of dijkstra
            // run dijkstra until deficit is found
            // increase potentials for all visited nodes
            // clear dijkstra buckets
            // increase flow on maximum value of flow possible (blocking flow)
            // assert target is not an excess (unit capacity) - remove source of a flow from active nodes
            // rerun dijkstra since all visited nodes were updated
            fHeap<LargeCost> dijkH;
            vector<pair<int,LargeCost>> visited;
            vector<int> parent(_res_node_num,-1); //parent node and an edge to a child
            LargeCost mindist[_res_node_num];
            fill(mindist, mindist+_res_node_num, std::numeric_limits<LargeCost>::max()); //those who were not popped yet, but are in the bucket already
            if (_active_nodes.size() == 0) break;
            for (deque<int>::iterator it = _active_nodes.begin(); it != _active_nodes.end(); ++it) {
                assert(_active_nodes.size() > 0);
                assert(_excess[*it] > 0);
                dijkH.enqueue(*it,0);
                mindist[*it] = 0;
            }

            LargeCost curcost;
            int u, tip;
            while (true) {
                dijkH.dequeue(tip, curcost);
                //take one node and add every neighbor to another bucket
                if (_excess[tip] < 0) goto increasePotentials;
                visited.push_back(pair<int,LargeCost>(tip,curcost));
                int last_out = (tip < _res_node_num -1) ? _first_out[tip+1]:_res_arc_num;
                LargeCost rc, l;
                LargeCost pi_tip = _pi[tip];
                for (int a = _first_out[tip]; a != last_out; a++) {
                    if (_res_cap[a] == 0) continue; // only feasible arcs
                    assert(_res_cap[a] > 0);
                    u = _target[a];
                    rc = _cost[a] - pi_tip + _pi[u]; //exactly +delta_pi because we increase potential
                    l = rc + _epsilon;
                    assert(l >= 0); //epsilon-optimality
                    //todo add visited array
                    if (mindist[u] > curcost + l) {
                        //check if exists in a heap
                        if (!dijkH.isExisted(u)) {
                            dijkH.enqueue(u, curcost + l);
                        } else {
                            dijkH.updatequeue(u, curcost + l);
                        }
                        parent[u] = a;
                        mindist[u] = curcost + l;
                    }
                }
            }
            increasePotentials:
            // increase potentials for all visited nodes
            for (int i = 0; i < visited.size(); i++) {
                assert(curcost >= visited[i].second);
                _pi[visited[i].first] += curcost - visited[i].second; //curbucket holds distance to deficit
            }

            //increase flow along path to deficit
            assert(_excess[tip] < 0); //we reached deficit
            while (parent[tip] != -1) {
                int edge = parent[tip];
                int parent_node = _source[edge];
                assert(_cost[edge] - _pi[parent_node] + _pi[tip] < 0); // admissible
                assert(_res_cap[edge] > 0);
                assert(_res_cap[_reverse[edge]] == 0);
                _res_cap[edge]--;
                _res_cap[_reverse[edge]]++;
                _excess[parent_node]--;
                _excess[tip]++;
                tip = parent_node;
            }

            assert(_excess[tip] >= 0);
            if (_excess[tip] == 0) {
                deque<int>::iterator it = std::find(_active_nodes.begin(),_active_nodes.end(),tip);
                _active_nodes.erase(it);
            }

            /* old implementation */
            // Find an augmenting path from the start node
//            int tip = start;
//            while (int(path.size()) < max_length && _excess[tip] >= 0) {
//                int u;
//                LargeCost rc, min_red_cost = std::numeric_limits<LargeCost>::max();
//                LargeCost pi_tip = _pi[tip];
//
//                int last_out = (tip < _res_node_num-1)?_first_out[tip+1]:_res_arc_num;
//                for (int a = _next_out[tip]; a != last_out; ++a) {
//                    if (_res_cap[a] > 0) {
//                        u = _target[a];
//                        rc = _cost[a] + pi_tip - _pi[u];
//                        if (rc < 0) {
//                            path.push_back(a);
//                            _next_out[tip] = a;
//                            if (path_arc[a]) {
//                                goto augment;   // a cycle is found, stop path search
//                            }
//                            tip = u;
//                            path_arc[a] = true;
//                            goto next_step;
//                        }
//                        else if (rc < min_red_cost) {
//                            min_red_cost = rc;
//                        }
//                    }
//                }
//
//                // Relabel tip node
//                last_out = _next_out[tip];
//                for (int a = _first_out[tip]; a != last_out; ++a) {
//                    if (_res_cap[a] > 0) {
//                        rc = _cost[a] + pi_tip - _pi[_target[a]];
//                        if (rc < min_red_cost) {
//                            min_red_cost = rc;
//                        }
//                    }
//                }
//                _pi[tip] -= min_red_cost + _epsilon;
//                _next_out[tip] = _first_out[tip];
//                ++relabel_cnt;
//
//                // Step back
//                if (tip != start) {
//                    int pa = path.back();
//                    path_arc[pa] = false;
//                    tip = _source[pa];
//                    path.pop_back();
//                }
//
//                next_step: ;
//            }

            // Augment along the found path (as much flow as possible)
//            augment:
//            Value delta;
//            int pa, u, v = start;
//            for (int i = 0; i != int(path.size()); ++i) {
//                pa = path[i];
//                u = v;
//                v = _target[pa];
//                path_arc[pa] = false;
//                delta = std::min(_res_cap[pa], _excess[u]);
//                _res_cap[pa] -= delta;
//                _res_cap[_reverse[pa]] += delta;
//                _excess[u] -= delta;
//                _excess[v] += delta;
//                if (_excess[v] > 0 && _excess[v] <= delta) {
//                    _active_nodes.push_back(v);
//                }
//            }
//            path.clear();

            // Global update heuristic
//          if (relabel_cnt >= next_global_update_limit) {
//            globalUpdate();
//                    next_global_update_limit += global_update_skip;
//          }
        }

    }
}
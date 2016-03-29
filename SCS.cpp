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

            //iteration reset
            dijkH.clear();
            updateH.clear();
            watched.resize(_res_node_num,false);
            parent.resize(_res_node_num,-1);
            mindist.resize(_res_node_num,std::numeric_limits<LargeCost>::max());

            if (_active_nodes.size() == 0) break;
//            for (deque<int>::iterator it = _active_nodes.begin(); it != _active_nodes.end(); ++it) {
//                assert(_active_nodes.size() > 0);
//                assert(_excess[*it] > 0);
//                dijkH.enqueue(*it,0);
//                mindist[*it] = 0;
//            }
            int q = _active_nodes.front();
            dijkH.enqueue(q,0);
            mindist[q] = 0;
            watched[q] = 0;
            heap_checkAndUpdateEdgeMin(globalH,q);

            // main loop (process_id)
            LargeCost curcost;
            int u, tip = -1;
            while (tip == -1 || (globalH.size() > 0 && mindist[tip]>globalH.getTopValue() - taumax)) {

                do //  || (globalH.size() > 0 && heap_getTopValue(dijkH)>heap_getTopValue(globalH))) //dijkH check because next in globalH can be not from current "thread", but must be added because it will be needed later
                {
                    int fromid;
                    LargeCost tmp;
                    globalH.dequeue(fromid, tmp);

                    //create new edge
                    uintT eid = _graph.fullE[fromid][QryCnt[fromid]-1];
                    //todo rewrite adding new nodes
                    uintT local_eid = _res_arc_num;
                    //todo why flow = 1??

                    //for spatial data this does not work
                    _first_out[fromid].push_back(local_eid);
                    int toid = _graph.get_pair(eid, fromid);
                    _first_out[toid].push_back(local_eid+1);//reverse edge
                    _arc_idf[eid] = local_eid;
                    _arc_idb[eid] = local_eid+1;
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
                    _res_cap.push_back(_upper[local_eid]);
                    _res_cap.push_back(0);

                    //todo reserve space for all vectors in init

                    _reverse.push_back(local_eid+1);
                    _reverse.push_back(local_eid);

                    heap_checkAndUpdateEdgeMin(globalH,fromid); //MUST BE HERE
                    _res_arc_num += 2;

                    /*
                     * Dijkstra Updates
                     */
                    updateHeaps(eid, fromid, toid);

                    int curid;
                    while (updateH.size() > 0) {
                        updateH.dequeue(curid, tmp);
                        //            isUpdated[curid] = 0;
                        for (vector<int>::iterator it = _first_out[curid].begin(); it < _first_out[curid].end(); it++) {
                            toid = _target[*it];
                            updateHeaps(eid, curid, toid);
                        }
                    }
                    /*
                     * End of dijkstra updates. If switch off => add clearing dijkstra in dijksra beginning
                     */
                } while (dijkH.size() == 0 && tip == -1); //todo why here tip==-1 needed??

                //running Dijkstra
                assert(dijkH.size() > 0);
                tip = dijkH.getTopIdx();
                while (true) {
                    //take one node and add every neighbor to another bucket
                    if (_excess[tip] < 0) goto increasePotentials;
                    visited.push_back(pair<int, LargeCost>(tip, curcost));
                    LargeCost rc, l;
                    LargeCost pi_tip = _pi[tip];
                    for (vector<int>::iterator a = _first_out[tip].begin(); a != _first_out[tip].end(); a++) {
                        if (_res_cap[*a] == 0) continue; // only feasible arcs
                        assert(_res_cap[*a] > 0);
                        u = _target[*a];
                        rc = _cost[*a] - pi_tip + _pi[u]; //exactly +delta_pi because we increase potential
                        l = rc + _epsilon;
                        assert(l >= 0); //epsilon-optimality
                        //todo add visited array
                        if (mindist[u] > curcost + l) {
                            //check if exists in a heap
                            if (!dijkH.isExisted(u)) {
                                dijkH.enqueue(u, curcost + l);
                            } else {
                                dijkH.updatequeue(u, curcost + l);
                                heap_checkAndUpdateEdgeMin(globalH, u);
                            }
                            parent[u] = *a;
                            mindist[u] = curcost + l;
                        }
                    }
                    heap_checkAndUpdateEdgeMin(globalH, tip);
                    watched[tip] = true;
                    dijkH.dequeue();
                    if (dijkH.size() == 0) {
                        tip = -1;
                        break;
                    }
                    tip = dijkH.getTopIdx();
                }

                //TODO keep up with taumax
                int gettaumax = 0;
                for (int i = 0; i<_res_node_num; i++)
                {
                    if (mindist[i]<globalH.getTopValue() && gettaumax < _pi[i])
                    {
                        gettaumax = _pi[i];
                    }
                }
                taumax = gettaumax;

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
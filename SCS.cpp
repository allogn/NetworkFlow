//
// Created by alvis on 20.03.16.
//

#include "SCS.h"

void SCS::startAugment(long max_length) {
    // Paramters for heuristics
//            const long PRICE_REFINEMENT_LIMIT = 2;
//            const double GLOBAL_UPDATE_FACTOR = 1.0;
//            const long global_update_skip = static_cast<long>(GLOBAL_UPDATE_FACTOR *
//                                                            (_res_node_num + _sup_node_num * _sup_node_num));
//            long next_global_update_limit = global_update_skip;

    // Perform cost scaling phases
//    longVector path;
//    BoolVector path_arc(_res_arc_num, false);
//    long relabel_cnt = 0;
    _active_nodes.clear();
    for (long u = 0; u != _res_node_num; ++u) {
        if (_excess[u] > 0) _active_nodes.push_back(u);
    }

    long eps_phase_cnt = 0;
    for (; _epsilon >= 1 || _active_nodes.size() > 0; _epsilon = _epsilon < _alpha && _epsilon > 1 ?
                                     1 : _epsilon / _alpha) {
        ++eps_phase_cnt;
        bool newEdgeAdded = false;
        // Price refinement heuristic
//        if (eps_phase_cnt >= PRICE_REFINEMENT_LIMIT) {
//          if (priceRefinement()) continue;
//        }
//        cout << "epsilon " << _epsilon << endl;
        // Initialize current phase
        initPhase();

        // Perform partial augment and relabel operations
        while (true) {
            // Select an active node (FIFO selection)
            _active_nodes.clear();
            for (long u = 0; u != _res_node_num; ++u) {
                if (_excess[u] > 0) _active_nodes.push_back(u);
            }
            if (_active_nodes.size() == 0) break;
            long start = _active_nodes.front();
//            cout << "starting " << start << endl;

            // Dijkstra
            // put all excess nodes in first bucket of dijkstra
            // run dijkstra until deficit is found
            // increase potentials for all visited nodes
            // clear dijkstra buckets
            // increase flow on maximum value of flow possible (blocking flow)
            // assert target is not an excess (unit capacity) - remove source of a flow from active nodes
            // rerun dijkstra since all visited nodes were updated

            fHeap<LargeCost> dijkH;
            fHeap<LargeCost> globalH;
            fHeap<LargeCost> updateH;
            vector<long> visited;
            vector<long> parent(_res_node_num, -1); //parent node and an edge to a child
            LargeCost mindist[_res_node_num];

            long w = _graph.get_next_neighbour_weight(start) * _res_node_num * _alpha;
            if (w >= 0)
                globalH.enqueue(start, w);

            //do not update every time

            LargeCost curcost;
            long u, tip;
            long taumax = 0;

            dijkH.clear();
            updateH.clear();
            dijkH.enqueue(start,0);
            visited.clear();
            std::fill(parent.begin(),parent.end(),-1);
            fill(mindist, mindist + _res_node_num,
                 std::numeric_limits<LargeCost>::max()); //those who were not popped yet, but are in the bucket already
            mindist[start] = 0;

            do {

                bool once = false;
                if  (globalH.size() > 0) {//}} && (!once || dijkH.getTopValue() > globalH.getTopValue() - taumax)) { //todo add this
                    newEdgeAdded = true; // for new epsilon calculation
                    once = true;
                    long fromid;
                    LargeCost tmp;
                    globalH.dequeue(fromid,tmp);
                    long eid = _graph.insert_nn_to_edge_list(fromid);

                    long w = _graph.get_next_neighbour_weight(fromid) * _res_node_num * _alpha;
                    if (w >= 0) {
                        globalH.enqueue(fromid, mindist[fromid] + w);
                    }

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

//                    cout << "added " << fromid << "->" << toid << endl;
                    //update dijkH heap to continue dijkstra from the polong it ended
                    LargeCost rc = _cost[local_eid] - _pi[fromid] + _pi[toid]; //exactly +delta_pi because we increase potential
                    LargeCost l = rc + _epsilon;
                    assert(l >= 0); //epsilon-optimality
                    assert(mindist[fromid] < std::numeric_limits<LargeCost>::max()); //because if it was placed in global => it was visited
                    if (mindist[toid] > mindist[fromid] + l) {
//                        if (!dijkH.isExisted(fromid)) {
//                            dijkH.enqueue(fromid, mindist[fromid]);
//                        }// else {
                         //   dijkH.updatequeue(fromid, mindist[fromid]);
                       // }
                        //check if exists in a heap
                        if (!dijkH.isExisted(toid)) {
                            dijkH.enqueue(toid, mindist[fromid] + l);
                        } else {
                            dijkH.updatequeue(toid, mindist[fromid] + l);
                        }
                        parent[toid] = local_eid;
                        mindist[toid] = mindist[fromid] + l;
                        if (globalH.isExisted(toid)) {
                            assert(_graph.get_next_neighbour_weight(toid) >= 0); //if in heap - then this is true
                            globalH.updatequeue(toid, mindist[toid] + _graph.get_next_neighbour_weight(toid) * _res_node_num * _alpha);
                        }
//                        cout << "updated distance to " << toid << ": " << mindist[toid] << endl;
                    }
                }

                while (true) {
                    //take one node and add every neighbor to another bucket

                    tip = dijkH.getTopIdx();
//                    cout << "tip " << tip << endl;
                    curcost = dijkH.getTopValue();
                    if (_excess[tip] < 0) {
//                        cout << "end of dijkstra (negative excess)" << endl;
                        goto checkEdges;
                    }
                    visited.push_back(tip);
                    dijkH.dequeue();
                    if (!globalH.isExisted(tip)) {
                        long w = _graph.get_next_neighbour_weight(tip) * _res_node_num * _alpha;
                        if (w >= 0) {
                            globalH.enqueue(tip, mindist[tip] + w);
                        }
                    }

                    LargeCost rc, l;
                    LargeCost pi_tip = _pi[tip];
                    assert(curcost == mindist[tip]);
                    for (vector<long>::iterator a = _first_out[tip].begin(); a != _first_out[tip].end(); a++) {
                        if (_res_cap[*a] == 0) continue; // only feasible arcs
                        assert(_res_cap[*a] > 0);
                        u = _target[*a];
                        rc = _cost[*a] - pi_tip + _pi[u]; //exactly +delta_pi because we increase potential
                        l = rc + _epsilon;
                        assert(l >= 0); //epsilon-optimality
                        if (mindist[u] > curcost + l) {
                            //check if exists in a heap
                            if (!dijkH.isExisted(u)) {
                                dijkH.enqueue(u, curcost + l);
                            } else {
                                dijkH.updatequeue(u, curcost + l);
                            }
                            parent[u] = *a;
                            mindist[u] = curcost + l;
                            if (globalH.isExisted(u)) {
                                assert(_graph.get_next_neighbour_weight(u) >= 0); //if in heap - then this is true
                                globalH.updatequeue(u, mindist[u] + _graph.get_next_neighbour_weight(u) * _res_node_num * _alpha);
                            }
                        }
                    }

                    if (dijkH.size() == 0) {
                        //no path
                        tip = -1;
                        goto checkEdges;
                    }
                }
                checkEdges:;

                if (globalH.size() > 0) {
                    long tv = globalH.getTopValue();
                    for (long i = 0; i < visited.size(); i++) {
                        long nodeid = visited[i];
                        if (mindist[nodeid] < tv && _pi[nodeid] > taumax) {
                            taumax = _pi[nodeid];
                        }
                    }
                }

                assert((globalH.size() == 0 && tip != -1) || globalH.size() > 0);
            } while (tip == -1 || (globalH.size() > 0 && mindist[tip] > globalH.getTopValue() - taumax + _epsilon/ _alpha));


#ifndef NDEBUG
            fHeap<LargeCost> dijkH2;
            dijkH2.enqueue(start,0);
            vector<long> visited2;
            vector<long> parent2(_res_node_num,-1);
            LargeCost mindist2[_res_node_num];
            fill(mindist2, mindist2 + _res_node_num,
                 std::numeric_limits<LargeCost>::max()); //those who were not popped yet, but are in the bucket already
            mindist2[start] = 0;
            LargeCost curcost2;
            long tip2;

            while (true) {
                //take one node and add every neighbor to another bucket
                assert(dijkH2.size() > 0);
                tip2 = dijkH2.getTopIdx();
                curcost2 = dijkH2.getTopValue();
                if (_excess[tip2] < 0) {
                    goto checkEdges2;
                }
                dijkH2.dequeue();
                visited2.push_back(tip2);
                LargeCost rc, l;
                LargeCost pi_tip = _pi[tip2];
                assert(curcost2 == mindist2[tip2]);
                for (vector<long>::iterator a = _first_out[tip2].begin(); a != _first_out[tip2].end(); a++) {
                    if (_res_cap[*a] == 0) continue; // only feasible arcs
                    assert(_res_cap[*a] > 0);
                    u = _target[*a];
                    rc = _cost[*a] - pi_tip + _pi[u]; //exactly +delta_pi because we increase potential
                    l = rc + _epsilon;
                    assert(l >= 0); //epsilon-optimality
                    if (mindist2[u] > curcost2 + l) {
                        //check if exists in a heap
                        if (!dijkH2.isExisted(u)) {
                            dijkH2.enqueue(u, curcost2 + l);
                        } else {
                            dijkH2.updatequeue(u, curcost2 + l);
                        }
                        parent2[u] = *a;
                        mindist2[u] = curcost2 + l;
                    }
                }
            }
            checkEdges2:;
            assert(curcost == curcost2);
            for (long k = 0; k < _res_node_num; k++) {
                if (curcost > mindist[k])
                    assert(mindist[k] == mindist2[k]);
            }
#endif

            //todo add visited here as well
            for (long i = 0; i < _res_node_num; i++) {
                if (curcost > mindist[i])
                {
                    _pi[i] += curcost - mindist[i];
                } //curbucket holds distance to deficit
            }
            assert(test_epsilon());
//            cout << "increasing potentials" << endl;
            //increase flow along path to deficit
            assert(_excess[tip] < 0); //we reached deficit
            while (parent[tip] != -1) {
                long edge = parent[tip];
                long parent_node = _source[edge];
                assert(_cost[edge] - _pi[parent_node] + _pi[tip] < 0); // admissible
                assert(_res_cap[edge] > 0);
                assert(_res_cap[_reverse[edge]] == 0);
                _res_cap[edge]--;
                _res_cap[_reverse[edge]]++;
                _excess[parent_node]--;
                _excess[tip]++;
                tip = parent_node;
            }
            assert(test_epsilon());
            assert(tip == start);
        }
        assert(test_epsilon());

        _active_nodes.clear();
        for (long u = 0; u != _res_node_num; ++u) {
            if (_excess[u] > 0) _active_nodes.push_back(u);
        }

//        if (newEdgeAdded) {
//            long max_w = 0;
//            for (long i = 0; i < _res_node_num; i++) {
//                long nw = _graph.get_next_neighbour_weight(i);
//                max_w = (nw > max_w) ? nw : max_w;
//            }
//            max_w = max_w * _alpha * _res_node_num;
//            _epsilon = (_epsilon > max_w) ? _epsilon : max_w;
//        }

    }

}
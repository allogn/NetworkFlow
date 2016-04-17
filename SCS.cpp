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
            vector<long> visited;
            vector<long> parent(_res_node_num, -1); //parent node and an edge to a child
            LargeCost mindist[_res_node_num];
            if (QryCnt[start] < _graph.fullE[start].size())
                globalH.enqueue(start, 0);

            //do not update every time

            LargeCost curcost;
            long u, tip;
            long taumax = 0;

            dijkH.clear();
            dijkH.enqueue(start,0);
            visited.clear();
            std::fill(parent.begin(),parent.end(),-1);
            fill(mindist, mindist + _res_node_num,
                 std::numeric_limits<LargeCost>::max()); //those who were not popped yet, but are in the bucket already
            mindist[start] = 0;

            do {

                if (globalH.size() > 0) {
                    long fromid;
                    LargeCost tmp;
                    globalH.dequeue(fromid,tmp);
                    long eid = _graph.fullE[fromid][QryCnt[fromid]];
                    QryCnt[fromid]++;
                    if (QryCnt[fromid] < _graph.fullE[fromid].size()) {
                        globalH.enqueue(fromid, mindist[fromid] + _graph.E[_graph.fullE[fromid][QryCnt[fromid]]].weight);
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

                    // calculate new epsilon based on two potential values and cost (note the minus sign!)
//                    LargeCost new_epsilon = std::max(-(lc + _pi[toid] - _pi[fromid]), _epsilon);
//                    LargeCost new_epsilon = std::max(lc, _epsilon);
//                    new_epsilon = std::max(-(lc + _pi[toid] - _pi[fromid]), new_epsilon); //remove?
//                    if (new_epsilon != _epsilon) cout << "new epsilon " << new_epsilon << endl;
//                    _epsilon = new_epsilon; //todo maybe this can be modified
//todo add while tophoeap differ condition
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
                    //todo add visited array
                    assert(mindist[fromid] < std::numeric_limits<LargeCost>::max()); //because if it was placed in global => it was visited
                    if (mindist[toid] > mindist[fromid] + l) {
                        //check if exists in a heap
                        if (!dijkH.isExisted(toid)) {
                            dijkH.enqueue(toid, mindist[fromid] + l);
                        } else {
                            dijkH.updatequeue(toid, mindist[fromid] + l);
                        }
                        parent[toid] = local_eid;
                        mindist[toid] = mindist[fromid] + l;
                        if (globalH.isExisted(toid)) {
                            assert(QryCnt[toid] < _graph.fullE[toid].size()); //if in heap - then this is true
                            globalH.updatequeue(toid, mindist[toid] + _graph.E[_graph.fullE[toid][QryCnt[toid]]].weight);
                        }
//                        cout << "updated distance to " << toid << ": " << mindist[toid] << endl;
                    }
                }

                while (true) {
                    //take one node and add every neighbor to another bucket

                    if (dijkH.size() == 0) {
                        //no path
                        tip = -1;
                        goto checkEdges;
                    }

                    tip = dijkH.getTopIdx();
//                    cout << "tip " << tip << endl;
                    curcost = dijkH.getTopValue();
                    if (_excess[tip] < 0) {
//                        cout << "end of dijkstra (negative excess)" << endl;
                        goto checkEdges;
                    }
                    dijkH.dequeue();
                    if (!globalH.isExisted(tip) && QryCnt[tip] < _graph.fullE[tip].size()) {
                        globalH.enqueue(tip, mindist[tip] + _graph.E[_graph.fullE[tip][QryCnt[tip]]].weight);
                    }
                    visited.push_back(tip);
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
                        //todo add visited array
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
                                assert(QryCnt[u] < _graph.fullE[u].size()); //if in heap - then this is true
                                globalH.updatequeue(u, mindist[u] + _graph.E[_graph.fullE[u][QryCnt[u]]].weight);
                            }
                        }
                    }
                }
                //todo maybe we can clear dijkstra since everything else is definetly farther than corrent result?
                checkEdges:;

                if (globalH.size() > 0) {
                    for (long i = 0; i < _res_node_num; i++) {
                        if (mindist[i] < globalH.getTopValue() && _pi[i] > taumax) {
                            taumax = _pi[i];
                        }
                    }
                }

                assert((globalH.size() == 0 && tip != -1) || globalH.size() > 0);
            } while (tip == -1 || (globalH.size() > 0 && mindist[tip] > globalH.getTopValue() - taumax));


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
//                    cout << "tip " << tip << endl;
                curcost2 = dijkH2.getTopValue();
                if (_excess[tip2] < 0) {
//                        cout << "end of dijkstra (negative excess)" << endl;
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
                    //todo add visited array
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


            // increase potentials for all visited nodes
//            for (long i = 0; i < visited.size(); i++) {
//                assert(QryCnt[visited[i]] == _graph.fullE[visited[i]].size());
//                if (curcost > mindist[visited[i]])
//                {
//                    _pi[visited[i]] += curcost - mindist[visited[i]];
//                } //curbucket holds distance to deficit
//            }

//            long totalFull = 0;
//            double others;
            for (long i = 0; i < _res_node_num; i++) {
                if (curcost > mindist[i])
                {
                    _pi[i] += curcost - mindist[i];
                } //curbucket holds distance to deficit

//                if (_graph.fullE[i].size() > 0) {
//                    if (_first_out[i].size() == _graph.fullE[i].size()) {
//                        totalFull++;
//                    } else {
//                        others += (double)_first_out[i].size()*100./(double)_graph.fullE[i].size();
//                    }
//                }
            }
//            others /= _res_node_num - totalFull;
//            cout << "Full nodes: " << totalFull << " out of " << _res_node_num << endl;
//            cout << "Others filled for " << others*100. << "%" << endl;
//            assert(test_epsilon());

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
    }
}
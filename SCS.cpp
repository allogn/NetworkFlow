//
// Created by alvis on 20.03.16.
//

#include "SCS.h"

void SCS::startAugment(int max_length) {
    watched.resize(_res_node_num, false);
    parent.resize(_res_node_num, -1);
    mindist.resize(_res_node_num, std::numeric_limits<LargeCost>::max());
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
    for (int u = 0; u != _res_node_num; ++u) {
        if (_excess[u] > 0) _active_nodes.push_back(u);
    }
    for (; _epsilon >= 1 || _active_nodes.size() > 0; _epsilon = _epsilon < _alpha && _epsilon > 1 ?
                                                                 1 : _epsilon / _alpha) {

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
            //iteration reset
            dijkH.clear();
            globalH.clear();
            updateH.clear();
            std::fill(watched.begin(), watched.end(), false);
            std::fill(parent.begin(), parent.end(), -1);
            std::fill(mindist.begin(), mindist.end(), std::numeric_limits<LargeCost>::max());
            visited.clear();

            if (_active_nodes.size() == 0) {
                break;
            }
//            for (deque<int>::iterator it = _active_nodes.begin(); it != _active_nodes.end(); ++it) {
//                assert(_active_nodes.size() > 0);
//                assert(_excess[*it] > 0);
//                dijkH.enqueue(*it,0);
//                mindist[*it] = 0;
//            }
            int q = _active_nodes.front();
            dijkH.enqueue(q, 0);
            mindist[q] = 0;
            watched[q] = 0;
//            assert(globalH.size() == 0);
//            if (QryCnt[q] < _graph.fullE[q].size()) {
//                globalH.enqueue(q, _graph.E[_graph.fullE[q][QryCnt[q]]].weight);
//            }
//
//            // main loop (process_id)
            int u, tip = -1;
//            while (tip == -1 || (globalH.size() > 0)) {//} && mindist[tip] > globalH.getTopValue() - taumax)) {
//
//                // adding new edge
////                do  //  || (globalH.size() > 0 && heap_getTopValue(dijkH)>heap_getTopValue(globalH))) //dijkH check because next in globalH can be not from current "thread", but must be added because it will be needed later
////                {
//                    int fromid;
//                    LargeCost tmp;
//                    assert(globalH.size() > 0);
//                    globalH.dequeue(fromid, tmp);
//                    assert(tmp > 0);
//                    assert(!globalH.isExisted(fromid));
//
//                    //create new edge
//                    uintT eid = _graph.fullE[fromid][QryCnt[fromid]];//todo check this in SIA
//
//                    QryCnt[fromid]++;
//                    if (QryCnt[fromid] < _graph.fullE[fromid].size()) {
//                        globalH.enqueue(fromid, _graph.E[_graph.fullE[fromid][QryCnt[fromid]]].weight);
//                    }
//
//                    //todo rewrite adding new nodes
//                    uintT local_eid = _res_arc_num;
//                    //todo why flow = 1??
//
//                    //for spatial data this does not work
//                    _first_out[fromid].push_back(local_eid);
//                    int toid = _graph.get_pair(eid, fromid);
//                    _first_out[toid].push_back(local_eid + 1);//reverse edge
//                    _arc_idf.push_back(local_eid);
//                    _arc_idb.push_back(local_eid + 1);
//                    _forward.push_back(true);
//                    _forward.push_back(false);
//                    _source.push_back(fromid);
//                    _source.push_back(toid);
//                    _target.push_back(toid);
//                    _target.push_back(fromid);
//                    //skipped lower
//                    _upper.push_back(_graph.E[eid].capacity);
//                    _upper.push_back(_graph.E[eid].capacity);
//                    _scost.push_back(_graph.E[eid].weight);
//                    _scost.push_back(-_graph.E[eid].weight);
//                    LargeCost lc =
//                            static_cast<LargeCost>(_scost[local_eid]) * _res_node_num * _alpha; //COST MODIFICATION
//                    _cost.push_back(lc);
//                    _cost.push_back(-lc);
//
//                    // calculate new epsilon based on two potential values and cost (note the minus sign!)
//                    LargeCost new_epsilon = std::max(-(lc + _pi[toid] - _pi[fromid]), _epsilon);
//
//                    _epsilon = new_epsilon; //todo maybe this can be modified
//                    _res_cap.push_back(_upper[local_eid]);
//                    _res_cap.push_back(0);
//
//                    //todo reserve space for all vectors in init
//
//                    _reverse.push_back(local_eid + 1);
//                    _reverse.push_back(local_eid);
//
//                    _res_arc_num += 2;
//
//                    /*
//                     * Dijkstra Updates
//                     */
//                    assert(_target[local_eid] == toid && _source[local_eid] == fromid);
//                    //watched: case when prefinal node : not all neighbours are considered => not watched, not dequeued from dijkH todo add watched to algorithm description
//                    assert(_res_cap[local_eid] > 0);
//                    // but in globalH and in DijkH (!). so, if in dijkH => everything is fine (will be updated later)
//                    if ((watched[fromid] == true) && updateMinDist(local_eid,fromid,toid)) {
//                        heap_checkAndUpdateEdgeMin(globalH, toid);
//                        if (!dijkH.isExisted(toid)) {
//                            if (watched[toid] == 1) heap_checkAndUpdateMin(updateH, toid, mindist[toid]);
//                            else dijkH.enqueue(toid, mindist[toid]);
//                        } else {
//                            dijkH.updatequeue(toid, mindist[toid]);
//                        }
//                    }
//
//                    int curid;
//                    while (updateH.size() > 0) {
//                        updateH.dequeue(curid, tmp);
//                        //            isUpdated[curid] = 0;
//                        for (vector<int>::iterator it = _first_out[curid].begin(); it < _first_out[curid].end(); it++) {
//                            toid = _target[*it];
//                            assert(_source[*it] == curid);
//                            assert(_target[*it] == toid && _source[*it] == curid);
//                            // but in globalH and in DijkH (!). so, if in dijkH => everything is fine (will be updated later)
//                            if (_res_cap[*it] == 0 && watched[fromid] == true && updateMinDist(*it,curid,toid)) {
//                                heap_checkAndUpdateEdgeMin(globalH, toid);
//                                //no isUpdated because enheap only if updated dist
//                                if (!dijkH.isExisted(toid)) {
//                                    if (watched[toid] == 1) heap_checkAndUpdateMin(updateH, toid, mindist[toid]);
//                                    else dijkH.enqueue(toid, mindist[toid]);
//                                } else {
//                                    dijkH.updatequeue(toid, mindist[toid]);
//                                }
//                            }
//                        }
//                    }
                    /*
                     * End of dijkstra updates. If switch off => add clearing dijkstra in dijksra beginning
                     */
//                } while (dijkH.size() == 0); //todo why here tip==-1 needed??

//                assert(test_dijkstra_heap(q));

                // at least one time must be added because otherwise infinite loop:
                // 0 -> 8, 8 is still in dijkH because final node
                // but threshold is not satisfied, so trying to add mode edges, but skipping the step because non-null dijkstra

                //running Dijkstra
                assert(dijkH.size() > 0);
                tip = dijkH.getTopIdx();
                while (true) {
                    //take one node and add every neighbor to another bucket
                    if (_excess[tip] < 0) {
//                        dijkH.dequeue();
                        break; // continue to taumax calculation and checking threshold
                    }
                    visited.push_back(tip);
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
                        if (mindist[u] > mindist[tip] + l) {
                            //check if exists in a heap
                            if (!dijkH.isExisted(u)) {
                                dijkH.enqueue(u, mindist[tip] + l);
                            } else {
                                dijkH.updatequeue(u, mindist[tip] + l);
                                heap_checkAndUpdateEdgeMin(globalH, u);
                            }
                            parent[u] = *a;
                            mindist[u] = mindist[tip] + l;
                        }
                    }
                    heap_checkAndUpdateEdgeMin(globalH, tip);
                    watched[tip] = true;
                    dijkH.dequeue();
                    if (dijkH.size() == 0) {
                        tip = -1;
                        cout << "break because no way" << endl;
                        break;
                    }
                    tip = dijkH.getTopIdx();
                }
//                assert(test_if_shortest(q, tip));

                /*
                 * the problem: 0->8 -- 8 is the best, but new edge is added. 0->10 added, but 0->8 is still the best
                 * Q: is this correct that new added edge does not change the best option?
                 * dijkstra is still correct, but new reached node must not be better than current one
                 * solution: never delete last node from the dijkstra! (see line 193)
                 */

                //TODO keep up with taumax
//                int gettaumax = 0;
//                for (int i = 0; i < _res_node_num; i++) {
//                    if (mindist[i] < globalH.getTopValue() && gettaumax < _pi[i]) {
//                        gettaumax = _pi[i];
//                    }
//                }
//                taumax = gettaumax;

//            }
//            increasePotentials:
            // increase potentials for all visited nodes

            /*
             * before globalH it used to have dist to tip the largest between all
             * after it is not the case because we can add an edge to a deficit that will not update
             * most of other distances to a farthest points because it has no path to them.
             * at the same time new edge goes, for example, directly to a deficit. so, here we must update
             * only those who has small enough path.
             *
             * here is a problem of redundant work: we traverse too much al loose that information
             */
//            cout << endl;
            for (int i = 0; i < visited.size(); i++) {
                if (mindist[tip] >= mindist[visited[i]]) {
//                    assert(mindist[tip] >= visited[i].second);
                    _pi[visited[i]] += mindist[tip] - mindist[visited[i]]; //curbucket holds distance to deficit
                }
            }

            //increase flow along path to deficit
            assert(_excess[tip] < 0); //we reached deficit
            while (parent[tip] != -1) {
                int edge = parent[tip];
                assert(edge < _target.size());
                int parent_node = _source[edge];
                assert(_target[edge] == tip);
                assert(_cost[edge] - _pi[parent_node] + _pi[tip] < 0); // admissible
                assert(_res_cap[edge] > 0);
                assert(_res_cap[_reverse[edge]] == 0);
                _res_cap[edge]--;
                _res_cap[_reverse[edge]]++;
                _excess[parent_node]--;
                _excess[tip]++;
                tip = parent_node;
            }

//            assert(_excess[tip] >= 0);
//            if (_excess[tip] == 0) {
//                deque<int>::iterator it = std::find(_active_nodes.begin(),_active_nodes.end(),tip);
//                _active_nodes.erase(it);
//            }

            _active_nodes.clear();
            for (int u = 0; u != _res_node_num; ++u) {
                if (_excess[u] > 0) _active_nodes.push_back(u);
            }
        }

        if (_active_nodes.size() > 0) assert(test_epsilon_optimality());

    }
}
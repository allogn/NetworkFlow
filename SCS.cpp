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
            vector<pair<int, LargeCost>> visited;
            vector<int> parent(_res_node_num, -1); //parent node and an edge to a child
            LargeCost mindist[_res_node_num];
            if (QryCnt[start] < _graph.fullE[start].size())
                globalH.enqueue(start, 0);

            LargeCost curcost;
            int u, tip;
            do {
                dijkH.clear();
                dijkH.enqueue(start,0);
                visited.clear();
                std::fill(parent.begin(),parent.end(),-1);
                fill(mindist, mindist + _res_node_num,
                     std::numeric_limits<LargeCost>::max()); //those who were not popped yet, but are in the bucket already
                mindist[_active_nodes.front()] = 0;

                int fromid;
                LargeCost tmp;
                globalH.dequeue(fromid,tmp);
                int eid = _graph.fullE[fromid][QryCnt[fromid]];
                QryCnt[fromid]++;
                if (QryCnt[fromid] < _graph.fullE[fromid].size()) {
                    globalH.enqueue(fromid, mindist[fromid] + _graph.E[_graph.fullE[fromid][QryCnt[fromid]]].weight);
                }

                uintT local_eid = _res_arc_num;
                //for spatial data this does not work
                _first_out[fromid].push_back(local_eid);
                int toid = _graph.get_pair(eid, fromid);
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
//                LargeCost new_epsilon = std::max(-(lc + _pi[toid] - _pi[fromid]), _epsilon);
//                _epsilon = new_epsilon; //todo maybe this can be modified
                _res_cap.push_back(_upper[local_eid]);
                _res_cap.push_back(0);
                _reverse.push_back(local_eid + 1);
                _reverse.push_back(local_eid);

                _res_arc_num += 2;

                while (true) {
                    dijkH.dequeue(tip, curcost);
                    //take one node and add every neighbor to another bucket
                    if (_excess[tip] < 0) goto checkEdges;
                    if (!globalH.isExisted(tip) && QryCnt[tip] < _graph.fullE[tip].size()) {
                        globalH.enqueue(tip, mindist[tip] + _graph.E[_graph.fullE[tip][QryCnt[tip]]].weight);
                    }
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
                            }
                            parent[u] = *a;
                            mindist[u] = curcost + l;
                        }
                    }
                }
                checkEdges:;
            } while (globalH.size() > 0);


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
                deque<int>::iterator it = std::find(_active_nodes.begin(), _active_nodes.end(), tip);
                _active_nodes.erase(it);
            }
        }

    }
}
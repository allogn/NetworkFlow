//
// Created by alvis on 20.03.16.
//

#ifndef NETWORKFLOW_SCS_H
#define NETWORKFLOW_SCS_H

/*
 * Simplified Cost Scaling
 *
 * Inherits code from Lemon with embedding of SIA graph structure
 */

#include <vector>
#include <deque>
#include <queue>
#include <limits>
#include <stack>
#include <map>

#include "Common.h"
#include "Graph.h"
#include "TimerTool.h"

class SCS {

    /*
     * Working with integer cost type
     * and long long Large
     */

    typedef long long LargeCost;
    typedef int Value;
    typedef int Cost;
    //TODO change this in other algorithms

private:


    typedef std::vector<int> IntVector;
    typedef std::vector <Value> ValueVector;
    typedef std::vector <Cost> CostVector;
    typedef std::vector <LargeCost> LargeCostVector;
    typedef std::vector<char> BoolVector;
    // Note: vector<char> is used instead of vector<bool>
    // for efficiency reasons

public:

    Timer timer;

    // Data related to the underlying digraph
    const Graph &_graph;
    int _node_num;
    int _arc_num;
    int _res_node_num;
    int _res_arc_num;

    // Parameters of the problem
    bool _has_lower;
    Value _sum_supply;

    // Data structures for storing the digraph
    IntVector _node_id;
    IntVector _arc_idf;
    IntVector _arc_idb;
    vector<vector<int>> _first_out;
    BoolVector _forward;
    IntVector _source;
    IntVector _target;
    IntVector _reverse;

    // Node and arc data
    ValueVector _lower;
    ValueVector _upper;
    CostVector _scost;
    ValueVector _supply;

    ValueVector _res_cap;
    LargeCostVector _cost;
    LargeCostVector _pi;
    ValueVector _excess;
    IntVector _next_out;
    std::deque<int> _active_nodes;

    // Data for scaling
    LargeCost _epsilon;
    int _alpha;

//    mmHeap dijkH;

public:

    SCS(const Graph &graph) :
            _graph(graph) {

        // Reset data structures
        reset();
    }

    SCS &reset() {
        // Resize vectors
        _node_num = _graph.n;
        _arc_num = _graph.m;
        _res_node_num = _node_num;
        _res_arc_num = 2 * _arc_num;

        _first_out.resize(_res_node_num, vector<int>());
        _forward.resize(_res_arc_num);
        _source.resize(_res_arc_num);
        _target.resize(_res_arc_num);
        _reverse.resize(_res_arc_num);

        _lower.resize(_res_arc_num);
        _upper.resize(_res_arc_num);
        _scost.resize(_res_arc_num);
        _supply.resize(_res_node_num);

        _res_cap.resize(_res_arc_num);
        _cost.resize(_res_arc_num);
        _pi.resize(_res_node_num);
        _excess.resize(_res_node_num);
        _next_out.resize(_res_node_num);

        _arc_idf.resize(_res_arc_num);
        _arc_idb.resize(_res_arc_num);

        // Copy the graph
        int j = 0;
        for (uintT i = 0; i < _graph.n; i++) {
            _supply[i] = _graph.V[i].supply;
        }
        for (uintT i = 0; i < _graph.n; i++) {
            for (uintT a = 0; a < _graph.completeE[i].size(); a++) {
                _first_out[i].push_back(j);
                uintT eid = _graph.completeE[i][a];
                if (_graph.is_forward(eid,i)) {
                    _arc_idf[eid] = j;
                } else {
                    _arc_idb[eid] = j;
                }
                _forward[j] = _graph.is_forward(eid,i);
                _source[j] = i;
                _target[j] = _graph.get_pair(eid, i);

                _lower[j] = _graph.E[eid].lower;
                _upper[j] = _graph.E[eid].capacity;
                _scost[j] = _graph.E[eid].weight;
                if (!_forward[j]) _scost[j] *= -1;

                if (_forward[j]) _res_cap[j] = _upper[j];
                _cost[j] = _scost[j] * _res_node_num * _alpha;

                j++;
            }
        }
        _has_lower = false;

        for (int eid = 0; eid < _graph.m; eid++) {
            int fi = _arc_idf[eid];
            int bi = _arc_idb[eid];
            _reverse[fi] = bi;
            _reverse[bi] = fi;
        }

        return *this;
    }

private:



    // Initialize the algorithm
    ProblemType init() {
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
//    if (_has_lower) {
//        for (int i = 0; i != _res_node_num; ++i) {
//            last_out = (i < _res_node_num-1)?_first_out[i+1]:_res_arc_num;
//            for (int j = _first_out[i]; j != last_out; ++j) {
//                if (_forward[j]) {
//                    Value c = _scost[j] < 0 ? _upper[j] : _lower[j];
//                    if (c >= MAX) return UNBOUNDED;
//                    _excess[i] -= c;
//                    _excess[_target[j]] += c;
//                }
//            }
//        }
//    } else {
//        for (int i = 0; i != _res_node_num; ++i) {
//            last_out = (i < _res_node_num-1)?_first_out[i+1]:_res_arc_num;
//            for (int j = _first_out[i]; j != last_out; ++j) {
//                if (_forward[j] && _scost[j] < 0) {
//                    Value c = _upper[j];
//                    if (c >= MAX) return UNBOUNDED;
//                    _excess[i] -= c;
//                    _excess[_target[j]] += c;
//                }
//            }
//        }
//    }
//            Value ex, max_cap = 0;
//            for (int i = 0; i != _res_node_num; ++i) {
//                ex = _excess[i];
//                _excess[i] = 0;
//                if (ex < 0) max_cap -= ex;
//            }
//            for (int j = 0; j != _res_arc_num; ++j) {
//                if (_upper[j] >= MAX) _upper[j] = max_cap;
//            }

        // Initialize the large cost vector and the epsilon parameter
//        _epsilon = 0;
        LargeCost lc;
        for (int i = 0; i != _res_node_num; ++i) {
            for (vector<int>::iterator j = _first_out[i].begin(); j != _first_out[i].end(); ++j) {
                lc = static_cast<LargeCost>(_scost[*j]) * _res_node_num * _alpha; //COST MODIFICATION
                _cost[*j] = lc;
                if (lc > _epsilon) _epsilon = lc;
            }
        }
        _epsilon /= _alpha;

//        LargeCost lc;
//        for (int i = 0; i != _graph.E.size(); ++i) {
//            lc = static_cast<LargeCost>(_graph.E[i].weight) * _res_node_num * _alpha; //COST MODIFICATION
//            if (lc > _epsilon) _epsilon = lc;
//        }
//        _epsilon /= _alpha;

        // initialize _res_cap with supply value for each node with positive supply for arbitrary edge
        for (int a = 0; a < _res_arc_num; a++) {
            if (_forward[a]) {
                _res_cap[a] = _upper[a];
            } else {
                _res_cap[a] = 0;
            }
        }


        // init spatial structures
        QryCnt.resize(_res_node_num,0);
        globalH.clear();
        taumax = 0;

        return OPTIMAL;
    }

    // Check if the upper bound is greater than or equal to the lower bound
    // on each forward arc.
    bool checkBoundMaps() {
        for (int j = 0; j != _res_arc_num; ++j) {
            if (_forward[j] && _upper[j] < _lower[j]) return false;
        }
        return true;
    }

    // Initialize a cost scaling phase
    void initPhase() {
        // Saturate arcs not satisfying the optimality condition
//        for (int u = 0; u != _res_node_num; ++u) {
//            LargeCost pi_u = _pi[u];
//            for (vector<int>::iterator a = _first_out[u].begin(); a != _first_out[u].end(); ++a) {
//                Value delta = _res_cap[*a];
//                if (delta > 0) {
//                    int v = _target[*a];
//                    if (_cost[*a] - pi_u + _pi[v] < 0) {
//                        _excess[u] -= delta;
//                        _excess[v] += delta;
//                        _res_cap[*a] = 0;
//                        _res_cap[_reverse[*a]] += delta;
//                    }
//                }
//            }
//        }

        // Find active nodes (i.e. nodes with positive excess)
        _active_nodes.clear();
        for (int u = 0; u != _res_node_num; ++u) {
            if (_excess[u] > 0) _active_nodes.push_back(u);
        }

        // Initialize the next arcs
//        for (int u = 0; u != _res_node_num; ++u) {
//            _next_out[u] = _first_out[u];
//        }
    }

    /// Execute the algorithm performing augment and relabel operations
    void startAugment(int max_length);

public:
    uintT totalCost;

    // Execute the algorithm and transform the results
    void runSCS() {
        _alpha = 16;
        reset();
        init();
        double total = timer.getTime();
        startAugment(_res_node_num - 1);
        timer.save_time("Total time", total);
        totalCost = 0;
        for (uintT n = 0; n < _res_node_num; n++) {
            for (vector<int>::iterator it = _first_out[n].begin(); it != _first_out[n].end(); it++ ) {
                if (!_forward[*it]) {
                    assert(_scost[*it] <= 0);
                    totalCost += _res_cap[*it] * (-_scost[*it]);
                }
            }
        }
    }


    /*
     * from SIA
     */
    vector<int> QryCnt;
    vector<LargeCost> mindist;
    vector<int> visited;
    vector<bool> watched; //traversed all children
    vector<int> parent; //parent node and an edge to a child
    LargeCost taumax;

    fHeap<LargeCost> dijkH;
    fHeap<LargeCost> globalH;
    fHeap<LargeCost> updateH;
    inline int heap_checkAndUpdateEdgeMin(fHeap<LargeCost>& heap, int fromid) // update if new value is less
    {
        //todo add sorting of fullE

        if (QryCnt[fromid] >= _graph.fullE[fromid].size()) return 0;
        //if fromid is in the heap - update distance to the nearest not added node since mindist to fromid may change
        if (heap.isExisted(fromid))
        {
            int weight = _graph.E[_graph.fullE[fromid][QryCnt[fromid]]].weight;
            heap.updatequeue(fromid, weight + mindist[fromid]);
            return 1;
        }
        heap.enqueue(fromid, _graph.E[_graph.fullE[fromid][QryCnt[fromid]]].weight + mindist[fromid]);
        return 0;
    }
    int heap_checkAndUpdateMin(fHeap<LargeCost>& heap, int id, int new_value) // update if new value is less
    {
        if (heap.isExisted(id))
        {
            heap.updatequeue(id,new_value);
            return 1;
        }
        else
            heap.enqueue(id,new_value);
        return 0;
    }
    inline int updateMinDist(int local_eid, int fromid, int toid)
    {
        assert(_target[local_eid] == toid && _source[local_eid] == fromid);
        long cost = _cost[local_eid] - _pi[fromid] + _pi[toid] + _epsilon;
        if (mindist[toid]>mindist[fromid]+cost) {
            mindist[toid] = mindist[fromid] + cost;
            parent[toid] = local_eid;
            return 1;
        }
        return 0;
    }
    void updateHeaps(int local_eid, int fromid, int toid)
    {
        assert(_target[local_eid] == toid && _source[local_eid] == fromid);
        if (watched[fromid] == 0) return; //case when prefinal node : not all neighbours are considered => not watched,
        if (_res_cap[local_eid] == 0) return;
        // but in globalH and in DijkH (!). so, if in dijkH => everything is fine (will be updated later)
        if (updateMinDist(local_eid,fromid,toid)) {
            //no isUpdated because enheap only if updated dist
            if (!dijkH.isExisted(toid)) {
                //isUpdated omitted here
                if (watched[toid] == 1) heap_checkAndUpdateMin(updateH, toid, mindist[toid]);
                else dijkH.enqueue(toid, mindist[toid]);
            } else {
                dijkH.updatequeue(toid, mindist[toid]);
            }
        }
    }

    bool test_epsilon_optimality() {
        for (int i = 0; i < _res_node_num; i++) {
            for (vector<int>::iterator it = _first_out[i].begin(); it != _first_out[i].end(); ++it) {
                if (_res_cap[*it] > 0) {
                    assert(_cost[*it] + _pi[_source[*it]] - _pi[_target[*it]] > - _epsilon);
                }
            }
        }
        return true;
    }

    bool test_if_shortest(int source, int target) {
        LargeCost mindist[_res_node_num];
        vector<int> sourceEdge(_res_node_num, 0);
        vector<int> localParent(_res_node_num, -1);
        fill(mindist, mindist+_res_node_num, std::numeric_limits<LargeCost>::max());
        mindist[source] = 0;
        queue<int> q;
        q.push(source);
        while (q.size() > 0) {
            int curnode = q.front();
            q.pop();
            for (int i = 0; i < _first_out[curnode].size(); i++) {
                int nb = _target[_first_out[curnode][i]];
                if (_res_cap[_first_out[curnode][i]] == 0) continue;
                LargeCost cp = _cost[_first_out[curnode][i]] + _pi[nb] - _pi[curnode];
                LargeCost newdist = mindist[curnode] + cp + _epsilon;
                assert(newdist >= 0);
                if (mindist[nb] > newdist) {
                    q.push(nb);
                    mindist[nb] = newdist;
                    sourceEdge[nb] = 1;
                    localParent[nb] = curnode;
                }
            }
        }
        assert(_excess[target] < 0);
        assert(mindist[target] < std::numeric_limits<LargeCost>::max());
        if (target == -1) {
            assert(mindist[target] == std::numeric_limits<LargeCost>::max());
            return true;
        }
        for (int i = 0; i < _res_node_num; i++) {
            if (_excess[i] < 0) {
                assert(mindist[target] <= mindist[i]);
            }
        }
        return true;
    }

    bool test_dijkstra_heap(int source) {
        LargeCost mindist[_res_node_num];
        vector<int> localParent(_res_node_num, -1);
        fill(mindist, mindist+_res_node_num, std::numeric_limits<LargeCost>::max());
        mindist[source] = 0;
        queue<int> q;
        q.push(source);
        while (q.size() > 0) {
            int curnode = q.front();
            q.pop();
            for (int i = 0; i < _first_out[curnode].size(); i++) {
                int nb = _target[_first_out[curnode][i]];
                if (_res_cap[_first_out[curnode][i]] == 0) continue;
                if (dijkH.isExisted(curnode)) continue;
                LargeCost cp = _cost[_first_out[curnode][i]] + _pi[nb] - _pi[curnode];
                LargeCost newdist = mindist[curnode] + cp + _epsilon;
                assert(newdist >= 0);
                if (mindist[nb] > newdist) {
                    q.push(nb);
                    mindist[nb] = newdist;
                    localParent[nb] = curnode;
                }
            }
        }
        /*
         * For all nodes that were watched + all in dijkstra heap - distances must be equal
         * + all paths in the graph must lead to nodes in dijkstra heap
         */
        stack<int> queue;
        queue.push(source);
        vector<bool> visited(_res_node_num, false);
        while(!queue.empty()) {
            int curnode = queue.top();
            visited[curnode] = true;
            queue.pop();
            if (!dijkH.isExisted(curnode)) {
                assert(mindist[curnode] == this->mindist[curnode]); //if two nodes are in dijkH and there are path
                // from one node to another - then distance to one of them not guaranteed to be minimum
                // (equal to dial's solution) because in dial's that edge may not have been considered,
                // but it will be considered here
                for (int j = 0; j < _first_out[curnode].size(); j++) {
                    if (_res_cap[_first_out[curnode][j]] > 0 && visited[_target[_first_out[curnode][j]]] == false)
                        queue.push(_target[_first_out[curnode][j]]);
                }
            }
        }
    }

};


#endif //NETWORKFLOW_SCS_H

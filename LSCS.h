//
// Created by alvis on 20.03.16.
//

#ifndef NETWORKFLOW_LSCS_H
#define NETWORKFLOW_LSCS_H

/*
 * Lemon Simplified Cost Scaling
 *
 * Inherits code from Lemon with embedding of SIA graph structure
 */

#include <vector>
#include <deque>
#include <limits>
#include <map>

#include "Common.h"
#include "Graph.h"
#include "TimerTool.h"

class LSCS {

    /*
     * Working with integer cost type
     * and long long Large
     */

    typedef long long LargeCost;
    typedef int Value;
    typedef int Cost;
    //TODO change this in other algorithms

private:


    typedef std::vector<long> IntVector;
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
    long _node_num;
    long _res_node_num;
    long _res_arc_num;

    // Parameters of the problem
    bool _has_lower;
    Value _sum_supply;

    // Data structures for storing the digraph
    IntVector _arc_idf;
    IntVector _arc_idb;
    vector<vector<long>> _first_out;
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
    std::deque<long> _active_nodes;

    // Data for scaling
    LargeCost _epsilon;
    int _alpha;

    vector<int> QryCnt;

    //profiling
    vector<bool> _ever_visited;

public:

    LSCS(const Graph &graph) :
            _graph(graph) {

        // Reset data structures
        reset();
    }

    LSCS &reset() {
        // Resize vectors
        _node_num = _graph.n;
        _res_node_num = _node_num;
        _res_arc_num = 0;//2 * _arc_num;

        _first_out.resize(_res_node_num, vector<long>());
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

        QryCnt.resize(_res_node_num,0);

        _ever_visited.resize(_res_node_num,false);

        // Copy the graph
//        int j = 0;
        for (long i = 0; i < _graph.n; i++) {
            _supply[i] = _graph.V[i].supply;
        }
//        for (long i = 0; i < _graph.n; i++) {
//            for (long a = 0; a < _graph.completeE[i].size(); a++) {
//                _first_out[i].push_back(j);
//                long eid = _graph.completeE[i][a];
//                if (_graph.is_forward(eid,i)) {
//                    _arc_idf[eid] = j;
//                } else {
//                    _arc_idb[eid] = j;
//                }
//                _forward[j] = _graph.is_forward(eid,i);
//                _source[j] = i;
//                _target[j] = _graph.get_pair(eid, i);
//
//                _lower[j] = _graph.E[eid].lower;
//                _upper[j] = _graph.E[eid].capacity;
//                _scost[j] = _graph.E[eid].weight;
//                if (!_forward[j]) _scost[j] *= -1;
//
//                if (_forward[j]) _res_cap[j] = _upper[j];
//
//                j++;
//            }
//        }
        _has_lower = false;

//        for (int eid = 0; eid < _graph.m; eid++) {
//            int fi = _arc_idf[eid];
//            int bi = _arc_idb[eid];
//            _reverse[fi] = bi;
//            _reverse[bi] = fi;
//        }

        return *this;
    }

private:



    // Initialize the algorithm
    ProblemType init() {
        if (_res_node_num <= 1) return INFEASIBLE;

        // Check the sum of supply values
        _sum_supply = 0;
        for (long i = 0; i != _res_node_num; ++i) {
            _sum_supply += _supply[i];
        }
        if (_sum_supply > 0) return INFEASIBLE;

        // Check lower and upper bounds
        LEMON_DEBUG(checkBoundMaps(),
                    "Upper bounds must be greater or equal to the lower bounds");


        // Initialize vectors
        for (long i = 0; i != _res_node_num; ++i) {
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
        _epsilon = 0;
        LargeCost lc;
        for (long i = 0; i != _res_node_num; ++i) {
//            for (vector<int>::iterator j = _first_out[i].begin(); j != _first_out[i].end(); ++j) {
//                lc = static_cast<LargeCost>(_scost[*j]) * _res_node_num * _alpha; //COST MODIFICATION
//                _cost[*j] = lc;
//                if (lc > _epsilon) _epsilon = lc;
//            }
            for (long j = 0; j < _graph.fullE[i].size(); j++) {
                lc = static_cast<LargeCost>(_graph.E[_graph.fullE[i][j]].weight) * _res_node_num * _alpha;
                if (lc > _epsilon) _epsilon = lc;
            }
        }
        _epsilon /= _alpha;

        // initialize _res_cap with supply value for each node with positive supply for arbitrary edge
        for (long a = 0; a < _res_arc_num; a++) {
            if (_forward[a]) {
                _res_cap[a] = _upper[a];
            } else {
                _res_cap[a] = 0;
            }
        }

        return OPTIMAL;
    }

    // Check if the upper bound is greater than or equal to the lower bound
    // on each forward arc.
    bool checkBoundMaps() {
        for (long j = 0; j != _res_arc_num; ++j) {
            if (_forward[j] && _upper[j] < _lower[j]) return false;
        }
        return true;
    }

    // Initialize a cost scaling phase
    void initPhase() {
        // Saturate arcs not satisfying the optimality condition
        for (int u = 0; u != _res_node_num; ++u) {
            LargeCost pi_u = _pi[u];
            for (vector<long>::iterator a = _first_out[u].begin(); a != _first_out[u].end(); ++a) {
                Value delta = _res_cap[*a];
                if (delta > 0) {
                    long v = _target[*a];
                    if (_cost[*a] + pi_u - _pi[v] < 0) {
                        _excess[u] -= delta;
                        _excess[v] += delta;
                        _res_cap[*a] = 0;
                        _res_cap[_reverse[*a]] += delta;
                    }
                }
            }
        }

        // Find active nodes (i.e. nodes with positive excess)
        for (long u = 0; u != _res_node_num; ++u) {
            if (_excess[u] > 0) _active_nodes.push_back(u);
        }

        // Initialize the next arcs
        for (long u = 0; u != _res_node_num; ++u) {
            _next_out[u] = 0;
        }
    }

    /// Execute the algorithm performing augment and relabel operations
    void startAugment(long max_length);
    void addArc(long fromid);

public:
    long totalCost;

    // Execute the algorithm and transform the results
    void runLSCS() {
        _alpha = 16;
        reset();
        init();
        double total = timer.getTime();
        startAugment(_res_node_num - 1);
        timer.save_time("Total time", total);
        totalCost = 0;
        for (long a = 0; a < _target.size(); a++) {
            if (!_forward[a]) {
                assert(_scost[a] <= 0);
                totalCost += _res_cap[a] * (-_scost[a]);
            }
        }
    }

    void save_profile_data(string filename, int experiment_id) {
        ofstream outf(filename,ios::app);

        //mean and std of how much edges were added
        double perc = 0;
        double percSqr = 0;
        double totalNodes = 0;
        for (long i = 0; i < _res_node_num; i++) {
            long added = 0;
            for (long j = 0; j < _first_out[i].size(); j++) {
                if (_forward[_first_out[i][j]]) {
                    added++;
                }
            }
            double val = (double)added/(double)_graph.fullE[i].size();
            perc += val;
            assert(val >= 0 && val <= 1);
            percSqr += val*val;
        }
        long visited = std::count_if(_ever_visited.begin(), _ever_visited.end(), [](bool v) { return v; });
        outf << experiment_id << ",Visited fraction," << (double)visited/(double)_res_node_num << "\n";
        outf << experiment_id << ",Saturation mean," << perc/(double)_res_node_num << "\n";
        outf << experiment_id << ",Saturation std," << sqrt(percSqr/(double)_res_node_num - (perc/(double)_res_node_num)*(perc/(double)_res_node_num)) << "\n";

        outf.close();
    }

};


#endif //NETWORKFLOW_LSCS_H

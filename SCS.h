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
#include <limits>

#include "Common.h"
#include "Graph.h"

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

private:

    // Data related to the underlying digraph
    const Graph &_graph;
    int _node_num;
    int _arc_num;
    int _res_node_num;
    int _res_arc_num;

    // Parameters of the problem
    bool _has_lower;
    Value _sum_supply;
    int _sup_node_num;

    // Data structures for storing the digraph
    IntVector _node_id;
    IntVector _arc_idf;
    IntVector _arc_idb;
    IntVector _first_out;
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

    IntVector _buckets;
    IntVector _bucket_next;
    IntVector _bucket_prev;
    IntVector _rank;
    int _max_rank;

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
        _res_arc_num = 2 * (_arc_num + _node_num);

        _first_out.resize(_res_node_num);
        _forward.resize(_res_arc_num);
        _source.resize(_res_arc_num);
        _target.resize(_res_arc_num);
        _reverse.resize(_res_arc_num);
        _node_id.resize(_res_arc_num);

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
        int j = 0, k = 2 * _arc_num + _node_num;
        for (uintT i = 0; i < _graph.n; i++) {
            _node_id[i] = i;
            _supply[i] = _graph.V[i].supply;
        }
        for (uintT i = 0; i < _graph.n; i++) {
            _first_out[i] = j;
            for (uintT a = 0; a < _graph.completeE[i].size(); a++) {
                uintT eid = _graph.completeE[i][a];
                if (_graph.is_forward(eid,i)) {
                    _arc_idf[eid] = j;
                } else {
                    _arc_idb[eid] = j;
                }
                _forward[j] = _graph.is_forward(eid,i);
                _source[j] = i;
                _target[j] = _node_id[_graph.get_pair(eid, i)];

                _lower[j] = _graph.E[eid].lower;
                _upper[j] = _graph.E[eid].capacity;
                _scost[j] = _graph.E[eid].weight;
                if (!_forward[j]) _scost[j] *= -1;

                if (_forward[j]) _res_cap[j] = _upper[j];

                j++;
            }
        }
        _has_lower = true;

        return *this;
    }

private:



    // Initialize the algorithm
    ProblemType init();

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
        for (int u = 0; u != _res_node_num; ++u) {
            int last_out =(u<_res_node_num-1)?_first_out[u + 1]:_res_arc_num;
            LargeCost pi_u = _pi[u];
            for (int a = _first_out[u]; a != last_out; ++a) {
                Value delta = _res_cap[a];
                if (delta > 0) {
                    int v = _target[a];
                    if (_cost[a] + pi_u - _pi[v] < 0) {
                        _excess[u] -= delta;
                        _excess[v] += delta;
                        _res_cap[a] = 0;
                        _res_cap[_reverse[a]] += delta;
                    }
                }
            }
        }

        // Find active nodes (i.e. nodes with positive excess)
        for (int u = 0; u != _res_node_num; ++u) {
            if (_excess[u] > 0) _active_nodes.push_back(u);
        }

        // Initialize the next arcs
        for (int u = 0; u != _res_node_num; ++u) {
            _next_out[u] = _first_out[u];
        }
    }

    /// Execute the algorithm performing augment and relabel operations
    void startAugment(int max_length);

public:
    uintT totalCost;

    // Execute the algorithm and transform the results
    void runSCS() {
        init();

        startAugment(_res_node_num - 1);

        totalCost = 0;
        for (uintT a = 0; a < _graph.m; a++) {
            int i = _arc_idb[a];
            totalCost += _res_cap[i] * _scost[i];
        }
    }

};


#endif //NETWORKFLOW_SCS_H

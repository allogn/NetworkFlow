//
// Created by alvis on 20.03.16.
//

#ifndef NETWORKFLOW_SCS_H
#define NETWORKFLOW_SCS_H

#include <vector>
#include <deque>
#include <limits>

#include "Graph.h"
#include "TimerTool.h"

class SCS {
public:
    typedef int Value;
    typedef int Cost;
    typedef long LargeCost;
    typedef std::vector<int> IntVector;
    typedef std::vector<Value> ValueVector;
    typedef std::vector<LargeCost> LargeCostVector;
    typedef std::vector<char> BoolVector;

    // Data related to the underlying digraph
    Graph &_graph;
    int _res_node_num;
    int _res_arc_num;

    LargeCostVector _cost;
    LargeCostVector _pi;
    ValueVector _excess;
    IntVector _next_out;
    std::deque<int> _active_nodes;

    long totalCost = -1;

    // Data for scaling
    LargeCost _epsilon;
    int _alpha;

    Timer timer;

    SCS(Graph& g) :
            _graph(g) {
    }

    ProblemType runSCS(int factor = 16) {
        _alpha = factor;
        init();
        double total = timer.getTime();
        startAugment(_res_node_num - 1);
        timer.save_time("Total time", total);
        totalCost = 0;
        for (int i = 0; i < _res_node_num; i++) {
            for (int j = 0; j < _graph.V[i].E.size(); j++) {
                int eid = _graph.V[i].E[j];
                if (!_graph.is_forward(eid, i)) {
                    totalCost += _graph.E[eid].weight; //unit capacities implied
                }
            }
        }
        return OPTIMAL;
    }

    inline long get_cp(int nodeid, int edge) {
        int eid = _graph.V[nodeid].E[edge];
        int nb = _graph.get_pair(eid, nodeid);
        long cost = _cost[eid];//_graph.E[eid].weight;
        if (!_graph.is_forward(eid, nodeid))
            cost *= -1;
        cost = cost - _pi[nb] + _pi[nodeid]; // - Delta(Pi)
        return cost;
    };

    inline void raiseFlow(uintT nodeid, int edge) {
        //assuming unit capacity
        uintT eid = _graph.V[nodeid].E[edge];
        int nb = _graph.get_pair(eid, nodeid);
        //add to nb E list
        _graph.V[nb].E.push_back(eid);
        //delete from nodeid
        _graph.V[nodeid].E[edge] = *_graph.V[nodeid].E.rbegin();
        _graph.V[nodeid].E.pop_back();
        //change excesses
        _excess[nodeid]--;
        _excess[nb]++;
    }

    ProblemType init() {
        _res_node_num = _graph.n;
        _res_arc_num = _graph.m;
        _cost.resize(_res_arc_num);
        _pi.resize(_res_node_num);
        _excess.resize(_res_node_num);
        _next_out.resize(_res_node_num);

        if (_res_node_num <= 1) return INFEASIBLE;

        // Initialize vectors
        for (int i = 0; i != _res_node_num; ++i) {
            _pi[i] = 0;
            _excess[i] = _graph.V[i].supply;
        }

        // Initialize the large cost vector and the epsilon parameter
        _epsilon = 0;
        LargeCost lc;
        for (long i = 0; i < _res_arc_num; i++) {
            lc = static_cast<LargeCost>(_graph.E[i].weight) * _res_node_num * _alpha;
            _cost[i] = lc;
            if (lc > _epsilon) _epsilon = lc;
        }
        _epsilon /= _alpha;

        return OPTIMAL;
    }

    // Initialize a cost scaling phase
    void initPhase() {
        // Saturate arcs not satisfying the optimality condition
        for (int u = 0; u < _res_node_num; u++) {
            for (int e = 0; e < _graph.V[u].E.size(); e++) {
                //if this edge is in node's E then it has an excess
                //check admissibility
                if (get_cp(u, e) < 0) {
                    raiseFlow(u, e);
                }
            }
        }

        // Find active nodes (i.e. nodes with positive excess)
        for (int u = 0; u != _res_node_num; ++u) {
            if (_excess[u] > 0) _active_nodes.push_back(u);
        }

        // Initialize the next arcs
        for (int u = 0; u != _res_node_num; ++u) {
            _next_out[u] = 0;
        }
    }

    /// Execute the algorithm performing augment and relabel operations
    void startAugment(int max_length) {

        // Perform cost scaling phases
        IntVector path;
        IntVector node_path; //for back propagation
        BoolVector path_arc(_res_arc_num, false);
        int relabel_cnt = 0;
        int eps_phase_cnt = 0;
        for (; _epsilon >= 1; _epsilon = _epsilon < _alpha && _epsilon > 1 ?
                                         1 : _epsilon / _alpha) {
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
                int start = _active_nodes.front();

                // Find an augmenting path from the start node
                int tip = start;
                int last_pushed;
                while (int(path.size()) < max_length && _excess[tip] >= 0) {
                    LargeCost rc, min_red_cost = std::numeric_limits<LargeCost>::max();
                    LargeCost pi_tip = _pi[tip];

                    for (int a = _next_out[tip]; a < _graph.V[tip].E.size(); a++) {
                        long rc = get_cp(tip, a);
                        if (rc < 0) {
                            path.push_back(a);
                            node_path.push_back(tip);
                            last_pushed = tip;
                            _next_out[tip] = a;
                            long eid = _graph.V[tip].E[a];
                            if (path_arc[eid]) {
                                goto augment;   // a cycle is found, stop path search
                            }
                            tip = _graph.get_pair(eid,tip);
                            path_arc[eid] = true;
                            goto next_step;
                        }
                        else if (rc < min_red_cost) {
                            min_red_cost = rc;
                        }
                    }

                    // Relabel tip node
//                    assert(path.size() == 0 ||
//                                   _graph.get_pair(_graph.V[last_pushed].E[path.back()],last_pushed) == tip); //TODO removethis if true or fix if false
//                    if (tip != start) {
//                        int ra = _reverse[path.back()];
//                        min_red_cost =
//                                std::min(min_red_cost, _cost[ra] + pi_tip - _pi[_target[ra]]);
//                    }
                    for (int a = 0; a < _graph.V[tip].E.size(); a++) {
                        long cp = get_cp(tip, a);
                        if (cp < min_red_cost) {
                            min_red_cost = cp;
                        }
                    }
                    _pi[tip] -= min_red_cost + _epsilon;
                    _next_out[tip] = 0;
                    ++relabel_cnt;

                    // Step back
                    if (tip != start) {
                        int pa = path.back(); // last edge that lead from node_path.back() to tip
                        int nb = node_path.back();
                        int eid = _graph.V[nb].E[pa];
                        path_arc[eid] = false;
                        tip = nb;
                        path.pop_back();
                        node_path.pop_back();
                    }

                    next_step:;
                }
                // Augment along the found path (as much flow as possible)
                augment:
                Value delta;
                int eid, pa, u, v = start;
                for (int i = 0; i != int(path.size()); ++i) {
                    pa = path[i];
                    u = v;
                    eid = _graph.V[u].E[pa];
                    v = _graph.get_pair(eid, u);
                    path_arc[eid] = false;

                    raiseFlow(u, pa);
                    if (_excess[v] > 0) {
                        _active_nodes.push_back(v);//add target node if needed
                    }
                }
                path.clear();
                node_path.clear();
            }

        }

    }
};

#endif //NETWORKFLOW_SCS_H

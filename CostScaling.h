//
// Created by alvis on 18.03.16.
//

#ifndef NETWORKFLOW_COSTSCALING_H
#define NETWORKFLOW_COSTSCALING_H

#include <iostream>
#include <stack>
#include <queue>

#include "Graph.h"
#include "Common.h"
#include "TimerTool.h"
#include "NodeList.h"

class CostScaling {
    Graph* g;
    long* excesses; // stores flow amount for each node
    long* flow; // flow amount for each edge
    bool* admissible;
    double epsilon; // initialization in reset()
    double* p;
    long totalExcesses; // used in process: total nodes with sum of flows
    Timer timer;

    inline double get_cp(long eid, long nodeId) {
        long neighbor = g->get_pair(eid, nodeId);
        if (g->is_forward(eid, nodeId)) {
            return (double)g->E[eid].weight + p[neighbor] - p[nodeId];
        } else {
            return -(double)g->E[eid].weight + p[neighbor] - p[nodeId];
        }
    }
    inline long get_residual(long eid, long nodeId) {
        if (g->is_forward(eid, nodeId)) {
            return g->E[eid].capacity - flow[eid];
        } else {
            return flow[eid] - g->E[eid].lower;
        }
    }

    inline bool isAdmissible(long node, long outEdge) {
        long eid = g->V[node].E[outEdge];
        long neighbor = g->get_pair(eid, node);
        double c_p = get_cp(eid, node);
        long residual = get_residual(eid, node);
        if (c_p < 0 && residual > 0) {
            admissible[eid] = true;//TODO check admissible where is used and set
        }
        return (c_p < 0 && residual > 0);
    }
    void changeFlow(long nodeid, long edgeInNodeId, long value);
    void dfs(long* blockingFlow, long* blockingEdge, long nodeId);
    void raise_potentials();
    void raise_flows();
    void refine();

    /*
     * Edges are assigned to both nodes if possible to change flow in both directions
     * Decision if an edge is forward/backward is done based on the fromid and toid values
     */

public:
    long totalCost;

    CostScaling(Graph* graph) {
        g = graph;
        totalCost = 0;
        excesses = new long[g->n];
        flow = new long[g->m];
        admissible = new bool[g->m];
        p = new double[g->n];
        reset();
    }
    ~CostScaling() {
        delete[] excesses;
        delete[] flow;
        delete[] admissible;
        delete[] p;
    }
    void reset();
    void runCostScaling();

    // Unit Tests
    bool is_feasible(double threshold);
    bool is_flow();
    bool test_excesses();
};


#endif //NETWORKFLOW_COSTSCALING_H

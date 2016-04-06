//
// Created by alvis on 19.03.16.
//

#ifndef NETWORKFLOW_LOCALDOMINANT_H
#define NETWORKFLOW_LOCALDOMINANT_H

#include <iostream>
#include <queue>

#include "Common.h"
#include "Graph.h"
#include "TimerTool.h"

struct AssignE {
    int from;
    int to;
    int cost;
    int cap;
    AssignE (int _from, int _to, int _cost, int _cap): from(_from), to(_to), cost(_cost), cap(_cap) {}
};

class LocalDominant {
    Graph* g;

    std::queue<int> Q;
    long* mate;
    long* candidate;
    long total_weight_ld;
    vector<AssignE> assignEdgesLD;

    void process_vertex(long nid, long* mate, long* candidate);
public:
    long totalCost;
    Timer timer;

    LocalDominant(Graph* graph) {
        g = graph;
        totalCost = 0;
    }
    void runLocalDominant();
};


#endif //NETWORKFLOW_LOCALDOMINANT_H

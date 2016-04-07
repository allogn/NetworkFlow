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
    long from;
    long to;
    int cost;
    int cap;
    AssignE (long _from, long _to, int _cost, int _cap): from(_from), to(_to), cost(_cost), cap(_cap) {}
};

class LocalDominant {
    Graph* g;

    std::queue<long> Q;
    long* mate;
    long* candidate;
    vector<AssignE> assignEdgesLD;

    void process_vertex(long nid);
public:
    LocalDominant(Graph* graph) {
        g = graph;
        mate = new long[g->n];
        candidate = new long[g->n];
    }
    ~LocalDominant() {
        delete[] mate;
        delete[] candidate;
    }
    void runLocalDominant();
    long totalCost() {
        //calculate total cost and move result to AssignEdges
        long totalCost = 0;
        assert(g->n%2 == 0);
        //TODO works only for bipartite graphs now and unit capacity
        for (int i = 0; i < g->n/2; i++) {
            //from i to mate[i]
            for (int j = 0; j < g->completeE[i].size(); j++)
            {
                if (mate[i] == g->get_pair(g->completeE[i][j], i))
                {
                    int weight = g->E[g->completeE[i][j]].weight;
                    assignEdgesLD.push_back(AssignE(i,mate[i],weight,1)); //capacity = 1 @todo capacity
                    totalCost += weight;
                    break;
                }
            }
        }
        return totalCost;
    }
};


#endif //NETWORKFLOW_LOCALDOMINANT_H

//
// Created by alvis on 16.03.16.
//

#ifndef NETWORKFLOW_SIA_H
#define NETWORKFLOW_SIA_H

#include <stack>
#include <algorithm>

#include "Graph.h"
#include "nheap.h"
//TODO optimization parameter for compiler
class SIA {
    Graph* g;

    int noA, noB;
    int totalflow;
    unsigned long currentcost;
    long taumax;

    int* free;
    int* flow;
    int* nodeFlow;
    int* QryCnt;
    long* psi;

    //public
    long* mindist;
    int* mineid;
    int* watched;

    mmHeap dijkH;
    mmHeap globalH;
    mmHeap updateH;

    vector<int> worklist; //list of nodes needed for current execution

    void iteration_reset(int nodeid_best_psi, int source_id); // worklist reset after one iteration
    void augmentFlow(int lastid);
    int insertEdgeFromHeap();
    void processId(int source_id);
    int runDijkstra(int source_id);
    inline long getCost(int eid, int fromid, int toid)
    { return (fromid<noA) ? g->E[eid].weight-psi[fromid]+psi[toid] : -g->E[eid].weight-psi[fromid]+psi[toid]; };
    inline long getEdgeCost(int edge_w, int fromid)
    {
        long new_cost = edge_w + mindist[fromid];
        return new_cost;
    }
    inline int heap_checkAndUpdateEdgeMin(mmHeap& heap, int fromid) // update if new value is less
    {
        if (fromid >= noA) return 0;
        if (heap.isExisted(fromid))
        {
            int weight = g->E[g->fullE[fromid][QryCnt[fromid]-1]].weight;
            heap.updatequeue(fromid, getEdgeCost(weight,fromid));
            return 1;
        }

        if (QryCnt[fromid] < noB) {
            heap.enqueue(fromid,getEdgeCost(g->E[g->fullE[fromid][QryCnt[fromid]]].weight,fromid));
            QryCnt[fromid]++;
        }
        return 0;
    }
    int heap_checkAndUpdateMin(mmHeap& heap, int id, int new_value) // update if new value is less
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
    int updateMinDist(int eid, int fromid, int toid)
    {
        long cost = getCost(eid,fromid,toid);
        if (mindist[toid]>mindist[fromid]+cost) {
            mindist[toid] = mindist[fromid] + cost;
            mineid[toid] = eid;
            if (toid < noA)
            {
                heap_checkAndUpdateEdgeMin(globalH, toid);
            }
            return 1;
        }
        return 0;
    }
    void updateHeaps(int eid, int fromid, int toid)
    {
        if (watched[fromid] == 0) return; //case when prefinal node : not all neighbours are considered => not watched,
        // but in globalH and in DijkH (!). so, if in dijkH => everything is fine (will be updated later)
        if (updateMinDist(eid,fromid,toid)) {
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
    int reserveNode(int node_id) {
        bool exist = false;
        for (int i = 0; i < worklist.size(); i++)
        {
            if (worklist[i] == node_id)
            {
                exist = true;
                break;
            }
        }
        if (!exist) worklist.push_back(node_id);
        return 0;
    }

public:
    intT totalCost;

    SIA(Graph* graph) {
        g = graph;
        reset();
    }
    ~SIA() {
        delete[] free;
        delete[] flow;
        delete[] nodeFlow;
        delete[] psi;
        delete[] QryCnt;
        delete[] mindist;
        delete[] mineid;
        delete[] watched;
    }

    void reset();

    int runOSIA() {
        cout << "Running SIA..." << endl;

        assert(g->test_sorting());
        reset();

        while (totalflow < noA)
        {
            int check;
            for(int i = 0; i < noA; i++)
            {
                check=free[i];
                free[i]|=1;

                if (check == 0) {
                    processId(i);
                    break;
                }
            }
        }


        //calculate total cost
        for (int i = 0; i < noB; i++)
        {
            currentcost += g->E[g->V[i+noA].E[0]].weight;
        }

        totalCost = currentcost;
    }

    // Unit tests
    bool test_has_path(int nodeId); // tests if nodeId is full (has path from Q to P)
    bool test_mineid_path_exist(uintT target_node, uintT source_node);
};

#endif //NETWORKFLOW_SIA_H

//
// Created by alvis on 16.03.16.
//

#ifndef NETWORKFLOW_SIA_H
#define NETWORKFLOW_SIA_H

#include <stack>
#include <algorithm>

#include "Graph.h"
#include "nheap.h"
#include "TimerTool.h"

//TODO optimization parameter for compiler
class SIA {
    Graph* g;

    long noA, noB;
    long long totalflow;
    long taumax;

    int* flow; //residual flow in an edge
    long* excess;
    long long* psi;

    long long* mindist;
    long* mineid;
    bool* watched;

    mmHeap dijkH;
    mmHeap globalH;
    mmHeap updateH;

    inline void iteration_reset(long nodeid_best_psi);
    inline void augmentFlow(long lastid);
    inline long insertEdgeFromHeap();
    void processId(long source_id);
    inline long runDijkstra();
    inline long getCost(long eid, long fromid, long toid) {
        return (fromid<noA) ? g->E[eid].weight-psi[fromid]+psi[toid] : -g->E[eid].weight-psi[fromid]+psi[toid];
    };
    inline long getEdgeCost(long edge_w, long fromid) {
        return edge_w + mindist[fromid];
    }
    inline long heap_checkAndUpdateEdgeMin(mmHeap& heap, long fromid) // update if new value is less
    {
        if (fromid >= noA) return 0;
        long weight = g->get_next_neighbour_weight(fromid);

        if (weight != -1) { //not full
            if (heap.isExisted(fromid))
            {
                heap.updatequeue(fromid, getEdgeCost(weight,fromid));
                return 1;
            }
            heap.enqueue(fromid,getEdgeCost(weight,fromid));
        } else {
            assert(!heap.isExisted(fromid)); //not full if exists in the heap
        }
        return 0;
    }
    long heap_checkAndUpdateMin(mmHeap& heap, long id, long new_value) // update if new value is less
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
    inline bool updateMinDist(long eid, long fromid, long toid)
    {
        long cost = getCost(eid,fromid,toid);
        if (mindist[toid]>mindist[fromid]+cost) {
            mindist[toid] = mindist[fromid] + cost;
            mineid[toid] = eid;
            if (toid < noA)
            {
                heap_checkAndUpdateEdgeMin(globalH, toid);
            }
            return true;
        }
        return false;
    }
    void updateHeaps(long eid, long fromid, long toid)
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

public:
    long totalCost;
    Timer timer;

    SIA(Graph* graph) {
        g = graph;
        assert(g->test_sorting());
        reset();
    }
    ~SIA() {
        delete[] flow;
        delete[] excess;
        delete[] psi;
        delete[] mindist;
        delete[] mineid;
        delete[] watched;
    }

    void reset();

    void runOSIA() {
        reset();
        double total = timer.getTime();
        long q = 0;
        while (totalflow < noA) {
            processId(q);
            q++;
            if (q >= noA) q = 0;
        }
        timer.save_time("Total time", total);

        //calculate total cost
        long long currentcost = 0;
        for (long i = 0; i < noB; i++) {
            currentcost += g->E[g->V[i+noA].E[0]].weight;
        }
        totalCost = currentcost;
    }

    // Unit tests
    bool test_has_path(long nodeId); // tests if nodeId is full (has path from Q to P)
    bool test_mineid_path_exist(long target_node, long source_node);
};

#endif //NETWORKFLOW_SIA_H

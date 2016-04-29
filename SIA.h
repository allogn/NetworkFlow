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

    int _param;

    //profiling counters
    long total_iterations;
    long total_dijkstra;

    inline void iteration_reset(long nodeid_best_psi);
    inline void augmentFlow(long lastid);
    inline long insertEdgeFromHeap();
    void processId(long source_id);
    inline long runDijkstra();
    inline long getCost(long eid, long fromid, long toid) {
        return (fromid<noA) ? g->E[eid].weight-psi[fromid]+psi[toid] : -g->E[eid].weight-psi[fromid]+psi[toid];
    };
    inline long getEdgeCost(long edge_w, long fromid) {
        if (_param == 1) {
            return edge_w - psi[fromid];
        }
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
//            cout << "updated distance to " << toid << ": " << mindist[fromid]+cost << endl;
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

    SIA(Graph* graph, int param = 0) {
        g = graph;
        _param = param;
        assert(g->isSpatial || g->test_sorting());
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

    void save_profile_data(string filename, long experiment_id) {
        //calculate data about a graph
        ofstream outf(filename,ios::app);

        //mean and std of how much edges were added
        double perc = 0;
        double percSqr = 0;
        double totalNodes = 0;
        for (long i = 0; i < noA; i++) {
            long added = 0;
            //calculate how much inversed edges points to i-th node
            for (long j = noA; j < g->n; j++) {
                for (long k = 0; k < g->V[j].E.size(); k++) {
                    long eid = g->V[j].E[k];
                    if (g->E[eid].fromid == i || g->E[eid].toid == i) {
                        added++;
                    }
                }
            }
            added += g->V[i].E.size();
            double val = (double)added/(double)g->fullE[i].size();
            perc += val;
            assert(val >= 0 && val <= 1);
            percSqr += val*val;
        }
        outf << experiment_id << ",Saturation mean," << perc/(double)noA << "\n";
        outf << experiment_id << ",Saturation std," << sqrt(percSqr/(double)noA - (perc/(double)noA)*(perc/(double)noA)) << "\n";

        outf << experiment_id << ",Dijkstra executions," << total_dijkstra << "\n";
        outf << experiment_id << ",Iterations," << total_iterations << "\n";

        outf << experiment_id << ",Total cost," << totalCost << "\n";
        outf.close();
    }

    // Unit tests
    bool test_has_path(long nodeId); // tests if nodeId is full (has path from Q to P)
    bool test_mineid_path_exist(long target_node, long source_node);
    bool test_correct_flows();
};

#endif //NETWORKFLOW_SIA_H

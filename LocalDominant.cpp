//
// Created by alvis on 19.03.16.
//

#include "LocalDominant.h"

void LocalDominant::runLocalDominant() {
    double totalTime = getclock();

    long total = g->n;
    long mate[total];
    long candidate[total];
    fill(mate, mate+total, -1);
    fill(candidate, candidate+total, -1);

    for (int i = 0; i < total; i++)
    {
        process_vertex(i, mate, candidate);
    }

    long u, v;
    while (Q.size() > 0)
    {
        u = Q.front();
        Q.pop();
        for (long i = 0; i < g->fullE[u].size(); i++)
        {
            v = g->get_pair(g->fullE[u][i],u);
            if (mate[u] != v && candidate[v] == u)
            {
                process_vertex(v, mate, candidate);
            }
        }
    }

    timer.save_time("Total Local Dominant", totalTime);

    //calculate total cost and move result to AssignEdges
    totalCost = 0;
    assert(g->n%2 == 0);
    //TODO works only for bipartite graphs now and unit capacity
    for (int i = 0; i < g->n/2; i++) {
        //from i to mate[i]
        for (int j = 0; j < g->fullE[i].size(); j++)
        {
            if (mate[i] == g->get_pair(g->fullE[i][j], i))
            {
                long weight = g->E[g->fullE[i][j]].weight;
                assignEdgesLD.push_back(AssignE(i,mate[i],weight,1)); //capacity = 1 @todo capacity
                totalCost += weight;
                break;
            }
        }

    }
}

void LocalDominant::process_vertex(long nid, long *mate, long *candidate) {
    long w;
    long min_wt = -1;
    long min_wt_id = -1;
    long s = -1;
    for (long i = 0; i < g->fullE[nid].size(); i++)
    {
        long eid = g->fullE[nid][i];
        w = g->E[eid].weight;
        s = g->get_pair(eid, nid);
        if ((mate[s] == -1) && (min_wt > w || min_wt == -1))
        {
            min_wt = w;
            min_wt_id = s;
        }
    }

    candidate[nid] = min_wt_id;
    if (min_wt_id > -1 && candidate[candidate[nid]] == nid)
    {
        /* Check if locally dominant */

        mate[nid] = candidate[nid];
        mate[candidate[nid]] = nid;

        Q.push(nid);
        Q.push(candidate[nid]);
    }
}

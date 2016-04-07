//
// Created by alvis on 19.03.16.
//

#include "LocalDominant.h"

void LocalDominant::runLocalDominant() {

    long total = g->n;
    fill(mate, mate+total, -1);
    fill(candidate, candidate+total, -1);

    for (int i = 0; i < total; i++)
    {
        process_vertex(i);
    }

    long u, v;
    while (Q.size() > 0)
    {
        u = Q.front();
        Q.pop();
        for (long i = 0; i < g->completeE[u].size(); i++)
        {
            v = g->get_pair(g->completeE[u][i],u);
            if (mate[u] != v && candidate[v] == u)
            {
                process_vertex(v);
            }
        }
    }
}

void LocalDominant::process_vertex(long nid) {
    long w;
    long min_wt = -1;
    long min_wt_id = -1;
    long s = -1;
    for (long i = 0; i < g->completeE[nid].size(); i++)
    {
        long eid = g->completeE[nid][i];
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

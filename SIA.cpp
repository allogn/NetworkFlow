//
// Created by alvis on 16.03.16.
//

#include "SIA.h"

void SIA::reset() {

    if (g->n % 2 != 0) {
        cout << "Error: SIA algorithm applied on odd number of nodes" << endl;
        exit(1);
    }
    noA = g->n/2;
    noB = noA;

    totalflow = 0;
    currentcost = 0;
    taumax = 0;

    free = new int[noA];
    fill(free,free+noA,0);

    nodeFlow = new int[noA+noB];
    fill(nodeFlow,nodeFlow+noA,1);
    fill(nodeFlow+noA,nodeFlow+noA+noB,-1);

    flow = new int[g->m];
    fill(flow,flow+g->m,0);

    QryCnt = new int[noA];
    fill(QryCnt,QryCnt+noA,0);

    psi = new long[noA+noB];
    fill(psi,psi+noA+noB,0);

    mindist = new long[noA+noB];
    fill(mindist,mindist+noA+noB,LONG_MAX);

    mineid = new int[noA+noB];
    fill(mineid,mineid+noA+noB,-1);

    watched = new int[noA+noB];
    fill(watched,watched+noA+noB,0);
}

void SIA::augmentFlow(int lastid) {
    int curid = lastid;
    int eid, toid;

    //no one can use shared nodes because they should be locked
    while(mineid[curid] != -1)
    {
        //reverse edges
        eid = mineid[curid];

        toid = curid;
        curid = g->get_pair(eid,curid);

        int cur_flow = (curid < noA)?flow[eid]-1:flow[eid]+1;
        flow[eid] = cur_flow;

        //edge goes from curid->toid. we must delete it if necessary and add toid->curid
        if ((cur_flow == 0 && curid < noA) || (cur_flow == 1 && curid >= noA))
        {
            //delete it from curid
            int i = 0;
            while (g->get_pair(g->V[curid].E[i++],curid) != toid);
            g->V[curid].E[--i] = *g->V[curid].E.rbegin();
            g->V[curid].E.pop_back();
        }

        //add edge if necessary
        int i = 0;
        while (i < g->V[toid].E.size() && g->get_pair(g->V[toid].E[i++],toid) != curid);
        if (i == g->V[toid].E.size())
        {
            g->V[toid].E.push_back(eid);
        }
    }
}

int SIA::runDijkstra() {

    //reset mindist and watched and mineid
    if (dijkH.size() == 0) {
        printf("Error: zero dijkH at the beginning of dijkstra execution");
        exit(1);
    }

    int current_node = dijkH.getTopIdx();// heap_dequeue(dijkH);
    int source_id = source_id;

    //access node without troubles - no one will access it simultaneously
    while(1)
    {
        //first node - from source always for each dijkstra

        //traverse all neighbours and update alpha
        //traverse all neighbours and update alpha
        for (int i = 0; i < g->V[current_node].E.size(); i++)
        {
            //check if target node is busy
            int edge_id = g->V[current_node].E[i];
            int node_id = g->get_pair(edge_id,current_node);

            //add current neighbour to heap or update the heap
            int isUpdated = updateMinDist(edge_id,current_node,node_id); //if new => dist=inf and updated. if new but dist!=inf => in DijkH exist and also will be updated if updAlpha
            //not updated -> not enqueue

//            if (watched[node_id] == 1) continue;

            if (isUpdated) dijkH.enqueue(node_id,mindist[node_id]);
        }

        if (nodeFlow[current_node] < 0)
        {
            return current_node;
        }

//            if (current_node < noA)
//                enqueueNextEdge(sc,pc,current_node);
        heap_checkAndUpdateEdgeMin(globalH,current_node);
        watched[current_node] = 1;

        dijkH.dequeue();
        if (dijkH.size() == 0)
        {
            current_node = -1;
            break;
        }
        current_node = dijkH.getTopIdx();
    }
}

int SIA::insertEdgeFromHeap()
{
    int fromid, tmp;
    globalH.dequeue(fromid, tmp);


    //create new edge
//    Edge e;
////    e.flow = 1; // TODO why flow = 1????
//    e.fromid = fromid;
//    e.toid = g->V[fromid].fullE[sc.QryCnt[fromid]-1].node_pair_id;
//    e.maxflow = 1;
//    e.weight = g->V[fromid].fullE[sc.QryCnt[fromid]-1].weight;
    uintT eid = g->fullE[fromid][QryCnt[fromid]-1];
//    uintT eid = g->get_next_neighbour(fromid);
    flow[eid] = 1;
//
//    int eid;
//    eid = g->E.size();
//    sc.edges.push_back(e);

//    g->V[g->E[eid].fromid].E.push_back(eid);
    heap_checkAndUpdateEdgeMin(globalH,g->E[eid].fromid); //MUST BE HERE
    return eid;
}

void SIA::iteration_reset(int nodeid_best_psi, int source_id) {

    long best_psi = mindist[nodeid_best_psi];
    int nodeid;

    dijkH.clear();
    updateH.clear();

    while (worklist.size() > 0)
    {
        //release all nodes. note - source node is not in this list
        nodeid = worklist.back();
        worklist.pop_back();

        //UPDATE POTENTIALS
        if (watched[nodeid] == 1 && mindist[nodeid] < best_psi)
        {
            psi[nodeid] = psi[nodeid] - mindist[nodeid] + best_psi; //@todo check if parallel reset if potentials are correct
            watched[nodeid] = 0; //here because multiple values in worklist
        }
        //END OF UPDATE POTENTIALS
    }

    int edgeid;

    //do not release source (if it is source => it is not attached), but update it's potential
    psi[source_id] = psi[source_id] - mindist[source_id] + best_psi;

    fill(mindist,mindist+noA+noB,LONG_MAX);
    fill(mineid,mineid+noA+noB,-1);
    fill(watched,watched+noA+noB,0);
}

void SIA::processId(int source_id)
{
    int target_node = -1;
    mindist[source_id] = 0;
    watched[source_id] = 1;
    heap_checkAndUpdateEdgeMin(globalH,source_id);
    dijkH.enqueue(source_id,0);

    while((target_node == -1) || (globalH.size() > 0 && mindist[target_node]>globalH.getTopValue() - taumax))
    {
        //add edge
        int eid;

        int once = 0;
        while ((dijkH.size() == 0 && target_node == -1) || once == 0)//  || (globalH.size() > 0 && heap_getTopValue(dijkH)>heap_getTopValue(globalH))) //dijkH check because next in globalH can be not from current "thread", but must be added because it will be needed later
        {
            once = 1;
            eid = insertEdgeFromHeap();
            int fromid = g->E[eid].fromid;
            int toid = g->E[eid].toid;
            //@todo
            bool exist = false;
            for (int i = 0; i < worklist.size(); i++)
            {
                if (worklist[i] == toid)
                {
                    exist = true;
                    break;
                }
            }
            if (!exist) worklist.push_back(toid); //@todo maybe delete this?

            /*
             * Dijkstra Updates
             */
            updateHeaps(eid, g->E[eid].fromid, g->E[eid].toid);

            //update Heaps @todo race condition! for parallel

            int curid, tmp;
            while (updateH.size() > 0) {
                updateH.dequeue(curid, tmp);
                //            isUpdated[curid] = 0;
                for (int i = 0; i < g->V[curid].E.size(); i++) {
                    eid = g->V[curid].E[i];
                    toid = g->get_pair(eid, curid);
                    updateHeaps(eid, curid, toid);
                }
            }
            /*
             * End of dijkstra updates. If switch off => add clearing dijkstra in dijksra beginning
             */
        }

        target_node = runDijkstra();

        int taumax = 0;
        for (int i = 0; i<noA; i++)
        {
            if (mindist[i]<globalH.getTopValue() && taumax < psi[i])
            {
                taumax = psi[i];
            }
        }
    };

    augmentFlow(target_node);

    /* Flow update BEFORE releasing! */
    nodeFlow[source_id]--;
    nodeFlow[target_node]++;

    //#pragma omp atomic update
    totalflow++;

    //release worklist
    iteration_reset(target_node, source_id);
}


/*
 * Unit tests
 */

bool SIA::test_has_path(int nodeId) {

}
//
// Created by alvis on 16.03.16.
//

#include "SIA.h"

void SIA::reset() {
    globalH.clear();
    dijkH.clear();
    updateH.clear();

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
    assert(dijkH.size() != 0);
    int current_node = dijkH.getTopIdx();

    //access node without troubles - no one will access it simultaneously
    while(1)
    {

        if (nodeFlow[current_node] < 0)
        {
            return current_node;
        }//todo why here not before looping through edges?
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

            if (isUpdated) {
                if (!dijkH.isExisted(node_id)) {
                    dijkH.enqueue(node_id,mindist[node_id]);
                } else {
                    dijkH.updatequeue(node_id,mindist[node_id]);
                }
            } //todo if object exists already?
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
    return current_node;
}

int SIA::insertEdgeFromHeap()
{
    int fromid, tmp;
    globalH.dequeue(fromid, tmp);

    assert(QryCnt[fromid] < g->fullE[fromid].size());
    uintT eid = g->fullE[fromid][QryCnt[fromid]];
    QryCnt[fromid]++;
    flow[eid] = 1; //one because an edge from noA and noB has different meanings on flow, it is residual flow

    g->V[g->E[eid].fromid].E.push_back(eid);
    heap_checkAndUpdateEdgeMin(globalH,g->E[eid].fromid); //MUST BE HERE
    return eid;
}

void SIA::iteration_reset(int nodeid_best_psi) {
    long best_psi = mindist[nodeid_best_psi];
    dijkH.clear();
    updateH.clear();
    globalH.clear();
    for (int nodeid = 0; nodeid < g->n; nodeid++) {
        if (mindist[nodeid] < best_psi)
        {
            psi[nodeid] = psi[nodeid] - mindist[nodeid] + best_psi; //@todo check if parallel reset if potentials are correct
        }
    }
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

        do //  || (globalH.size() > 0 && heap_getTopValue(dijkH)>heap_getTopValue(globalH))) //dijkH check because next in globalH can be not from current "thread", but must be added because it will be needed later
        {
            eid = insertEdgeFromHeap();

            int toid = g->E[eid].toid;
            /*
             * Dijkstra Updates
             */
            updateHeaps(eid, g->E[eid].fromid, g->E[eid].toid);

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
        } while (dijkH.size() == 0 && target_node == -1);

        target_node = runDijkstra();

        int gettaumax = 0; //todo maintain tau
        for (int i = 0; i<noA; i++)
        {
            if (mindist[i]<globalH.getTopValue() && gettaumax < psi[i])
            {
                gettaumax = psi[i];
            }
        }
        taumax = gettaumax;
    };
    augmentFlow(target_node);

    /* Flow update BEFORE releasing! */
    nodeFlow[source_id]--;
    nodeFlow[target_node]++;

    totalflow++;

    //release worklist
    iteration_reset(target_node);
}


/*
 * Unit tests
 */

bool SIA::test_has_path(int nodeId) {
    stack<uintT> queue;
    bool visited[g->n];
    fill(visited, visited+g->n, false);
    for (uintT i = noA; i < g->n; i++) {
        queue.push(i);
    }
    while (queue.size() > 0) {
        uintT curid = queue.top();
        queue.pop();
        if (curid == nodeId) return true;
        if (visited[curid] == true) continue;
        visited[curid] = true;
        for (uintT i = 0; i < g->V[curid].E.size(); i++) {
            queue.push(g->get_pair(g->V[curid].E[i],curid));
        }
    }
    return false;
}

bool SIA::test_mineid_path_exist(uintT target_node, uintT source_node) {
    cout << "Testing if path from target node exists..." << endl;
    uintT curnode = target_node;
    while(curnode != source_node) {
        assert(mineid[curnode] != -1);
        uintT eid = mineid[curnode];
        uintT neighbor = g->get_pair(eid, curnode);
        //check if eid is in E of neighbour
        assert(find(g->V[neighbor].E.begin(), g->V[neighbor].E.end(), eid) != g->V[neighbor].E.end());
        curnode = neighbor;
    }
    return true;
}

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
    taumax = 0;

    excess = new long[noA+noB];
    fill(excess,excess+noA,1);
    fill(excess+noA,excess+noA+noB,-1);

    flow = new int[g->m];
    fill(flow,flow+g->m,1);

    psi = new long long[noA+noB];
    fill(psi,psi+noA+noB,0);

    mindist = new long long[noA+noB];
    fill(mindist,mindist+noA+noB,LONG_MAX);

    mineid = new long[noA+noB];
    fill(mineid,mineid+noA+noB,-1);

    watched = new bool[noA+noB];
    fill(watched,watched+noA+noB,0);

    total_iterations = 0;
    total_dijkstra = 0;
}

void SIA::augmentFlow(long lastid) {
    long curid = lastid;
    long eid, toid;
    while(mineid[curid] != -1)
    {
        //reverse edges
        eid = mineid[curid];
        toid = curid;
        curid = g->get_pair(eid,curid);

        int cur_flow = (curid < noA)?flow[eid]-1:flow[eid]+1;
        assert(cur_flow == 0 || cur_flow == 1);
        flow[eid] = cur_flow;

        //deleting edges with no residual cost and moving them to the list of a reverse edge
        if ((cur_flow == 0 && curid < noA) || (cur_flow == 1 && curid >= noA))
        {
            //delete edge from curid node edge list
            long i = 0;
            while (g->get_pair(g->V[curid].E[i++],curid) != toid);
            g->V[curid].E[--i] = *g->V[curid].E.rbegin();
            g->V[curid].E.pop_back();
        }

        //adding edge to a reverse edge if necessary
        long i = 0;
        while (i < g->V[toid].E.size() && g->get_pair(g->V[toid].E[i++],toid) != curid);
        if (i == g->V[toid].E.size()) {
            g->V[toid].E.push_back(eid);
        }
    }
}

long SIA::runDijkstra() {
    total_dijkstra++;
    assert(dijkH.size() != 0);
    long current_node = dijkH.getTopIdx();

    while(dijkH.size()>0 && (globalH.size()==0 || dijkH.getTopValue()<=globalH.getTopValue()))
    {
        //return of deficit is reached
        if (excess[current_node] < 0) {
            return current_node;
        }

        //traverse all neighbours and update alpha
        for (long i = 0; i < g->V[current_node].E.size(); i++) {
            //check if target node is busy
            long edge_id = g->V[current_node].E[i];
            long node_id = g->get_pair(edge_id,current_node);

            //updating minimum distance for Dijkstra and updating globalH inside
            bool isUpdated = updateMinDist(edge_id,current_node,node_id);
            if (isUpdated) {
                if (!dijkH.isExisted(node_id)) {
                    dijkH.enqueue(node_id,mindist[node_id]);
                } else {
                    dijkH.updatequeue(node_id,mindist[node_id]);
                }
            }
        }
        heap_checkAndUpdateEdgeMin(globalH,current_node);
        watched[current_node] = true;

        dijkH.dequeue();
        if (dijkH.size() == 0) {
            current_node = -1;
            break;
        }
        current_node = dijkH.getTopIdx();

    }
    return -1;
}

long SIA::insertEdgeFromHeap()
{
    long fromid;
    long long tmp;
    globalH.dequeue(fromid, tmp);

    long eid = g->insert_nn_to_edge_list(fromid);
    assert(eid >= 0); // not full
    flow[eid] = 1; //one because an edge from noA and noB has different meanings on flow, it is residual flow

    g->V[g->E[eid].fromid].E.push_back(eid);
    heap_checkAndUpdateEdgeMin(globalH,g->E[eid].fromid); //MUST BE HERE
    return eid;
}

void SIA::iteration_reset(long nodeid_best_psi) {
    long best_psi = mindist[nodeid_best_psi];
    dijkH.clear();
    updateH.clear();
    globalH.clear();
    for (long nodeid = 0; nodeid < g->n; nodeid++) {
        if (mindist[nodeid] < best_psi) {
            psi[nodeid] = psi[nodeid] - mindist[nodeid] + best_psi;
        }
    }
    fill(mindist,mindist+noA+noB,std::numeric_limits<long long>::max());
    fill(mineid,mineid+noA+noB,-1);
    fill(watched,watched+noA+noB,0);
}

void SIA::processId(long source_id)
{
    total_iterations = 0;
    double time_start = timer.getTime();
    long target_node = -1;
    mindist[source_id] = 0;
    watched[source_id] = 1;
    heap_checkAndUpdateEdgeMin(globalH,source_id);
    dijkH.enqueue(source_id,0);

    assert(test_correct_flows());
    while((target_node == -1) || (globalH.size() > 0 && mindist[target_node]>globalH.getTopValue() - taumax))
    {
        do //  || (globalH.size() > 0 && heap_getTopValue(dijkH)>heap_getTopValue(globalH)))
        {
            //adding new edge to the subgraph
            long eid = insertEdgeFromHeap();

            //updating Dijkstra heap by using additional heap
            updateHeaps(eid, g->E[eid].fromid, g->E[eid].toid);
            long curid;
            long long tmp;
            while (updateH.size() > 0) {
                updateH.dequeue(curid, tmp);
                for (long i = 0; i < g->V[curid].E.size(); i++) {
                    eid = g->V[curid].E[i];
                    updateHeaps(eid, curid, g->get_pair(eid, curid));
                }
            }
        } while (dijkH.size() == 0);

        double time_dijkstra = timer.getTime();
        target_node = runDijkstra();
        timer.save_time("Dijkstra", time_dijkstra);

        long long gettaumax = 0;
        for (long i = 0; i<noA; i++)
        {
            if (mindist[i]<globalH.getTopValue() && gettaumax < psi[i])
            {
                gettaumax = psi[i];
            }
        }
        taumax = gettaumax;
    };
    double time_augment = timer.getTime();
    augmentFlow(target_node);

    //assuming one iterations makes exactly one node not free
    excess[source_id]--;
    excess[target_node]++;
    totalflow++;

    //raise potentials and clear everything after iteration
    iteration_reset(target_node);

    //log time
    timer.save_time("Process node", time_start);
    timer.save_time("Augment and reset", time_augment);
}

void SIA::processIdApprox(long source_id)
{
    //maybe it WILL be exact if globalH resets each iteration

    total_iterations = 0;
    double time_start = timer.getTime();
    long target_node = -1;
    mindist[source_id] = 0;
    watched[source_id] = 1;
    heap_checkAndUpdateEdgeMin(globalH,source_id);
    dijkH.enqueue(source_id,0);
    assert(test_correct_flows());
    while((target_node == -1))// || (globalH.size() > 0 && mindist[target_node]>globalH.getTopValue() - taumax))
    {
        do //  || (globalH.size() > 0 && heap_getTopValue(dijkH)>heap_getTopValue(globalH)))
        {
            //adding new edge to the subgraph
            long eid = insertEdgeFromHeap();

            //updating Dijkstra heap by using additional heap
            updateHeaps(eid, g->E[eid].fromid, g->E[eid].toid);
        } while (dijkH.size() == 0 || (globalH.size() > 0 && dijkH.getTopValue()>globalH.getTopValue()));

        long curid;
        long long tmp;
        while (updateH.size() > 0) {
            updateH.dequeue(curid, tmp);
            for (long i = 0; i < g->V[curid].E.size(); i++) {
                long eid = g->V[curid].E[i];
                updateHeaps(eid, curid, g->get_pair(eid, curid));
            }
        }

        double time_dijkstra = timer.getTime();
        target_node = runDijkstra();
        timer.save_time("Dijkstra", time_dijkstra);

        long long gettaumax = 0;
        for (long i = 0; i<noA; i++)
        {
            if (mindist[i]<globalH.getTopValue() && gettaumax < psi[i])
            {
                gettaumax = psi[i];
            }
        }
        taumax = gettaumax;
    };
    double time_augment = timer.getTime();
    augmentFlow(target_node);

    //assuming one iterations makes exactly one node not free
    excess[source_id]--;
    excess[target_node]++;
    totalflow++;

    //raise potentials and clear everything after iteration
    iteration_reset(target_node);

    //log time
    timer.save_time("Process node", time_start);
    timer.save_time("Augment and reset", time_augment);
}


/*
 * Unit tests
 */

bool SIA::test_has_path(long nodeId) {
    cout << "Test id node has path" << endl;
    stack<long> queue;
    bool visited[g->n];
    fill(visited, visited+g->n, false);
    for (long i = noA; i < g->n; i++) {
        queue.push(i);
    }
    while (queue.size() > 0) {
        long curid = queue.top();
        queue.pop();
        if (curid == nodeId) return true;
        if (visited[curid]) continue;
        visited[curid] = true;
        for (long i = 0; i < g->V[curid].E.size(); i++) {
            queue.push(g->get_pair(g->V[curid].E[i],curid));
        }
    }
    return false;
}

bool SIA::test_mineid_path_exist(long target_node, long source_node) {
    cout << "Testing if path from target node exists..." << endl;
    long curnode = target_node;
    while(curnode != source_node) {
        assert(mineid[curnode] != -1);
        long eid = mineid[curnode];
        long neighbor = g->get_pair(eid, curnode);
        //check if eid is in E of neighbour
        assert(find(g->V[neighbor].E.begin(), g->V[neighbor].E.end(), eid) != g->V[neighbor].E.end());
        curnode = neighbor;
    }
    return true;
}

bool SIA::test_correct_flows() {
    for (long i = 0; i < g->V.size(); i++) {
        for (long j = 0; j < g->V[i].E.size(); j++) {
            long eid = g->V[i].E[j];
            assert((i < noA && flow[eid] == 1) || (i >= noA && flow[eid] == 0));
        }
    }

    return true;
}

//
// Created by alvis on 18.03.16.
//

#include "CostScaling.h"
//TODO check where admissible set
void CostScaling::reset() {
    // init epsilon
    long maxWeight = 0;
    for(long i = 0; i < g->n; i++) {
        for (long j = 0; j < g->V[i].E.size(); j++) {
            maxWeight = (maxWeight < g->E[g->V[i].E[j]].weight) ? g->E[g->V[i].E[j]].weight : maxWeight;
        }
    }
    epsilon = maxWeight;

    // other init
    long totalPosSup = 0;
    for(long i = 0; i < g->n; i++) {
        p[i] = 0;
        excesses[i] = g->V[i].supply;
        if (excesses[i] > 0)
            totalPosSup++;
    }
    totalExcesses = totalPosSup;


    for(long i = 0; i < g->m; i++) {
            flow[i] = 0;
            admissible[i] = false;
    }
}

/*
 * edgeInNodeId here is id in the list of outgoing edges in the node
 *
 * if nodeId is fromId (source node) => value is added to the flow, substracted otherwise
 */
void CostScaling::changeFlow(long nodeId, long edgeInNodeId, long value) {
    assert(value > 0);
    if (value == 0) return;

    long eid = g->V[nodeId].E[edgeInNodeId];
    long neighbor = g->get_pair(eid, nodeId);
    long valueForAddition = (g->is_forward(eid, nodeId))?value:-value;
    assert(flow[eid] + valueForAddition <= g->E[eid].capacity
           && flow[eid] + valueForAddition >= g->E[eid].lower);

    // if flow is equal to capacity or lower - then it should be added to edge list of neighbor
    // because after flow addition the backward edge will be not saturated anymore
//    if (flow[eid] == g->E[eid].capacity || flow[eid] == g->E[eid].lower) {
        assert(((flow[eid] == g->E[eid].capacity) && !g->is_forward(eid, nodeId)) ||
                ((flow[eid] == g->E[eid].lower) && g->is_forward(eid, nodeId)));
        g->V[neighbor].E.push_back(eid);
        //TODO assert it does not exist there
        //TODO is not unit capacities - it may exist there, so add only if it does not
//    }

    flow[eid] += valueForAddition;
    // if flow reached capacity or lower - remove from the source
    //TODO if not unit capacities - add this check. for unit capacities it is always the truth
//    if (flow[eid] == g->E[eid].capacity || flow[eid] == g->E[eid].lower) {
        //assuming unit capacity here
        g->V[nodeId].E[edgeInNodeId] = *g->V[nodeId].E.rbegin();
        g->V[nodeId].E.pop_back();
    assert(((flow[eid] == g->E[eid].capacity) && !g->is_forward(eid, neighbor)) ||
           ((flow[eid] == g->E[eid].lower) && g->is_forward(eid, neighbor)));
//    }

    long old_val_node = excesses[nodeId];
    long old_val_neig = excesses[neighbor];

    excesses[nodeId] -= value;
    excesses[neighbor] += value;

    // calculate new total excess amount
    // assuming unit capacities here - stupid way
    short count = 0;
    if (old_val_node <= 0 && excesses[nodeId] > 0) count++;
    if (old_val_node > 0 && excesses[nodeId] <= 0) count--;
    if (old_val_neig <= 0 && excesses[neighbor] > 0) count++;
    if (old_val_neig > 0 && excesses[neighbor] <= 0) count--;

    totalExcesses += count;
    assert(test_excesses());
}

// returns nearest excess to deficit which was found at the end
// in blocking flow array save ID of an edge for each node throughout the path
void CostScaling::dfs(long* blockingFlow, long* blockingEdge, long nodeId) {
    //blockingFlow here also is a visited array to avoid loops
    queue<long> node_queue;
    node_queue.push(nodeId);

    while(!node_queue.empty()) {
        long node = node_queue.front();
        node_queue.pop();

        for (long i = 0; i < g->V[node].E.size(); i++) {
            //for each out edge check if it is admissible
            //if so - update distance and run dfs
            long neighbor = g->get_pair(g->V[node].E[i], node);
            if (isAdmissible(node,i) && blockingFlow[neighbor] == std::numeric_limits<long>::max()
                && excesses[neighbor] <= 0) {
                blockingFlow[neighbor] = node;
                blockingEdge[neighbor] = i;
                node_queue.push(neighbor);

                if (excesses[neighbor] < 0) {
                    return;
                }
            }
        }
    }
}

void CostScaling::raise_potentials() {

    // boolean visited array
    bool visited[g->n];

    //calculate shortest distances
    NodeList buckets(g->n);
    for (long i = 0; i < g->n; i++) {
        visited[i] = false;
        if (excesses[i] > 0) {
            buckets.addToBucket(0, i);
        }
    }

    //calculate shortest path in Residual (!) graph
    long long length_to_deficit;
    while (true) {
        long i = buckets.getNearest();

        if (visited[i]) {
            buckets.popNearest();
            continue;
        }
        visited[i] = true;
        buckets.saveAndPopNearest(); // save to visited inside buckets

        if (excesses[i] < 0) {
            //run algorithm until first deficit is reached
            length_to_deficit = buckets.smallest_dist;
            break;
        }

        //iterate through residual edges
        for (long j = 0; j < g->V[i].E.size(); j++) {
            long eid = g->V[i].E[j];
            long neighbour = g->get_pair(eid, i);
            if (! visited[neighbour]) { //@todo change if statements (with edge) and maybe faster
                if (get_residual(eid, i) > 0) { //TODO remove this check!!!! it implies that residual is not empty because it IS in E list!!!!
                    double c_p = get_cp(eid, i);
                    assert(floor(c_p/epsilon) + 1 >= 0);
                    long long length = floor(c_p/epsilon) + 1;
                    long long newDist = buckets.smallest_dist + length;
                    buckets.addToBucket(newDist, neighbour);
                }
            }
        }
    }

//    assert(test_shortest_path(buckets));

    //delete everything else from buckets without moving to "visited"
    buckets.clearBuckets(visited);

    //raise potentials for all visited vertices
    while (buckets.visited != NULL) {
        long i = buckets.getIdVisited();
        long long dist = buckets.getDistVisited();
        buckets.popVisited();
        p[i] += (length_to_deficit - dist)*epsilon;
    }

}

void CostScaling::raise_flows() {
    long* blockingSearch = new long[g->n];
    long* blockingSearchEdges = new long[g->n];
    for (long i = 0; i < g->n; i++) {
        blockingSearch[i] = std::numeric_limits<long>::max();
        blockingSearchEdges[i] = std::numeric_limits<long>::max();
    }
    for (long i = 0; i < g->n; i++) {
        if (excesses[i] > 0) {
            dfs(blockingSearch, blockingSearchEdges, i);
        }
    }

    //starting from the lastest excess returned from DFS
    //raise flow throughout all the path
    //each blockingFlow element contains ID of the edge that should be used next
    //if ID is greater than outDegree - then inverse edge should be used

    for (long i = 0; i < g->n; i++) {
        if (excesses[i] < 0 && blockingSearch[i] < std::numeric_limits<long>::max()) {
            long curnode = i;
            long nextnode;
            while (true) {
                nextnode = blockingSearch[curnode];
                long edge = blockingSearchEdges[curnode];
                assert(isAdmissible(nextnode, edge));
                if (excesses[nextnode] > 0) {
                    changeFlow(nextnode, edge, 1);
                    break;
                } else {
                    changeFlow(nextnode, edge, 1);
                }
                curnode = nextnode;
            }
        }
    }

    assert(test_excesses());
    free(blockingSearch);
}

void CostScaling::refine() {
    assert(is_feasible(-2*epsilon));
    assert(test_excesses());
    //loop through arcs, saturate those with c_p less than 0 and add excesses

    /*
     * Important: iteration by E.size() using for loop is wrong here - the flow must NOT change during such iteration
     * since the size and content of E is changing with a flow
     * todo check throughout the algorithm
     */
    //TODO we can avoid this looping by saving which edges have c_p<0 from shortest path of previous iteration
    for (long i = 0; i < g->n; i++) {
        long j = 0;
        while (j < g->V[i].E.size()) {
            if (isAdmissible(i, j)) {
                changeFlow(i, j, 1); //assuming unit capacity and edge is not saturated (otherwise must have been moved)
                //assuming if flow is changed - edge is moved - then j must not be increased because another
                // node is now in the place of the old one
            } else {
                j++;
            }
        }
    }

    assert(is_feasible(-epsilon));
    assert(test_excesses());

    while (totalExcesses > 0) {
//        long last_increased = 1;
//        while (last_increased > 0) {
        raise_potentials();
        assert(is_feasible(-epsilon));
//            last_increased = increase_subgraph(epsilon, GA);
//        }
//        assert(is_feasible(-epsilon,GA));
        raise_flows();
    }

    assert(is_flow());
}

void CostScaling::runCostScaling() {
    cout << "Running CostScaling..." << endl;

    reset();
    double total = timer.getTime();
    while (epsilon > (double)1./g->n) {
        refine();
        epsilon = epsilon / 2.;
    }
    timer.save_time("Total CostScaling", total);

    long long self_totalcost = 0;
    for (long eid = 0; eid < g->m; eid++) {
        self_totalcost += g->E[eid].weight*flow[eid];
    }
    totalCost = self_totalcost;
}

/*
 * Unit Tests
 */

bool CostScaling::is_feasible(double threshold) {
    for (long i = 0; i < g->n; i++) {
        for (long j = 0; j < g->V[i].E.size(); j++) {
            long eid = g->V[i].E[j];
            if (get_residual(eid, i) > 0) {
                long long c_p = get_cp(eid, i);
                assert(c_p >= threshold);
            }
        }
    }
    return true;
}

bool CostScaling::is_flow() {
    long* testExcesses = new long[g->n];
    for (long i = 0; i < g->n; i++) {
        testExcesses[i] = g->V[i].supply;
    }
    for (long i = 0; i < g->m; i++) {
        testExcesses[g->E[i].fromid]-=flow[i];
        testExcesses[g->E[i].toid]+=flow[i];
    }
    for (long i = 0; i < g->n; i++) {
        assert(testExcesses[i] == 0);
    }
    delete[] testExcesses;
    return true;
}

/*
 * Important note: node->E list contains only non-saturated edges, but all edges adjacent to the node
 * are deleted from this list. There is no list with in edges in current implementation.
 */
bool CostScaling::test_excesses() {
    long* testExcesses = new long[g->n];
    for (long i = 0; i < g->n; i++) {
        testExcesses[i] = g->V[i].supply;
    }
    long ctotalExcesses = 0;
    for (long i = 0; i < g->m; i++) {
        testExcesses[g->E[i].fromid]-=flow[i];
        testExcesses[g->E[i].toid]+=flow[i];
    }
    for (long i = 0; i < g->n; i++) {
        assert(testExcesses[i] == excesses[i]);
        if (testExcesses[i] > 0)
            ctotalExcesses++;
    }
    assert(totalExcesses == ctotalExcesses);
    delete[] testExcesses;
    return true;
}
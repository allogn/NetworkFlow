//
// Created by alvis on 22.03.16.
//
#include <iostream>
#include <stack>
#include "Graph.h"

using namespace std;

int main(int argc, const char** argv) {
    cout << "Generating .gr from .adj and adding weights and capacities" << endl;

    if (argc < 3) {
        cout << "Provide input and output" << endl;
        exit(1);
    }

    string inputfile(argv[1]);
    string outputfile(argv[2]);

    Graph g;
    g.load_adj_graph(inputfile); // this loads graph with unit capacities and unit weight

    // make weights random
    long seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    uniform_int_distribution<uintT> uniGen(1, 1000);
//    normal_distribution<double> gausGen(param1, param2);
//    exponential_distribution<double> expGen((double)param1/100.);
    parallel_for(long i = 0; i < g.m; i++) {
        g.E[i].weight = uniGen(generator);
    }
    g.init_neighbors();
    // set up source and target and test if problem is feasible. select another if not
    uniform_int_distribution<long> uniGenNode(0, g.n-1);
    bool feasible = false;
    long node1;
    long node2;
    do {
        cout << "Choosing another source and target..." << endl;
        //choose another nodes for source and target
        node1 = uniGenNode(generator);
        do {
            node2 = uniGenNode(generator);
        } while (node2 == node1);

        //simple DFS
        stack<long> q;
        bool visited[g.n];
        fill(visited, visited+g.n, false);
        q.push(node1);
        while (q.size() > 0) {
            long curnode = q.top();
            q.pop();
            if (visited[curnode]) continue;
            visited[curnode] = true;
            if (curnode == node2) {
                feasible = true;
                break;
            }
            for (long i = 0; i < g.fullE[curnode].size(); i++) {
                q.push(g.get_pair(g.fullE[curnode][i], curnode));
            }
        }
    } while (!feasible);
    g.V[node1].supply = 1;
    g.V[node2].supply = -1;
    g.save_graph(outputfile);

    cout << inputfile << " to " << outputfile << " done." << endl;
}
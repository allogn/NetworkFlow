//
// Created by alvis on 16.03.16.
//

#include <iostream>

#include "Common.h"
#include "Graph.h"

int main(int argc, const char** argv) {

    Graph g;
    g.generate_full_bipartite_graph(100, 0, 1000);
    g.test_graph_structure();
    g.init_neighbors();
    g.sort_neighbors();
    g.test_sorting();

    g.test_save_load();

    cout << "Testing done" << endl;
}
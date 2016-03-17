//
// Created by alvis on 16.03.16.
//

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>

#include "Common.h"
#include "Graph.h"
#include "SIA.h"

int main(int argc, const char** argv) {

    Graph g;
    g.generate_full_bipartite_graph(100, 0, 1000);
    g.test_graph_structure();
    g.init_neighbors();
    g.sort_neighbors();
    g.test_sorting();

    g.test_save_load();

    // load answer sheet and compare results
    cout << "Compare results with answer sheet..." << endl;
    string ansfname = "AnswerSheet.txt";
    ifstream infile(ansfname);
    if (FILE *file = fopen(ansfname.c_str(), "r")) {
        fclose(file);
    } else {
        cout << "File " << ansfname << " does not exist" << endl;
        exit(1);
    }
    string graphfile;
    long answer;

    while(infile >> graphfile >> answer) {
        cout << "Checking file " << graphfile << " with answer " << answer << "..." << endl;
        g.clear_graph();
        string path = "../../data/NetworkFlowTests/bipartite/";
        path.append(graphfile);
        path[path.size()-1] = 'r';
        path[path.size()-2] = 'g';
        g.load_graph(path);
        g.init_neighbors();
        g.sort_neighbors();
        SIA SIAsolv(&g);
        SIAsolv.runOSIA();
        assert(SIAsolv.totalCost == answer);
    }

    cout << "Testing done" << endl;
}
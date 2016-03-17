//
// Created by alvis on 17.03.16.
//

#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "Common.h"
#include "Graph.h"
#include "SIA.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, const char** argv) {
    std::cout << "Graph Generator (c) Alvis Logins 2016" << std::endl;

    // parsing parameters
    int algorithm;
    int graph_type;
    int distr;
    string outf;
    string outfBlossom;
    uintT n,m,min_weight,max_weight;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce help message")
            ("graph,g", po::value<int>(&graph_type)->required(), "graph type 0:bipartite")
            ("output,o", po::value<string>(&outf)->required(), "save generated graph")
            ("blossom", po::value<string>(&outfBlossom)->required(), "save generated graph for blossom")
            ("nodes,n", po::value<uintT>(&n)->required(), "number of nodes")
            ("edges,m", po::value<uintT>(&m), "number of edges (default: clique)")
            ("distr,d", po::value<int>(&distr)->default_value(0), "type of distribution for weights (0:u,1:g,2:e)")
            ("minweight,l", po::value<uintT>(&min_weight)->default_value(0), "parameter1 for distribution")
            ("maxweight,u", po::value<uintT>(&max_weight)->default_value(1000), "parameter2 for distribution");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    po::notify(vm);


    Graph g;
    switch(graph_type) {
        case 0:
            g.generate_full_bipartite_graph(n, min_weight, max_weight, distr);
            break;
    }

    g.save_graph(outf);
    g.save_graph_blossom(outfBlossom);
    cout << "Done." << endl;
}
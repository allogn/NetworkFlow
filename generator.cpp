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
    std::cout << ":: Graph Generator ::" << std::endl;

    // parsing parameters
    int algorithm;
    int graph_type;
    int distr;
    long sources;
    long targets;
    string outf;
    string outfBlossom;
    long n,m,min_weight,max_weight;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce help message")
            ("graph,g", po::value<int>(&graph_type)->required(), "graph type (0:bipart,1:dense-sink)")
            ("output,o", po::value<string>(&outf)->required(), "save generated graph")
            ("blossom", po::value<string>(&outfBlossom)->required(), "filename for generated graph for blossom")
            ("nodes,n", po::value<long>(&n)->required(), "number of nodes")
            ("sources,s", po::value<long>(&sources)->default_value(1), "number of sources (if not bipartite)")
// capacity of each target = 1, capacity of a source depends on targets, equally distributed
            ("targets,t", po::value<long>(&targets)->default_value(1), "number of targets (if not bipartite)")
            ("edges,m", po::value<long>(&m), "number of edges (default: clique)")
            ("distr,d", po::value<int>(&distr)->default_value(0), "type of distribution for weights (0:u,1:g,2:e)")
            ("par1,l", po::value<long>(&min_weight)->default_value(0), "parameter1 for distribution")
            ("par2,u", po::value<long>(&max_weight)->default_value(1000), "parameter2 for distribution");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    po::notify(vm);


    Graph g;
    ofstream grfile(outf);
    switch(graph_type) {
        case 0:
            grfile << "# Type FullBipartite\n";
            switch(distr) {
                case 0:
                    grfile << "# Distribution Uniform\n";
                    grfile << "# MinCost " << min_weight << "\n";
                    grfile << "# MaxCost " << max_weight << "\n";
                    break;
                case 1:
                    grfile << "# Distribution Gaussian\n";
                    grfile << "# Mean " << min_weight << "\n";
                    grfile << "# Std " << max_weight << "\n";
                    break;
                case 2:
                    grfile << "# Distribution Exponential\n";
                    grfile << "# Lambda " << min_weight << "\n";
            }
            g.generate_full_bipartite_graph(n, min_weight, max_weight, distr);
            break;
        case 1:
            grfile << "# Type Clique\n";
            switch(distr) {
                case 0:
                    grfile << "# Distribution Uniform\n";
                    grfile << "# MinCost " << min_weight << "\n";
                    grfile << "# MaxCost " << max_weight << "\n";
                    break;
                case 1:
                    grfile << "# Distribution Gaussian\n";
                    grfile << "# Mean " << min_weight << "\n";
                    grfile << "# Std " << max_weight << "\n";
                    break;
                case 2:
                    grfile << "# Distribution Exponential\n";
                    grfile << "# Lambda " << min_weight << "\n";
            }
            g.generate_clique(n, min_weight, max_weight, distr);
            //redistribute sources and sinks
            grfile << "# Sources " << sources << "\n";
            grfile << "# Targets " << targets << "\n";
            vector<long> excesses(n,0);
            if (sources > targets) {
                cout << "Not enough sources" << endl;
                exit(1);
            }
            for (long i = 0; i < sources-1; i++) {
                excesses[i] = targets/sources;
            }
            excesses[sources-1] = targets/sources + targets % sources;
            for (long i = sources; i < sources+targets; i++) {
                excesses[i] = -1;
            }
            std::random_shuffle<vector<long>::iterator>(excesses.begin(),excesses.end());
            for (long i = 0; i < n; i++) {
                g.V[i].supply = excesses[i];
            }
            break;
    }
    grfile.close();
    if (graph_type == 0) g.save_graph_blossom(outfBlossom); // only for bipartite
    g.save_graph(outf);
    cout << "Done." << endl;
}
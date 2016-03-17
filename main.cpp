#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "Common.h"
#include "Graph.h"
#include "SIA.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, const char** argv) {
    std::cout << "Network Flows (c) Alvis Logins 2016" << std::endl;


    // check for parallelism
    // should be enabled at least for fast graph generation
#if defined(CILK)
    std::cout << "Cilk Enabled" << std::endl;
// intel cilk+
#elif defined(CILKP)
    std::cout << "CilkPlus Enabled" << std::endl;
// openmp
#elif defined(OPENMP)
    std::cout << "OpenMP enabled" << std::endl;
// c++
#else
    std::cout << "No parallelism enabled" << std::endl;
#endif
    //TODO openmp enabling

    // parsing parameters
    int algorithm;
    int graph_type;
    uintT size;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce help message")
            ("algorithm,a", po::value<int>(&algorithm)->default_value(0), "0:ALL, 1:SIA, 2:ASIA, 3:COSTSCALING")
            ("graph,g", po::value<int>(&graph_type)->default_value(0),"graph type 0:bipartite")
            ("output,o", po::value<string>(), "save generated graph")
            ("size,s", po::value<uintT>(&size)->default_value(100), "size of generated graph")
            ("input,i", po::value<string>(), "input graph")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    Graph g;
    switch(graph_type) {
        case 0:
            g.generate_full_bipartite_graph(size);
            break;
    }
    g.init_neighbors();
    g.sort_neighbors();
    g.test_graph_structure();
    g.test_sorting();

    switch(algorithm) {
        case 0:
            SIA SIAsolv(&g);
            SIAsolv.runOSIA();
            cout << "Total cost of SIA: " << SIAsolv.totalCost << endl;
            break;
    }

    cout << "Done." << endl;
    return 0;
}
//
// Created by alvis on 06.04.16.
//


#include <iostream>
#include <string>
#include <boost/program_options.hpp>
#include <random>
#include <chrono>
#include <fstream>

using namespace std;
namespace po = boost::program_options;

int main(int argc, const char** argv) {
    std::cout << ":: Spatial Graph Generator ::" << std::endl;

    // parsing parameters
    int distr;
    string outf;
    long n, excesses;
    int cl_num;
    bool bipartite;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce help message")
            ("output,o", po::value<string>(&outf)->required(), "save generated points")
            ("size,s", po::value<long>(&n)->required(), "number of nodes")
            ("clust,c", po::value<int>(&cl_num)->default_value(5), "number of clusters")
            ("excesses,e", po::value<long>(&excesses)->default_value(-1), "number of nodes with positive supply (-1 for bipartite)")
            ("distr,d", po::value<int>(&distr)->default_value(0), "type of distribution for weights (0:u,1:gauss,2:clust)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    po::notify(vm);


    long seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    uniform_real_distribution<double> uniGen(0, 1);
    normal_distribution<double> gausGen(0.5, 0.4);

    //distribute excesses uniformly
    if (excesses == -1) {
        if (n%2 != 0) {
            cout << "Size must be even" << endl;
            exit(1);
        }
        excesses = n/2;
    }
    vector<int> excessMap(n,0); //1 - excess, -1 - deficit
    for (long i = 0; i < excesses; i++) {
        excessMap[i] = 1;
        excessMap[i+excesses] = -1;
    }
    std::random_shuffle(excessMap.begin(), excessMap.end());

    ofstream grfile(outf);
    grfile << "# Type Random Points\n";
    switch (distr) {
        case 0:
            grfile << "# Distribution Uniform\n";
            for(long i = 0; i < n; i++) {
                grfile << uniGen(generator) << " " << uniGen(generator) << " " << excessMap[i] << endl;
            }
            break;
        case 1:
            grfile << "# Distribution Gaussian\n";
            for(long i = 0; i < n; i++) {
                double x, y;
                grfile << gausGen(generator) << " " << gausGen(generator) << " " << excessMap[i] << endl;
            }
            break;
        case 2:
            grfile << "# Distribution Clusters\n";
            uniform_real_distribution<double> centerGen(0.2, 0.8);
            uniform_real_distribution<double> stdGen(0.05, 0.3);
            for (int cl = 0; cl < cl_num; cl++) {
                normal_distribution<double> clGen(centerGen(generator), stdGen(generator));
                for(long i = 0; i < n/cl_num; i++) {
                    //each cluster is equal size
                    grfile << clGen(generator) << " " << clGen(generator) << " " << excessMap[cl*(n/cl_num)+i] << endl;
                }
                if (cl == cl_num -1) {
                    for(long i = 0; i < n % cl_num; i++) {
                        grfile << clGen(generator) << " " << clGen(generator) << " " << excessMap[(n - n % cl_num)+i] << endl;
                    }
                }
            }
            break;
    }
}
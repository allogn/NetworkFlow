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

// parsing parameters
int distr;
string outf;
long size, sources, targets;
int cl_num;
double clsize;
int stdistr;
bool bipartite;
int repeats;
int sdens;

void generate_graph(string output) {
    long seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    uniform_real_distribution<double> uniGen(0, 1);
    normal_distribution<double> gausGen(0.5, 0.4);

    //distribute excesses uniformly
    if (sources == -1) {
        if (size % 2 != 0) {
            cout << "Size must be even" << endl;
            exit(1);
        }
        sources = size/2;
        targets = size/2;
    }
    if (targets + sources > size) {
        cout << "Size is too small" << endl;
        exit(1);
    }
    //redistribute sources and sinks
    ofstream grfile(output);
    grfile << "# Sources " << sources << "\n";
    grfile << "# Targets " << targets << "\n";
    grfile << "# Supply " << sdens << "\n";
    vector<long> excessMap(size,0);
    if (sources > targets) {
        cout << "Not enough sources" << endl;
        exit(1);
    }

    switch (stdistr) {
        case 0:
            //fully random
            grfile << "# SourceDistr Random\n";
            for (long i = 0; i < sources-1; i++) {
                excessMap[i] = targets/sources*sdens;
            }
            excessMap[sources-1] = targets/sources*sdens + (targets % sources)*sdens;
            for (long i = sources; i < sources+targets; i++) {
                excessMap[i] = -sdens;
            }
            std::random_shuffle<vector<long>::iterator>(excessMap.begin(),excessMap.end());
            break;
        case 1:
            //one source one target - is not checked
            grfile << "# SourceDistr Pair\n";
            excessMap[1] = -sdens;
            excessMap[0] = sdens;
            break;
        case 2:
            //four sources, n targets
            grfile << "# SourceDistr Quartet\n";
            excessMap[0] = sdens*(targets/4);
            excessMap[1] = sdens*(targets/4);
            excessMap[2] = sdens*(targets/4);
            excessMap[3] = sdens*(targets/4 + targets % 4);
            for (long i = 4; i < 4+targets; i++) {
                excessMap[i] = -sdens;
            }
            std::random_shuffle<vector<long>::iterator>(excessMap.begin()+4,excessMap.end());
            break;
        case 3:
            //one cluster with sources in the center, clusters with targets around - is not checked
            grfile << "# SourceDistr CentralCluster\n";
            if (size / cl_num < sources || size / cl_num < targets / (cl_num - 1)) {
                cout << "Too small graph to distribute sources and targets" << endl;
                exit(1);
            }
            for (long i = 0; i < sources; i++) {
                excessMap[i] = sdens * targets/sources;
            }
            excessMap[sources-1] += sdens * (targets % sources);
            std::random_shuffle<vector<long>::iterator>(excessMap.begin(),excessMap.begin()+size/cl_num);

            for (long i = size/cl_num; i < targets + size/cl_num; i++) {
                excessMap[i] = -sdens;
            }
            std::random_shuffle<vector<long>::iterator>(excessMap.begin()+size/cl_num,excessMap.end());

            break;
        case 4:
            //one source in the center of each cluster, parameter -c is used for sources - is not checked
            grfile << "# SourceDistr ClusterCenter\n";
            for (long i = 0; i < cl_num-1; i++) {
                excessMap[size/cl_num*i] = sdens*targets/cl_num;
                for (long j = 0; j < targets/cl_num; j++) {
                    excessMap[size/cl_num*i + j + 1] = - sdens;
                }
                std::random_shuffle<vector<long>::iterator>(excessMap.begin()+size/cl_num*i+1,excessMap.begin()+size/cl_num*(i+1));
            }
            excessMap[size/cl_num*(cl_num-1)] = sdens*targets/cl_num + sdens*(targets % cl_num);
            for (long j = 0; j < targets/cl_num + targets % cl_num; j++) {
                excessMap[size/cl_num*(cl_num-1) + j + 1] = - sdens;
            }
            std::random_shuffle<vector<long>::iterator>(excessMap.begin()+size/cl_num*(cl_num - 1)+1,excessMap.end());
            break;
        case 5:
            //random distribution in each cluster that sums up to zero --is not checked
            grfile << "# SourceDistr ClusterEqual\n";
            for (long i = 0; i < cl_num-1; i++) {
                for (long j = size / cl_num * i; j < size / cl_num * i + sources/cl_num; j++) {
                    excessMap[j] = sdens*(targets/sources);
                }
                excessMap[size / cl_num * i + sources/cl_num] += sdens*(targets % sources);
                for (long j = size / cl_num * i + sources/cl_num; j < size / cl_num * i + sources/cl_num + targets/cl_num; j++) {
                    excessMap[j] = -sdens;
                }
                std::random_shuffle<vector<long>::iterator>(excessMap.begin()+size/cl_num*i,excessMap.begin()+size/cl_num*(i+1));
            }
            for (long j = size / cl_num * (cl_num - 1); j < size / cl_num * (cl_num-1) + sources/cl_num + sources%cl_num; j++) {
                excessMap[j] = sdens*(targets/sources);
            }
            excessMap[size / cl_num * (cl_num-1) + sources/cl_num] += sdens*(targets % sources);
            for (long j = size / cl_num * (cl_num-1) + sources/cl_num + sources%cl_num;
                 j < size / cl_num * (cl_num-1) + sources/cl_num + sources%cl_num + targets/cl_num + targets%cl_num; j++) {
                excessMap[j] = -sdens;
            }
            std::random_shuffle<vector<long>::iterator>(excessMap.begin()+size/cl_num*(cl_num-1),excessMap.end());
            break;
    }
    assert(std::accumulate(excessMap.begin(),excessMap.end(),0) == 0);
    long pos = 0;
    long neg = 0;
    for (long i = 0; i < excessMap.size(); i++) {
        if (excessMap[i] > 0) pos++;
        if (excessMap[i] < 0) neg++;
    }
    assert(pos == sources);
    assert(neg == targets);


    grfile << "# Type RandomPoints\n";
    switch (distr) {
        case 0:
            grfile << "# Distribution Uniform\n";
            grfile << size << endl;
            switch (stdistr) {
                case 0:
                    for(long i = 0; i < size; i++) {
                        grfile << uniGen(generator) << " " << uniGen(generator) << " " << excessMap[i] << endl;
                    }
                    break;
                case 1:
                    grfile << "0.2 0.2 " << excessMap[0] << endl;
                    grfile << "0.8 0.8 " << excessMap[1] << endl;
                    for(long i = 2; i < size; i++) {
                        grfile << uniGen(generator) << " " << uniGen(generator) << " " << excessMap[i] << endl;
                    }
                    break;
                case 2:
                    grfile << "0.2 0.2 " << excessMap[0] << endl;
                    grfile << "0.8 0.8 " << excessMap[1] << endl;
                    grfile << "0.2 0.8 " << excessMap[2] << endl;
                    grfile << "0.8 0.2 " << excessMap[3] << endl;
                    for(long i = 4; i < size; i++) {
                        grfile << uniGen(generator) << " " << uniGen(generator) << " " << excessMap[i] << endl;
                    }
                    break;
                case 3:
                    grfile << "0.5 0.5 " << excessMap[0] << endl;
                    for(long i = 1; i < size; i++) {
                        grfile << uniGen(generator) << " " << uniGen(generator) << " " << excessMap[i] << endl;
                    }
                    break;
                default:
                    cout << "Incorrect distribution combinations" << endl;
                    exit(1);
                    break;
            }
            break;
        case 1:
            grfile << "# Distribution Gaussian\n";
            grfile << size << endl;
            switch (stdistr) {
                case 0:
                    for(long i = 0; i < size; i++) {
                        double x, y;
                        grfile << gausGen(generator) << " " << gausGen(generator) << " " << excessMap[i] << endl;
                    }
                    break;
                case 1:
                    grfile << "0.2 0.2 " << excessMap[0] << endl;
                    grfile << "0.8 0.8 " << excessMap[1] << endl;
                    for(long i = 2; i < size; i++) {
                        double x, y;
                        grfile << gausGen(generator) << " " << gausGen(generator) << " " << excessMap[i] << endl;
                    }
                    break;
                case 2:
                    grfile << "0.2 0.2 " << excessMap[0] << endl;
                    grfile << "0.8 0.8 " << excessMap[1] << endl;
                    grfile << "0.2 0.8 " << excessMap[2] << endl;
                    grfile << "0.8 0.2 " << excessMap[3] << endl;
                    for(long i = 4; i < size; i++) {
                        double x, y;
                        grfile << gausGen(generator) << " " << gausGen(generator) << " " << excessMap[i] << endl;
                    }
                    break;
                case 3:
                    grfile << "0.5 0.5 " << excessMap[0] << endl;
                    for(long i = 1; i < size; i++) {
                        double x, y;
                        grfile << gausGen(generator) << " " << gausGen(generator) << " " << excessMap[i] << endl;
                    }
                    break;
                default:
                    cout << "Incorrect distribution combinations" << endl;
                    exit(1);
                    break;
            }
            break;
        case 2:
            grfile << "# Distribution Clusters\n";
            grfile << "# ClustersNum " << cl_num << "\n";
            grfile << size << endl;
            uniform_real_distribution<double> centerGen(0.2, 0.8);
            uniform_real_distribution<double> stdMultGen(0.5, 2);
            switch (stdistr) {
                case 0:
                    for (int cl = 0; cl < cl_num; cl++) {
                        normal_distribution<double> clGenX(centerGen(generator), stdMultGen(generator)*clsize);
                        normal_distribution<double> clGenY(centerGen(generator), stdMultGen(generator)*clsize);
                        for(long i = 0; i < size/cl_num; i++) {
                            //each cluster is equal size
                            grfile << clGenX(generator) << " " << clGenY(generator) << " " << excessMap[cl*(size/cl_num)+i] << endl;
                        }
                        if (cl == cl_num -1) {
                            for(long i = 0; i < size % cl_num; i++) {
                                grfile << clGenX(generator) << " " << clGenY(generator) << " " << excessMap[(size - size % cl_num)+i] << endl;
                            }
                        }
                    }
                    break;
                case 1:
                    //one s one t
                    grfile << "0.2 0.2 " << excessMap[0] << endl;
                    grfile << "0.8 0.8 " << excessMap[1] << endl;
                    for (int cl = 0; cl < cl_num; cl++) {
                        normal_distribution<double> clGen(centerGen(generator), stdMultGen(generator)*clsize);
                        for(long i = 0; i < (size-2)/cl_num; i++) {
                            //each cluster is equal size
                            grfile << clGen(generator) << " " << clGen(generator) << " " << excessMap[cl*(size/cl_num)+i+2] << endl;
                        }
                        if (cl == cl_num -1) {
                            for(long i = 0; i < (size-2) % cl_num; i++) {
                                grfile << clGen(generator) << " " << clGen(generator) << " " << excessMap[(size - size % cl_num)+i+2] << endl;
                            }
                        }
                    }
                    break;
                case 2:
                    grfile << "0.2 0.2 " << excessMap[0] << endl;
                    grfile << "0.8 0.8 " << excessMap[1] << endl;
                    grfile << "0.2 0.8 " << excessMap[2] << endl;
                    grfile << "0.8 0.2 " << excessMap[3] << endl;
                    for (int cl = 0; cl < cl_num; cl++) {
                        normal_distribution<double> clGen(centerGen(generator), stdMultGen(generator)*clsize);
                        for(long i = 0; i < (size-4)/cl_num; i++) {
                            //each cluster is equal size
                            grfile << clGen(generator) << " " << clGen(generator) << " " << excessMap[cl*(size/cl_num)+i+4] << endl;
                        }
                        if (cl == cl_num -1) {
                            for(long i = 0; i < (size-4) % cl_num; i++) {
                                grfile << clGen(generator) << " " << clGen(generator) << " " << excessMap[(size - size % cl_num)+i+4] << endl;
                            }
                        }
                    }
                    break;
                case 3:
                    {
                        normal_distribution<double> cclGenX(0.5, stdMultGen(generator)*clsize);
                        normal_distribution<double> cclGenY(0.5, stdMultGen(generator)*clsize);
                        for (long i = 0; i < size/cl_num; i++) {
                            grfile << cclGenX(generator) << " " << cclGenY(generator) << " " << excessMap[i] << endl;
                        }
                    }

                    for (int cl = 1; cl < cl_num; cl++) {
                        normal_distribution<double> clGenX(centerGen(generator), stdMultGen(generator)*clsize);
                        normal_distribution<double> clGenY(centerGen(generator), stdMultGen(generator)*clsize);
                        for(long i = 0; i < size/cl_num; i++) {
                            //each cluster is equal size
                            grfile << clGenX(generator) << " " << clGenY(generator) << " " << excessMap[cl*(size/cl_num)+i] << endl;
                        }
                        if (cl == cl_num-1) {
                            for(long i = 0; i < size % cl_num; i++) {
                                grfile << clGenX(generator) << " " << clGenY(generator) << " " << excessMap[(size - size % cl_num)+i] << endl;
                            }
                        }
                    }
                    break;
                case 4:
                    for (int cl = 0; cl < cl_num; cl++) {
                        double center_x = centerGen(generator);
                        double center_y = centerGen(generator);
                        double std = stdMultGen(generator)*clsize;
                        normal_distribution<double> clGenX(center_x, std);
                        normal_distribution<double> clGenY(center_y, std);
                        grfile << center_x << " " << center_y << " " << excessMap[cl*(size/cl_num)] << endl;
                        for(long i = 1; i < size/cl_num; i++) {
                            //each cluster is equal size
                            grfile << clGenX(generator) << " " << clGenY(generator) << " " << excessMap[cl*(size/cl_num)+i+1] << endl;
                        }
                        if (cl == cl_num -1) {
                            for(long i = 1; i < size % cl_num; i++) {
                                grfile << clGenX(generator) << " " << clGenY(generator) << " " << excessMap[(size - size % cl_num)+i+1] << endl;
                            }
                        }
                    }
                    break;
                case 5:
                    for (int cl = 0; cl < cl_num; cl++) {
                        normal_distribution<double> clGen(centerGen(generator), stdMultGen(generator)*clsize);
                        for(long i = 0; i < (size-4)/cl_num; i++) {
                            //each cluster is equal size
                            grfile << clGen(generator) << " " << clGen(generator) << " " << excessMap[cl*(size/cl_num)+i+4] << endl;
                        }
                        if (cl == cl_num -1) {
                            for(long i = 0; i < (size-4) % cl_num; i++) {
                                grfile << clGen(generator) << " " << clGen(generator) << " " << excessMap[(size - size % cl_num)+i+4] << endl;
                            }
                        }
                    }
                    break;
            }
            break;
    }
}

int main(int argc, const char** argv) {
    std::cout << ":: Spatial Graph Generator running ::" << std::endl;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce help message")
            ("output,o", po::value<string>(&outf)->required(), "output file (without extention, default: sgr)")
            ("instances,r", po::value<int>(&repeats)->default_value(1), "number of generated files")
            ("size,n", po::value<long>(&size)->required(), "number of nodes")
            ("clust,c", po::value<int>(&cl_num)->default_value(5), "number of clusters")
            ("stdistr,p", po::value<int>(&stdistr)->default_value(0), "type of source/target distribution\n"
                    "0 : fully random, amount of s/t set by -s -t params. s>t, supply(s) = s/t\n"
                    "1 : 1 source, 1 target in the opposite corners (almost corners)\n"
                    "2 : 4 sources in the opposite corners, random target in the middle\n"
                    "3 : some sources in the central cluster, -t targets around\n"
                    "4 : clusters: one source in the center of each cluster, -t targets in EACH cluster randomly around\n"
                    "5 : clusters: random -s -t in each cluster, supply of each cluster equal to zero\n")
            ("density", po::value<int>(&sdens)->default_value(1), "excess amount for each source")
            ("clsize", po::value<double>(&clsize)->default_value(0.2), "spatial size of each cluster ([0.,1.])")
            ("excesses,s", po::value<long>(&sources)->default_value(-1), "number of nodes with positive supply (-1 for bipartite)")
            ("deficits,t", po::value<long>(&targets)->default_value(-1), "number of nodes with negative supply (-1 for bipartite)")
            ("distr,d", po::value<int>(&distr)->default_value(0), "type of distribution of points (0:u,1:gauss,2:clust)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    po::notify(vm);
    for (int r = 1; r <= repeats; r++) {
        string name = outf;
        name = name.append("_").append(to_string(r)).append(".sgr");
        cout << name << endl;
        cout << "Generating " << r << "/" << repeats << " file..." << endl;
        generate_graph(name);
    }

    cout << "Done." << endl;
}
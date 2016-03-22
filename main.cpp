#include <iostream>
#include <boost/program_options.hpp>
#include <lemon/cost_scaling.h>
#include <time.h>
#include <sys/time.h>

#include "Common.h"
#include "Graph.h"
#include "SIA.h"
#include "CostScaling.h"
#include "LocalDominant.h"
#include "Lemon.h"
#include "SCS.h"
#include "TimerTool.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, const char** argv) {
    std::cout << "Network Flows (c) Alvis Logins 2016" << std::endl;
    Timer timer;
    double totalRunningTime = timer.getTime();

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
    int rounds;
    int experiment_id;
    string input_graph;
    string log_filename;
    uintT size;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce help message")
            ("algorithm,a", po::value<int>(&algorithm)->default_value(1),
                    "0 : SIA\n"
                    "1 : CostScaling\n"
                    "2 : LocalDominant\n"
                    "3 : Lemon Modified\n"
                    "4 : Simplified Cost Scaling (SCS)\n"
                    "5 : Original Lemon Cost Scaling\n")
            ("graph,g", po::value<int>(&graph_type)->default_value(0),"graph type 0:bipartite")
            ("output,o", po::value<string>(), "save generated graph")
            ("size,s", po::value<uintT>(&size)->default_value(100), "size of generated graph")
            ("input,i", po::value<string>(&input_graph)->required(), "input graph")
            ("rounds,r", po::value<int>(&rounds)->default_value(1), "rounds for an experiment")
            ("log,l", po::value<string>(&log_filename)->required(), "log file for timings")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    po::notify(vm);

    //init experiment ID
    struct timeval time;
    gettimeofday(&time,NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    experiment_id = rand();
    cout << "Experiment ID " << experiment_id << endl;

    //init log file
    ofstream logf(log_filename, ios::app);
    char buffer[30];
    time_t curtime;
    curtime=time.tv_sec;
    strftime(buffer,30,"%m-%d-%Y  %T.",localtime(&curtime));
    logf << experiment_id << "," << "DateTime" << "," << buffer << "\n";

    //save to log file algorithm name
    switch(algorithm) {
        case ALG_LEMON_MODIF:
            logf << experiment_id << ",Algorithm,Modified Lemon CostScaling" << "\n";
            break;
        case ALG_SCS:
            logf << experiment_id << ",Algorithm,Simplified CostScaling" << "\n";
            break;
        case ALG_SIA:
            logf << experiment_id << ",Algorithm,SIA" << "\n";
            break;
        case ALG_LEMON_ORIG:
            logf << experiment_id << ",Algorithm,Lemon CostScaling" << "\n";
            break;
        case ALG_LOCAL_DOMINANT:
            logf << experiment_id << ",Algorithm,Local Dominant" << "\n";
            break;
        default:
            cout << "Error: No such algorithm" << endl;
            exit(1);
    }
    logf.close(); //allow timers to write to the same file

    // check if file with graph exists
    if (FILE *file = fopen(input_graph.c_str(), "r")) {
        fclose(file);
    } else {
        cout << "File " << input_graph << " does not exist" << endl;
        exit(1);
    }

    //init Graph (Common for all algorithms)
    Graph g;
    g.load_graph(input_graph, log_filename, experiment_id);
    g.init_neighbors();

    for (int current_round = 0; current_round < rounds; current_round++ ) {
        cout << "== Round " << current_round << "/" << rounds << endl;
        switch(algorithm) {
            case 1:
                //prepare graph
                g.sort_neighbors();
                g.test_graph_structure();
                g.test_sorting();
                {
                    SIA SIAsolv(&g);
                    SIAsolv.runOSIA();
                    cout << "Total cost of SIA: " << SIAsolv.totalCost << endl;
                }
                break;
            case 2:
                g.add_all();
                {
                    CostScaling CostScalingSolv(&g);
                    CostScalingSolv.runCostScaling();
                    cout << "Total cost of CostScaling: " << CostScalingSolv.totalCost << endl;
                }
                break;
            case 3:
                g.add_all();
                {
                    LocalDominant LocalDominantSolv(&g);
                    LocalDominantSolv.runLocalDominant();
                    cout << "Total cost of LocalDominant: " << LocalDominantSolv.totalCost << endl;
                }
                break;
            case 4:
                g.add_all();
                {
                    Timer algtimer;
                    lemon::ListDigraph _graph;

                    lemon::ListDigraph::ArcMap<int> weight(_graph);
                    lemon::ListDigraph::ArcMap<int> flow(_graph);
                    lemon::ListDigraph::ArcMap<int> cap(_graph);
                    lemon::ListDigraph::ArcMap<int> lower(_graph);
                    lemon::ListDigraph::NodeMap<int> supply(_graph);

                    lemon::ListDigraph::Node* nodes;
                    nodes = (lemon::ListDigraph::Node*)malloc(sizeof(lemon::ListDigraph::Node)*g.n);
                    for (int i = 0; i < g.n; i++) {
                        nodes[i] = _graph.addNode();
                        supply[nodes[i]] = g.V[i].supply;
                    }

                    for (uintT i = 0; i < g.m; i++) {
                        lemon::ListDigraph::Arc e = _graph.addArc(nodes[g.E[i].fromid], nodes[g.E[i].toid]);
                        weight[e] = g.E[i].weight;
                        cap[e] = g.E[i].capacity;
                        lower[e] = g.E[i].lower;
                    }

                    lemon::ModifiedCostScaling<lemon::ListDigraph,int,int> cost_scaling(_graph);
                    cost_scaling.costMap(weight);
                    cost_scaling.upperMap(cap);
                    cost_scaling.lowerMap(lower);
                    cost_scaling.supplyMap(supply);
                    double total = timer.getTime();
                    cost_scaling.run();
                    algtimer.save_time("Total time", total);
                    cout << "Total Cost of Modified Lemon CostScaling: " << cost_scaling.totalCost() << endl;
                    cout << "Total Time of Modified Lemon CostScaling: " << algtimer.timings["Total time"].back() << endl;
                    algtimer.output(log_filename, experiment_id);
                    logf.open(log_filename, ios::app);
                    logf << experiment_id << ",Total cost," << cost_scaling.totalCost() << "\n";
                    logf.close();
                }
                break;
            case 5:
                g.add_all();
                {
                    SCS SCSsolv(g);
                    SCSsolv.runSCS();
                    cout << "Total cost of SCS: " << SCSsolv.totalCost << endl;
                }
                break;
            case 6:
                g.add_all();
                {
                    Timer algtimer;
                    lemon::ListDigraph _graph;

                    lemon::ListDigraph::ArcMap<int> weight(_graph);
                    lemon::ListDigraph::ArcMap<int> flow(_graph);
                    lemon::ListDigraph::ArcMap<int> cap(_graph);
                    lemon::ListDigraph::ArcMap<int> lower(_graph);
                    lemon::ListDigraph::NodeMap<int> supply(_graph);

                    lemon::ListDigraph::Node* nodes;
                    nodes = (lemon::ListDigraph::Node*)malloc(sizeof(lemon::ListDigraph::Node)*g.n);
                    for (int i = 0; i < g.n; i++) {
                        nodes[i] = _graph.addNode();
                        supply[nodes[i]] = g.V[i].supply;
                    }

                    for (uintT i = 0; i < g.m; i++) {
                        lemon::ListDigraph::Arc e = _graph.addArc(nodes[g.E[i].fromid], nodes[g.E[i].toid]);
                        weight[e] = g.E[i].weight;
                        cap[e] = g.E[i].capacity;
                        lower[e] = g.E[i].lower;
                    }

                    lemon::CostScaling<lemon::ListDigraph,int,int> cost_scaling(_graph);
                    cost_scaling.costMap(weight);
                    cost_scaling.upperMap(cap);
                    cost_scaling.lowerMap(lower);
                    cost_scaling.supplyMap(supply);
                    double total = timer.getTime();
                    cost_scaling.run();
                    algtimer.save_time("Total time", total);
                    cout << "Total Cost of Modified Lemon CostScaling: " << cost_scaling.totalCost() << endl;
                    cout << "Total Time of Modified Lemon CostScaling: " << algtimer.timings["Total time"].back() << endl;
                    algtimer.output(log_filename, experiment_id);
                    logf.open(log_filename, ios::app);
                    logf << experiment_id << ",Total cost," << cost_scaling.totalCost() << "\n";
                    logf.close();
                }
                break;
            default:
                cout << "Error: no such algorithm" << endl;
                exit(1);
        }
    }

    timer.save_time("Total execution time", totalRunningTime);
    timer.output(log_filename, experiment_id);
    cout << "Done." << endl;
    return 0;
}
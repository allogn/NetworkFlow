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
#include "LSCS.h"
#include "TimerTool.h"

#define LONG

using namespace std;
namespace po = boost::program_options;

int main(int argc, const char **argv) {
    std::cout << ":: Network Flows running ::" << std::endl;
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
    int rounds;
    int experiment_id;
    int param;
    double delta;
    string input_graph;
    string log_filename;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce help message")
            ("algorithm,a", po::value<int>(&algorithm)->required(),
             "0 : SIA\n"
                     "1 : CostScaling\n"
                     "2 : LocalDominant\n"
                     "3 : Lemon Modified\n"
                     "4 : Simplified Cost Scaling (SCS)\n"
                     "5 : Original Lemon Cost Scaling\n"
                     "6 : Lemon Simplified Cost Scaling (LSCS)\n"
                     "7 : Approximate SIA")
            ("input,i", po::value<string>(&input_graph)->required(), "input graph")
            ("param,p", po::value<int>(&param)->default_value(0), "parameter for an algorithm\n"
                    "For Lemon: 0 - no sorting, 1 - just sorting, 2 - sorting with pruning")
            ("rounds,r", po::value<int>(&rounds)->default_value(1), "rounds for an experiment")
            ("delta,d", po::value<double>(&delta)->default_value(0), "maximum connection distance (0:unbounded)")
            ("log,l", po::value<string>(&log_filename)->required(), "log file for timings");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    po::notify(vm);

    //init experiment ID
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    experiment_id = rand();
    cout << BASH_YELLOW << "Experiment ID " << experiment_id << BASH_NC << endl;

    //init log file
    ofstream logf(log_filename, ios::app);
    char buffer[30];
    time_t curtime;
    curtime = time.tv_sec;
    strftime(buffer, 30, "%m-%d-%Y  %T.", localtime(&curtime));
    logf << experiment_id << "," << "DateTime" << "," << buffer << "\n";
    logf << experiment_id << ",Input file," << input_graph << "\n";

    //save to log file algorithm name
    switch (algorithm) {
        case ALG_LEMON_MODIF:
            logf << experiment_id << ",Algorithm,DF-CSA" << "\n";
            break;
        case ALG_SCS:
            logf << experiment_id << ",Algorithm,S-CSA" << "\n";
            break;
        case ALG_SIA:
            logf << experiment_id << ",Algorithm,SIA" << "\n";
            break;
        case ALG_LEMON_ORIG:
            logf << experiment_id << ",Algorithm,Lemon CSA" << "\n";
            break;
        case ALG_LOCAL_DOMINANT:
            logf << experiment_id << ",Algorithm,Local Dominant" << "\n";
            break;
        case ALG_COST_SCALING:
            logf << experiment_id << ",Algorithm,CSA" << "\n";
            break;
        case ALG_LSCS:
            logf << experiment_id << ",Algorithm,S-DF-CSA" << "\n";
            break;
        case ALG_ASIA:
            logf << experiment_id << ",Algorithm,A-SIA" << "\n";
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
    if (input_graph.substr(input_graph.size()-3) == "sgr") {
        g.load_points(input_graph, log_filename, experiment_id);
        g.delta = delta;
    } else {
        g.load_graph(input_graph, log_filename, experiment_id);
    }
    g.init_neighbors();
    double result_cost;
    switch (algorithm) {
        case ALG_SIA: {
            if (!g.isSpatial) {
                //prepare graph
                cout << "Sorting edges..." << endl;
                g.sort_neighbors();
//                assert(g.test_graph_structure()); -- too long even for small graphs
//                assert(g.test_sorting());
            }

            //run SIA <rounds> times
            cout << "Staring SIA..." << endl;
            SIA SIAsolv(&g);
            for (int current_round = 1; current_round <= rounds; current_round++) {
                cout << "== Round " << current_round << "/" << rounds << " ==" << endl;
                g.clear_edge_list();
                SIAsolv.runOSIA();
                g.get_fill_status();
                cout << "Total cost of SIA: " << SIAsolv.totalCost << endl;
                cout << "Total time of SIA: " << SIAsolv.timer.timings["Total time"].back() << endl;
            }
            SIAsolv.save_profile_data(log_filename, experiment_id);
            SIAsolv.timer.output(log_filename, experiment_id);
        }
            break;
        case ALG_ASIA: {
            if (!g.isSpatial) {
                //prepare graph
                cout << "Sorting edges..." << endl;
                g.sort_neighbors();
//                assert(g.test_graph_structure()); -- too long even for small graphs
//                assert(g.test_sorting());
            }

            //run SIA <rounds> times
            cout << "Staring A-SIA..." << endl;
            SIA SIAsolv(&g, true);
            for (int current_round = 1; current_round <= rounds; current_round++) {
                cout << "== Round " << current_round << "/" << rounds << " ==" << endl;
                g.clear_edge_list();
                SIAsolv.runOSIA();
                cout << "Total cost of A-SIA: " << SIAsolv.totalCost << endl;
                cout << "Total time of A-SIA: " << SIAsolv.timer.timings["Total time"].back() << endl;
            }
            SIAsolv.save_profile_data(log_filename, experiment_id);
            SIAsolv.timer.output(log_filename, experiment_id);
        }
            break;
        case 1:
            for (int current_round = 1; current_round <= rounds; current_round++) {
                cout << "== Round " << current_round << "/" << rounds << " ==" << endl;
                CostScaling CostScalingSolv(&g);
                CostScalingSolv.runCostScaling();
                cout << "Total cost of CostScaling: " << CostScalingSolv.totalCost << endl;
            }
            break;
        case ALG_LOCAL_DOMINANT:
            for (int current_round = 1; current_round <= rounds; current_round++) {
                cout << "== Round " << current_round << "/" << rounds << " ==" << endl;
                LocalDominant LocalDominantSolv(&g);
                double totalTime = timer.getTime();
                LocalDominantSolv.runLocalDominant();
                timer.save_time("Total time", totalTime);

                long tc = LocalDominantSolv.totalCost();
                logf.open(log_filename, ios::app);
                logf << experiment_id << ",Total cost," << tc << "\n";
                logf.close();

                cout << "Total cost of LocalDominant: " << tc << endl;
                cout << "Total time of LocalDominant: " << timer.timings["Total time"].back() << endl;
            }
            break;
        case ALG_LEMON_MODIF:
            if (g.isSpatial) {
                g.fill_full_graph();
            }
            for (int current_round = 1; current_round <= rounds; current_round++) {
                cout << "== Round " << current_round << "/" << rounds << " ==" << endl;
                lemon::ListDigraph _graph;

                lemon::ListDigraph::ArcMap<int> weight(_graph);
                lemon::ListDigraph::ArcMap<int> flow(_graph);
                lemon::ListDigraph::ArcMap<int> cap(_graph);
                lemon::ListDigraph::ArcMap<int> lower(_graph);
                lemon::ListDigraph::NodeMap<int> supply(_graph);

                lemon::ListDigraph::Node *nodes;
                nodes = (lemon::ListDigraph::Node *) malloc(sizeof(lemon::ListDigraph::Node) * g.n);
                for (long i = g.n-1; i >=0; i--) {//inverse because ids in lemon for some reason assigned inversed
                    nodes[i] = _graph.addNode();
                    supply[nodes[i]] = g.V[i].supply;
                }

                for (long i = g.m-1; i >= 0; i--) {
                    lemon::ListDigraph::Arc e = _graph.addArc(nodes[g.E[i].fromid], nodes[g.E[i].toid]);
                    weight[e] = g.E[i].weight;
                    cap[e] = g.E[i].capacity;
                    lower[e] = g.E[i].lower;
                }

                lemon::ModifiedCostScaling<lemon::ListDigraph, int, int> cost_scaling(_graph, weight, param);
//                cost_scaling.costMap(weight);
                cost_scaling.upperMap(cap);
                cost_scaling.lowerMap(lower);
                cost_scaling.supplyMap(supply);
                cost_scaling.run();
                timer.save_time_total("Total time", cost_scaling.timer.timings["Total time"].back());
                cout << "Total Cost of Modified Lemon CostScaling: " << cost_scaling.totalCost() << endl;
                cout << "Total Time of Modified Lemon CostScaling: " << cost_scaling.timer.timings["Total time"].back() << endl;
                if (current_round == 1)
                    result_cost = cost_scaling.totalCost();
                else
                if (result_cost != cost_scaling.totalCost()) {
                    cout << "Error: two rounds gave different results" << endl;
                    exit(1);
                }
            }
            logf.open(log_filename, ios::app);
            logf << experiment_id << ",Total cost," << result_cost << "\n";
            switch(param) {
                case 0:
                    logf << experiment_id << ",Variation,No Sorting\n";
                    break;
                case 1:
                    logf << experiment_id << ",Variation,Just Sorting\n";
                    break;
                case 2:
                    logf << experiment_id << ",Variation,Sorting with pruning\n";
                    break;
            }
            logf.close();
            break;
        case ALG_SCS: {
            if (!g.isSpatial) {
                g.sort_neighbors();
            }
            SCS SCSsolv(g);
            if (rounds > 1) {
                cout << "rounds not working here yet" << endl; //todo
                exit(1);
            }
            for (int current_round = 1; current_round <= rounds; current_round++) {
                cout << "== Round " << current_round << "/" << rounds << " ==" << endl;
//                SCSsolv._graph.add_all(); // reset E lists --- this is for case if SIA implementation enabled
                SCSsolv.runSCS(param);
                g.get_fill_status();
                cout << "Total cost of SCS: " << SCSsolv.totalCost << endl;
                cout << "Total time of SCS: " << SCSsolv.timer.get_last_time("Total time") << endl;
            }
            SCSsolv.timer.output(log_filename, experiment_id);
            logf.open(log_filename, ios::app);
            logf << experiment_id << ",Total cost," << SCSsolv.totalCost << "\n";
            logf.close();
        }
            break;
        case ALG_LEMON_ORIG:
            if (g.isSpatial) {
                g.fill_full_graph();
            }
            for (int current_round = 1; current_round <= rounds; current_round++) {
                cout << "== Round " << current_round << "/" << rounds << " ==" << endl;
                lemon::ListDigraph _graph;

                lemon::ListDigraph::ArcMap<long> weight(_graph);
                lemon::ListDigraph::ArcMap<long> flow(_graph);
                lemon::ListDigraph::ArcMap<long> cap(_graph);
                lemon::ListDigraph::ArcMap<long> lower(_graph);
                lemon::ListDigraph::NodeMap<long> supply(_graph);

                lemon::ListDigraph::Node *nodes;
                nodes = (lemon::ListDigraph::Node *) malloc(sizeof(lemon::ListDigraph::Node) * g.n);
                for (long i = g.n-1; i >=0; i--) {
                    nodes[i] = _graph.addNode();
                    supply[nodes[i]] = g.V[i].supply;
                }

                for (long i = 0; i < g.m; i++) {
                    lemon::ListDigraph::Arc e = _graph.addArc(nodes[g.E[i].fromid], nodes[g.E[i].toid]);
                    weight[e] = g.E[i].weight;
                    cap[e] = g.E[i].capacity;
                    lower[e] = g.E[i].lower;
                }

                lemon::CostScaling<lemon::ListDigraph, long, long> cost_scaling(_graph);
                cost_scaling.costMap(weight);
                cost_scaling.upperMap(cap);
                cost_scaling.lowerMap(lower);
                cost_scaling.supplyMap(supply);
                double total = timer.getTime();
                int result = cost_scaling.run();
                if (result != cost_scaling.OPTIMAL)
                {
                    cout << "Problem unfeasible" << endl;
                    exit(1);
                }
                timer.save_time("Total time", total);
                cout << "Total Cost of Lemon CostScaling: " << cost_scaling.totalCost() << endl;
                cout << "Total Time of Lemon CostScaling: " << timer.timings["Total time"].back() << endl;
                if (current_round == 1)
                    result_cost = cost_scaling.totalCost();
                else
                    if (result_cost != cost_scaling.totalCost()) {
                        cout << "Error: two rounds gave different results" << endl;
                        exit(1);
                    }
            }
            logf.open(log_filename, ios::app);
            logf << experiment_id << ",Total cost," << result_cost << "\n";
            logf.close();
            break;
        case ALG_LSCS: {
            if (!g.isSpatial) {
                //prepare graph
                double sortingTime = timer.getTime();
                cout << "Sorting edges..." << endl;
                g.sort_neighbors();
                timer.save_time("Sorting", sortingTime);
            }
//            g.fill_full_graph();
            LSCS LSCSsolv(g);
            for (int current_round = 1; current_round <= rounds; current_round++) {
                cout << "== Round " << current_round << "/" << rounds << " ==" << endl;
                LSCSsolv.runLSCS();
                cout << "Total cost of LSCS: " << LSCSsolv.totalCost << endl;
                cout << "Total time of LSCS: " << LSCSsolv.timer.get_last_time("Total time") << endl;
            }
            //LSCSsolv.save_profile_data(log_filename, experiment_id);
            LSCSsolv.timer.output(log_filename, experiment_id);
            logf.open(log_filename, ios::app);
            logf << experiment_id << ",Total cost," << LSCSsolv.totalCost << "\n";
            logf.close();
        }
            break;
        default:
            cout << "Error: no such algorithm" << endl;
            exit(1);
    }

    g.timer.output(log_filename,experiment_id);
    timer.save_time("Total execution time", totalRunningTime);
    timer.output(log_filename, experiment_id);
    cout << "Done." << endl;
    return 0;
}
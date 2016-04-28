//
// Created by alvis on 16.03.16.
//

#include <iostream>
#include <fstream>
#include <assert.h>
#include <boost/program_options.hpp>
#include <lemon/cost_scaling.h>

#include "Common.h"
#include "Graph.h"
#include "SIA.h"
#include "CostScaling.h"
#include "Lemon.h"
#include "SCS.h"
#include "LSCS.h"

//#define NDEBUG

namespace po = boost::program_options;

void test_utility() {
    Graph g;

    g.generate_full_bipartite_graph(100, 0, 1000);
    g.test_graph_structure();
    g.init_neighbors();
    g.sort_neighbors();
    g.test_sorting();

    g.test_save_load();
}

int main(int argc, const char** argv) {

    int algorithm;
    po::options_description desc("Allowed options");
    desc.add_options()
            ("help,h", "produce help message")
            ("algorithm,a", po::value<int>(&algorithm)->default_value(-1),
                    "Default: Test utility functions\n"
                    "0 : SIA\n"
                    "1 : CostScaling\n"
                    "2 : LocalDominant\n"
                    "3 : Lemon Modified\n"
                    "4 : Simplified Cost Scaling (SCS)\n"
                    "5 : Original Lemon Cost Scaling\n"
                    "6 : Lemon Simplified Cost Scaling (LSCS)\n")
            ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }
    po::notify(vm);

    if (algorithm == -1) {
        test_utility();
        //todo not working load-save
    }

    // our star
    Graph g;


    /*
     * Check bipartite and compare with blossom
     */

    cout << BASH_YELLOW << "Compare results with answer sheet..." << BASH_NC << endl;
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
        string path = "../../data/tests/bipartite/";
        path.append(graphfile);
        path[path.size()-1] = 'r';
        path[path.size()-2] = 'g';

        if (algorithm == ALG_ASIA) {
            g.clear_graph();
            g.load_graph(path);
            g.init_neighbors();
            g.sort_neighbors();
            SIA SIAsolv(&g, true);
            SIAsolv.runOSIA();
            if (SIAsolv.totalCost != answer) {
                cout << "SIA wrong answer: " << SIAsolv.totalCost << " correct: " << answer << endl;
                exit(1);
            };
        }

        if (algorithm == ALG_COST_SCALING) {
            g.clear_graph();
            g.load_graph(path);
            g.init_neighbors();
            CostScaling CostScalingSolv(&g);
            CostScalingSolv.runCostScaling();
            if (CostScalingSolv.totalCost != answer) {
                cout << "CostScaling wrong answer: " << CostScalingSolv.totalCost << " correct: " << answer << endl;
                exit(1);
            };
        }

        if (algorithm == ALG_SCS) {
            g.clear_graph();
            g.load_graph(path);
            g.init_neighbors();
            g.sort_neighbors();
            SCS SCSsolv(g);
            SCSsolv.runSCS(0);
            if (SCSsolv.totalCost != answer) {
                cout << "Simplified CostScaling wrong answer: " << SCSsolv.totalCost << " correct: " << answer << endl;
                exit(1);
            };
        }


        if (algorithm == ALG_LSCS) {
            g.clear_graph();
            g.load_graph(path);
            g.init_neighbors();
            g.sort_neighbors();
            LSCS LSCSsolv(g);
            LSCSsolv.runLSCS();
            if (LSCSsolv.totalCost != answer) {
                cout << "Lemon Simplified CostScaling wrong answer: " << LSCSsolv.totalCost << " correct: " << answer << endl;
                exit(1);
            };
        }

        if (algorithm == ALG_LEMON_MODIF) {
            g.clear_graph();
            g.load_graph(path);
            g.init_neighbors();
            lemon::ListDigraph _graph;
            lemon::ListDigraph::ArcMap<int> weight(_graph);
            lemon::ListDigraph::ArcMap<int> flow(_graph);
            lemon::ListDigraph::ArcMap<int> cap(_graph);
            lemon::ListDigraph::ArcMap<int> lower(_graph);
            lemon::ListDigraph::NodeMap<int> supply(_graph);

            lemon::ListDigraph::Node* nodes;
            nodes = (lemon::ListDigraph::Node*)malloc(sizeof(lemon::ListDigraph::Node)*g.n);
            for (long i = 0; i < g.n; i++) {
                nodes[i] = _graph.addNode();
                supply[nodes[i]] = g.V[i].supply;
            }

            for (long i = 0; i < g.m; i++) {
                lemon::ListDigraph::Arc e = _graph.addArc(nodes[g.E[i].fromid], nodes[g.E[i].toid]);
                weight[e] = g.E[i].weight;
                cap[e] = g.E[i].capacity;
                lower[e] = g.E[i].lower;
            }

            lemon::ModifiedCostScaling<lemon::ListDigraph,int,int> cost_scaling(_graph, weight);
            cost_scaling.costMap(weight);
            cost_scaling.upperMap(cap);
            cost_scaling.lowerMap(lower);
            cost_scaling.supplyMap(supply);
            cost_scaling.run();
            if (cost_scaling.totalCost() != answer) {
                cout << "Modified Lemon CostScaling wrong answer: " << cost_scaling.totalCost() << " correct: " << answer << endl;
                exit(1);
            }
        }
    }


    /*
     * Check spatial random generated points
     */

    cout << BASH_YELLOW << "Spatial points checks..." << BASH_NC << endl;
    if (algorithm == ALG_LEMON_MODIF || algorithm == ALG_SCS) {
        int total = 100;
        cout << "Checking Spatial data " << total << " times" << endl;
        for (int i = 0; i < total; i++) {
            cout << i << "/" << total << endl;
            g.clear_graph();
            g.generate_random_points(200,i%100,i%100);
            g.init_neighbors();
            g.fill_full_graph();
            assert(g.m == g.n*(g.n-1));
            lemon::ListDigraph _graph;

            lemon::ListDigraph::ArcMap<int> weight(_graph);
            lemon::ListDigraph::ArcMap<int> flow(_graph);
            lemon::ListDigraph::ArcMap<int> cap(_graph);
            lemon::ListDigraph::ArcMap<int> lower(_graph);
            lemon::ListDigraph::NodeMap<int> supply(_graph);

            lemon::ListDigraph::Node* nodes;
            nodes = (lemon::ListDigraph::Node*)malloc(sizeof(lemon::ListDigraph::Node)*g.n);
            for (long i = 0; i < g.n; i++) {
                nodes[i] = _graph.addNode();
                supply[nodes[i]] = g.V[i].supply;
            }

            for (long i = 0; i < g.m; i++) {
                lemon::ListDigraph::Arc e = _graph.addArc(nodes[g.E[i].fromid], nodes[g.E[i].toid]);
                weight[e] = g.E[i].weight;
                cap[e] = g.E[i].capacity;
                lower[e] = g.E[i].lower;
            }

            lemon::ModifiedCostScaling<lemon::ListDigraph,int,int> MLCSsolv(_graph, weight);
            MLCSsolv.upperMap(cap);
            MLCSsolv.lowerMap(lower);
            MLCSsolv.supplyMap(supply);
            MLCSsolv.run();
            long long totalCost = MLCSsolv.totalCost();
            lemon::CostScaling<lemon::ListDigraph,int,int> LCSsolv(_graph);
            LCSsolv.costMap(weight);
            LCSsolv.upperMap(cap);
            LCSsolv.lowerMap(lower);
            LCSsolv.supplyMap(supply);
            LCSsolv.run();
            long long totalCostCorrect = LCSsolv.totalCost();

            if (totalCost != totalCostCorrect) {
                cout << "Wrong result: MLCS=" << totalCost << ", LCS=" << totalCostCorrect << endl;
                g.save_graph("temp.gr");
                cout << "Graph saved in temp.gr." << endl;
                exit(1);
            }

            //reinitialize
            g.E.clear();
            for (long j = 0; j < g.n; j++) {
                g.QryIt[j] = g.rtree.qbegin(bgi::nearest(g.coords[j], (unsigned int)g.n));
                ++g.QryIt[j];//skip the node itself
            }
            //add edges incrementally
            SCS SCSsolv(g);
            double time = Timer::getTime();
            SCSsolv.runSCS(0);
            time = -time + Timer::getTime();
            cout << "Time: " << time << endl;
            if (SCSsolv.totalCost != totalCostCorrect) {
                cout << "Simplified CostScaling wrong answer: " << SCSsolv.totalCost << " correct: " << totalCostCorrect << endl;
                exit(1);
            };
        }
    }



    /*
     * Generate a lot of random small graphs and compare results
     */
    cout << BASH_YELLOW << "Random graph checks..." << BASH_NC << endl;
    if (algorithm == ALG_LEMON_MODIF || algorithm == ALG_SCS) {
        int total = 100;
        cout << "Checking Lemon Modified " << total << " times" << endl;
        for (int i = 0; i < total; i++) {
            cout << i << "/" << total << endl;
            g.clear_graph();
            g.generate_full_bipartite_graph(100,0,1000);
            g.init_neighbors();

            lemon::ListDigraph _graph;

            lemon::ListDigraph::ArcMap<int> weight(_graph);
            lemon::ListDigraph::ArcMap<int> flow(_graph);
            lemon::ListDigraph::ArcMap<int> cap(_graph);
            lemon::ListDigraph::ArcMap<int> lower(_graph);
            lemon::ListDigraph::NodeMap<int> supply(_graph);

            lemon::ListDigraph::Node* nodes;
            nodes = (lemon::ListDigraph::Node*)malloc(sizeof(lemon::ListDigraph::Node)*g.n);
            for (long i = 0; i < g.n; i++) {
                nodes[i] = _graph.addNode();
                supply[nodes[i]] = g.V[i].supply;
            }

            for (long i = 0; i < g.m; i++) {
                lemon::ListDigraph::Arc e = _graph.addArc(nodes[g.E[i].fromid], nodes[g.E[i].toid]);
                weight[e] = g.E[i].weight;
                cap[e] = g.E[i].capacity;
                lower[e] = g.E[i].lower;
            }

            lemon::ModifiedCostScaling<lemon::ListDigraph,int,int> MLCSsolv(_graph, weight);
            MLCSsolv.upperMap(cap);
            MLCSsolv.lowerMap(lower);
            MLCSsolv.supplyMap(supply);
            MLCSsolv.run();
            long long totalCost = MLCSsolv.totalCost();
            lemon::CostScaling<lemon::ListDigraph,int,int> LCSsolv(_graph);
            LCSsolv.costMap(weight);
            LCSsolv.upperMap(cap);
            LCSsolv.lowerMap(lower);
            LCSsolv.supplyMap(supply);
            LCSsolv.run();
            long long totalCostCorrect = LCSsolv.totalCost();

            if (totalCost != totalCostCorrect) {
                cout << "Wrong result: MLCS=" << totalCost << ", LCS=" << totalCostCorrect << endl;
                g.save_graph("temp.gr");
                cout << "Graph saved in temp.gr." << endl;
                exit(1);
            }

            g.sort_neighbors();
            SCS SCSsolv(g);
            SCSsolv.runSCS(0);
            if (SCSsolv.totalCost != totalCostCorrect) {
                cout << "Simplified CostScaling wrong answer: " << SCSsolv.totalCost << " correct: " << totalCostCorrect << endl;
                exit(1);
            };
        }
    }



    cout << "Testing done" << endl;
}
//
// Created by alvis on 16.03.16.
//

#ifndef NETWORKFLOW_GRAPH_H
#define NETWORKFLOW_GRAPH_H

#include <vector>
#include <iostream>
#include <sys/time.h>
#include <fstream>

#include "Common.h"
#include "parallel.h"

using namespace std;

class Edge {
public:
    uintT fromid;
    uintT toid;
    intT weight;
    uintT capacity;
    uintT lower;
    uintT ID;
};

class Vertex {
public:
    uintT ID;
    intT supply;
    std::vector<uintT> E;

    Vertex(uintT _ID, intT _supply) {
        E.clear();
        ID = _ID;
        supply = _supply;
    }
};

class Graph {
public:
    uintT n;
    uintT m;

    std::vector<Vertex> V;
    std::vector<Edge> E;

    // additional data for incremental neighbour addition
    vector<vector<uintT>> fullE;
    vector<uintT> cursor;

    void clear_graph();

    inline uintT get_pair(uintT edgeId, uintT nodeId) {
        return (E[edgeId].fromid == nodeId)?E[edgeId].toid:E[edgeId].fromid;
    }

    // initialization for simple graph includes sorting and storing sorted edges
    // for spatial data - building NN data structure
    void init_neighbors() {
        fullE.resize(n);
        parallel_for(uintT i = 0; i < n; i++) {
            fullE[i].clear();
            // add only forward edges!
            for(uintT j = 0; j < m; j++) {
                if (E[j].fromid == i) {
                    fullE[i].push_back(j);
                }
            }
        }
        // pointer to the next element for insertion to E_sub
        cursor.resize(n,0);
    };
    void myqsort(std::vector<uintT>& edgelist, int start, int end)
    {
        if (start >= end) return;

        uintT temp;
        int pivot_ind = start;
        int pivot_val = E[edgelist[start]].weight;
        int i = start+1;
        while (i <= end)
        {
            if (E[edgelist[i]].weight < pivot_val)
            {
                temp = edgelist[i];
                edgelist[i] = edgelist[pivot_ind+1];
                edgelist[pivot_ind+1] = edgelist[pivot_ind];
                edgelist[pivot_ind++] = temp;
            }
            i++;
        }
        myqsort(edgelist,start,pivot_ind-1);
        myqsort(edgelist,pivot_ind+1,end);
        return;
    }
    void sort_neighbors() {
        cout << "Sorting edges..." << endl;
        for (int i = 0; i < n; i++)
        {
            myqsort(fullE[i], 0, fullE[i].size()-1);
        }
    }
    int get_next_neighbour(uintT nodeId) {
        if (cursor[nodeId] >= V[nodeId].E.size()) {
            cout << "Error: cursor is out of range" << endl;
            exit(1);
        }
        return fullE[nodeId][cursor[nodeId]++];
    };

    // graph generators
    void generate_full_bipartite_graph(uintT size, uintT min_weight = 0, uintT max_weight = UINT_T_MAX);

    void print_graph();
    void save_graph(string& filename);
    void load_graph(string& filename);
    void save_graph_info(string& filename, int experiment_id);

    // operators
    bool operator==(const Graph &other) const {
        // compare E, V, fullE
        for (uintT i = 0; i < n; i++) {
            if (V[i].E != other.V[i].E ||
                    V[i].ID != other.V[i].ID ||
                    V[i].supply != other.V[i].supply)
                return false;
            if (fullE[i] != other.fullE[i])
                return false;
        }
        for (uintT i = 0; i < m; i++) {
            if (E[i].ID != other.E[i].ID ||
                    E[i].fromid != other.E[i].fromid ||
                    E[i].toid != other.E[i].toid ||
                    E[i].capacity != other.E[i].capacity ||
                    E[i].lower != other.E[i].lower)
                return false;
        }
        return true;
    }

    // Unit Tests
    bool test_graph_structure();
    bool test_sorting();
    bool test_save_load();
};


#endif //NETWORKFLOW_GRAPH_H

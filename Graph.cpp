//
// Created by alvis on 16.03.16.
//

#include "Graph.h"

void Graph::clear_graph() {
    V.clear();
    E.clear();
    m = 0;
    n = 0;
    fullE.clear();
    cursor.clear();
}

void Graph::generate_full_bipartite_graph(uintT size, uintT min_weight, uintT max_weight) {
    std::cout << "Bipartite graph generation..." << std::endl;
    clear_graph();

    if (size % 2 != 0) {
        std::cout << "Error: size for bipartite graph must be even" << std::endl;
        exit(1);
    }
    n = size;
    m = size*size/4;
    E.reserve(m);
    V.reserve(n);

    // reset random generator
    struct timeval time;
    gettimeofday(&time,NULL);

    // microsecond has 1 000 000
    // Assuming you did not need quite that accuracy
    // Also do not assume the system clock has that accuracy.
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));

    parallel_for(uintT i = 0; i < n/2; i++) {
        for (uintT j = n/2; j < n; j++) {
            Edge e;
            e.capacity = 1;
            e.fromid = i;
            e.toid = j;
            e.ID = E.size();
            e.lower = 0;
            e.weight = (min_weight + rand() % (UINT_T_MAX - min_weight)) % max_weight;
            E.push_back(e);
        }
    }

    parallel_for(uintT i = 0; i < n; i++) {
        Vertex v(i,0);
        V.push_back(v);
    }
}

void Graph::print_graph() {
    std::cout << "Total " << n << " nodes and " << m << " edges" << std::endl;
    std::cout << "Nodes: " << std::endl;
    for (uintT i = 0; i < n; i++) {
        std::cout << "Node " << i << " edges: ";
        for (uintT j = 0; j < V[i].E.size(); j++) {
            std::cout << E[V[i].E[j]].ID << ",";
        }
        std::cout << std::endl;
    }
    std::cout << "Edges: " << std::endl;
    for (uintT i = 0; i < m; i++) {
        Edge e = E[i];
        std::cout << i << " (" << e.fromid << "->" << e.toid << ", weight " << e.weight << ")" << std::endl;
    }
}

bool Graph::test_graph_structure() {
    std::cout << "Testing graph structure..." << std::endl;
    assert(n == V.size());
    assert(m == E.size());

    assert(n > 0);
    assert(m > 0);

    // test that edges are unique (no two edges between two vertices)
    parallel_for(uintT i = 0; i < m; i++) {
        uintT curId = E[i].ID;
        assert(E[i].ID == i);
        uintT curFromId = E[i].fromid;
        for (uintT j = 0; j < m; j++) {
            if (E[j].fromid == curFromId) {
                assert(E[j].toid != E[i].fromid || E[j].ID == curId);
            }
            if (E[j].ID == curId)
                assert(j == i);
        }
    }

    // test if nodes are unique
    parallel_for(uintT i = 0; i < n; i++) {
        uintT curId = V[i].ID;
        assert(V[i].ID == i);
        for (uintT j = 0; j < n; j++) {
            if (V[j].ID == curId)
                assert(j == i);
        }
    }

    // test if every edge in the lists of edges is unique (but may not exist)
    bool exist[m];
    parallel_for(uintT i = 0; i < m; i++) exist[i] = false;
    for(uintT i = 0; i < n; i++) {
        for(uintT j = 0; j < V[i].E.size(); j++) {
            assert(V[i].E[j] < m);
            assert(exist[V[i].E[j]] == false);
            exist[V[i].E[j]] = true;
        }
    }
    return true;
}

bool Graph::test_sorting() {
    cout << "Testing sorting..." << endl;

    // test that every edge exists exactly once
    assert(fullE.size() == n);
    assert(fullE.size() > 0);
    bool exists[m];
    for (uintT i = 0; i < m; i++) exists[i] = false;
    for (uintT i = 0; i < n; i++) {
        for (uintT j = 0; j < fullE[i].size(); j++) {
            assert(fullE[i][j] < m);
            assert(exists[fullE[i][j]] == false);
            exists[fullE[i][j]] = true;
        }
    }
    for (uintT i = 0; i < m; i++) {
        assert(exists[i] == true);
    }

    for (uintT i = 0; i < n; i++) {
        if (fullE[i].size() > 0) {
            for (uintT j = 0; j < fullE[i].size() - 1; j++) {
                assert(j < fullE[i].size());
                assert(E[fullE[i][j]].weight <= E[fullE[i][j+1]].weight);
            }
        }
    }
    return true;
}

/*
 * My Graph Format: .gr
 * nodes edges
 * supply for nodeId1
 * supply for nodeId2
 * ...
 * (from)nodeId1 (to)nodeId2 weight capacity lower(bound)
 */
void Graph::save_graph(string& filename) {
    cout << "Saving graph to " << filename << "..." << endl;
    ofstream outf(filename);

    outf << n << " " << m << "\n";
    for (uintT i = 0; i < n; i++) {
        outf << V[i].supply << "\n";
    }
    for (uintT i = 0; i < m; i++) {
        outf << E[i].fromid << " " << E[i].toid << " " << E[i].weight << " " << E[i].capacity << " " << E[i].lower << "\n";
    }
    outf.close();
}

void Graph::load_graph(string& filename) {
    cout << "Loading graph from " << filename << "..." << endl;
    ifstream infile(filename);

    clear_graph();
    infile >> n >> m;
    intT supply;
    for (uintT i = 0; i < n; i++) {
        infile >> supply;
        V.push_back(Vertex(i, supply));
    }
    for (uintT i = 0; i < m; i++) {
        Edge e;
        infile >> e.fromid >> e.toid >> e.weight >> e.capacity >> e.lower;
    }
    infile.close();
}

bool Graph::test_save_load() {
    cout << "Testing save/load..." << endl;
    //save, then load and check if the graph is equal

    Graph g;
    g = *this;

    string fname = "test_save_load.gr";
    save_graph(fname);
    load_graph(fname);
    assert(g == *this);
    return true;
}

void Graph::save_graph_info(string &filename, int experiment_id) {
    cout << "Saving data to log file (ID " << experiment_id << endl;

    ofstream outf(filename);
    outf << experiment_id << ",nodes," << n << "\n";
    outf << experiment_id << ",edges," << m << "\n";
    outf.close();
}

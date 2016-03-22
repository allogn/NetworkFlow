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

/*
 * @input
 * distr - type of weight distribution
 * 0 : uniform (param1: min, param2: max)
 * 1 : gaussian (param1: mean, param2: std)
 * 2 : exponential (param1: lambda)
 */
void Graph::generate_full_bipartite_graph(uintT size, uintT param1, uintT param2, int distr) {
    clear_graph();

    if (size % 2 != 0) {
        std::cout << "Error: size for bipartite graph must be even" << std::endl;
        exit(1);
    }
    n = size;
    m = size*size/4;
    E.reserve(m);
    V.reserve(n);

    long seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    uniform_int_distribution<uintT> uniGen(param1, param2);
    normal_distribution<double> gausGen(param1, param2);
    exponential_distribution<double> expGen((double)param1/100.);

    parallel_for(uintT i = 0; i < n/2; i++) {
        for (uintT j = n/2; j < n; j++) {
            Edge e(E.size());
            e.capacity = 1;
            e.fromid = i;
            e.toid = j;
            e.lower = 0;
            switch(distr) {
                case 0:
                    e.weight = uniGen(generator);
                    break;
                case 1:
                    e.weight = (uintT)abs(gausGen(generator));
                    break;
                case 2:
                    e.weight = (uintT)abs(expGen(generator));
                    break;
                default:
                    cout << "Error: incorrect distribution" << endl;
                    exit(1);
            }
            E.push_back(e);
        }
    }

    parallel_for(uintT i = 0; i < n; i++) {
        if (i < n/2) {
            Vertex v(i,1);
            V.push_back(v);
        } else {
            Vertex v(i,-1);
            V.push_back(v);
        }
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
 * # Parameter Value
 * # Parameter2 1023
 * nodes edges
 * supply for nodeId1
 * supply for nodeId2
 * ...
 * (from)nodeId1 (to)nodeId2 weight capacity lower(bound)
 */
void Graph::save_graph(string filename) {
    cout << "Saving graph to " << filename << "..." << endl;
    ofstream outf(filename,ios::app); // appending because graph annotation and parameters can be before

    outf << n << " " << m << "\n";
    for (uintT i = 0; i < n; i++) {
        outf << V[i].supply << "\n";
    }
    for (uintT i = 0; i < m; i++) {
        outf << E[i].fromid << " " << E[i].toid << " " << E[i].weight << " " << E[i].capacity << " " << E[i].lower << "\n";
    }
    outf.close();
}

void Graph::load_graph(string filename, string log_filename = "", int experiment_id = 0) {
    cout << "Loading graph from " << filename << "..." << endl;
    if (FILE *file = fopen(filename.c_str(), "r")) {
        fclose(file);
    } else {
        cout << "File " << filename << " does not exist" << endl;
        exit(1);
    }
    ifstream infile(filename);

    if (log_filename != "") {
        if (FILE *file = fopen(log_filename.c_str(), "r")) {
            fclose(file);
        } else {
            cout << "File " << log_filename << " does not exist" << endl;
            exit(1);
        }
    }
    ofstream log(log_filename);

    //parse comments in the beginning of the file with parameters of data
    string line;
    while (infile >> line) {
        if (line[0] != '#') {
            infile.close();
            infile.open(filename);
            break;
        }
        char parameter[1000];
        char value[1000];
        sscanf(line.c_str(),"# %s %s",parameter,value);
        if (log_filename != "") log << experiment_id << "," << parameter << "," << value << "\n";
    };

    clear_graph();
    infile >> n >> m;
    if (log_filename != "") log << experiment_id << ",Edges," << m << "\n";
    if (log_filename != "") log << experiment_id << ",Nodes," << n << "\n";
    intT supply;
    for (uintT i = 0; i < n; i++) {
        infile >> supply;
        V.push_back(Vertex(i, supply));
    }
    for (uintT i = 0; i < m; i++) {
        Edge e(i);
        infile >> e.fromid >> e.toid >> e.weight >> e.capacity >> e.lower;
        E.push_back(e);
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
    load_graph(fname, "test.log", 0);
    assert(g == *this);
    return true;
}

void Graph::save_graph_info(string filename, int experiment_id) {
    cout << "Saving data to log file (ID " << experiment_id << endl;

    ofstream outf(filename);
    outf << experiment_id << ",nodes," << n << "\n";
    outf << experiment_id << ",edges," << m << "\n";
    outf.close();
}

void Graph::save_graph_blossom(string filename) {
    cout << "Saving blossom graph to " << filename << "..." << endl;
    ofstream outf(filename);
    outf << n << " " << m << "\n";
    for (int i = 0; i < m; i++)
    {
        outf << E[i].fromid << " " << E[i].toid << " " << E[i].weight << "\n";
    }
    outf.close();
}

void Graph::load_lgf_graph(string filename) {
    //TODO
//    ListDigraph g;
//    ListDigraph::ArcMap<int> weight(g);
//    ListDigraph::ArcMap<int> flow(g);
//    ListDigraph::NodeMap<int> potential(g);
//
//    digraphReader(g, argv[1])
//            .arcMap("weight", weight)
//            .run();
//
//    ListDigraph::ArcMap<int> cap(g);
//    for (ListDigraph::ArcIt it(g); it != INVALID; ++it) {
//        cap[it] = 1;
//    }
//
//    //set lower values on flow
//    //by default: first node is an excess, last - a deficit
//
//    ListDigraph::Node source = g.nodeFromId(0);
//    ListDigraph::Node target = g.nodeFromId(countNodes(g)-1);
//
//    ListDigraph::ArcMap<int> lower(g);
//    for (ListDigraph::ArcIt it(g); it != INVALID; ++it) {
//        lower[it] = 0;
//    }
}

void Graph::load_adj_graph(string filename) {
    char fname[1000];
    strcpy(fname, filename.c_str());

    _seq<char> S = readStringFromFile(fname);
    words W = stringToWords(S.A, S.n);
    if (W.Strings[0] != (string) "AdjacencyGraph") {
        cout << "Bad input file" << endl;
        abort();
    }

    long len = W.m -1;
    long n = atol(W.Strings[1]);
    long m = atol(W.Strings[2]);
    if (len != n + m + 2) {
        cout << "Bad input file" << endl;
        abort();
    }

    uintT* offsets = newA(uintT,n);
    uintE* edges = newA(uintE,m);

    /*
     * an array of edges contain two elements per edge: neighbor and id of extra info
     * in the array of ei;
     */

    {parallel_for(long i=0; i < n; i++) offsets[i] = atol(W.Strings[i + 3]);}
    {parallel_for(long i=0; i<m; i++) {
            edges[i] = atol(W.Strings[i+n+3]);
        }}
    //W.del(); // to deal with performance bug in malloc

    vertex* v = newA(vertex,n);

    {parallel_for (uintT i=0; i < n; i++) {
            uintT o = offsets[i];
            uintT l = ((i == n-1) ? m : offsets[i+1])-offsets[i];
            v[i].setOutDegree(l);
            v[i].setOutNeighbors(edges+o);
        }}

    uintT* tOffsets = newA(uintT,n);
    {parallel_for(long i=0;i<n;i++) tOffsets[i] = INT_T_MAX;}
    uintE* inEdges = newA(uintE,m);
    intPair* temp = newA(intPair,m);
    {parallel_for(long i=0;i<n;i++){
            uintT o = offsets[i];
            for(uintT j=0;j<v[i].getOutDegree();j++){
                temp[o+j] = make_pair(v[i].getOutNeighbor(j),i);
            }
        }}
    free(offsets);

    quickSort(temp,m,pairFirstCmp<uintE>());

    tOffsets[temp[0].first] = 0;
    inEdges[0] = temp[0].second;
    {parallel_for(long i=1;i<m;i++) {
            inEdges[i] = temp[i].second;
            if(temp[i].first != temp[i-1].first) {
                tOffsets[temp[i].first] = i;
            }
        }}

    free(temp);

    //fill in offsets of degree 0 vertices by taking closest non-zero
    //offset to the right
    sequence::scanIBack(tOffsets,tOffsets,n,minF<uintT>(),(uintT)m);

    {parallel_for(long i=0;i<n;i++){
            uintT o = tOffsets[i];
            uintT l = ((i == n-1) ? m : tOffsets[i+1])-tOffsets[i];
            v[i].setInDegree(l);
            v[i].setInNeighbors(inEdges+o);
        }}

    free(tOffsets);

    //now graph is loaded, structured according to Ligra
    //copying to internal structure
    this->n = n;
    this->m = m;
    long eid = 0;
    for (long i = 0; i < n; i++) {
        this->V.push_back(Vertex(i,0));
        for (long j = 0; j < v[i].outDegree; j++) {
            Edge e(eid++);
            e.fromid = i;
            e.toid = v[i].getOutNeighbor(j);
            e.capacity = 1;
            e.lower = 0;
            e.weight = 1;
            this->E.push_back(e);

        }
    }
}

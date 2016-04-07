//
// Created by alvis on 16.03.16.
//

#ifndef NETWORKFLOW_GRAPH_H
#define NETWORKFLOW_GRAPH_H

#include <vector>
#include <iostream>
#include <sys/time.h>
#include <fstream>
#include <random>
#include <chrono>
#include <lemon/list_graph.h>
#include <lemon/lgf_reader.h>
#include <lemon/lgf_writer.h>
#include <cstring>

//for spatial objects
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "Common.h"
#include "parallel.h"
#include "quickSort.h"
#include "utils.h"

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

//http://www.boost.org/doc/libs/1_60_0/libs/geometry/doc/html/geometry/spatial_indexes/rtree_examples/iterative_query.html
typedef bg::model::point<double, 2, bg::cs::cartesian> point;
typedef std::pair<point,long> value;
typedef bgi::rtree< value, bgi::linear<16> > rtree_t;

class Edge {
public:
    long fromid;
    long toid;
    int weight;
    int capacity;
    int lower;
    long ID;

    Edge(long id) {
        ID = id;
    }
};

class Vertex {
public:
    long ID;
    int supply;
    std::vector<long> E;

    Vertex(long _ID, int _supply) {
        E.clear();
        ID = _ID;
        supply = _supply;
    }
};

class Graph {
public:
    long n;
    long m;
    bool isSpatial;

    //stiring cursor to next edge for each node
    vector<long> QryCnt;
    vector<rtree_t::const_query_iterator> QryIt;

    std::vector<Vertex> V;
    std::vector<Edge> E;
    std::vector<point> coords;

    // additional data for incremental neighbour addition
    vector<vector<long>> fullE; //contains IDs of all outgoing edges
    vector<vector<long>> completeE; //contains IDs of all edges where a node participates (input and output)

    rtree_t rtree;

    void clear_graph() {
        V.clear();
        E.clear();
        m = 0;
        n = 0;
        QryCnt.clear();
        fullE.clear();
        completeE.clear();
    }
    void clear_edge_list() {
        for (long i = 0; i < n; i++) {
            QryCnt[i] = 0;
            V[i].E.clear();
        }
    }

    inline long get_pair(long edgeId, long nodeId) const {
        return (E[edgeId].fromid == nodeId)?E[edgeId].toid:E[edgeId].fromid;
    }
    inline bool is_forward(long edgeId, long nodeId) const {
        return E[edgeId].fromid == nodeId;
    }

    // initialization for simple graph includes storing edges
    // for spatial data - building NN data structure
    void init_neighbors() {
        if (isSpatial) {
            //init r-tree
            for (long i = 0; i < n; i++) {
                rtree.insert(std::make_pair(coords[i],i)); //todo make this r*-tree and bucket-based
                QryIt[i] = rtree.qbegin(bgi::nearest(coords[i], (unsigned int)n)); //todo here is a problem! no more nodes than 65000
                ++QryIt[i];//skip the node itself
            }
        } else {
            //E is already filled
            fullE.resize(n);
            completeE.resize(n); // needed for SCS and Lemon

            parallel_for(long i = 0; i < m; i++) {
                fullE[E[i].fromid].push_back(i); // only forward edges in fullE
                completeE[E[i].fromid].push_back(i);
                completeE[E[i].toid].push_back(i);
            }
            // pointer to the next element for insertion to E_sub
            QryCnt.resize(n,0);
        }
    };
    void myqsort(std::vector<long>& edgelist, long start, long end) //TODO change to quicksort.h
    {
        if (start >= end) return;

        long temp;
        long pivot_ind = start;
        long pivot_val = E[edgelist[start]].weight;
        long i = start+1;
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
        for (long i = 0; i < n; i++)
        {
            myqsort(fullE[i], 0, fullE[i].size()-1);
        }
    }

    //checks if everything was added for particular node
    inline bool isFull(long nodeId) {
        if (isSpatial) {
            return !(QryIt[nodeId] != rtree.qend());
        }
        return QryCnt[nodeId] >= fullE[nodeId].size();
    }

    //returns weight of next candidate for addition
    inline long get_next_neighbour_weight(long nodeId) {
        if (isFull(nodeId)) {
            return -1;
        }
        if (isSpatial) {
            double d = bg::distance(coords[nodeId], QryIt[nodeId]->first);
            return static_cast<long>(d);
        } else {
            return E[fullE[nodeId][QryCnt[nodeId]]].weight;
        }
    };
    //return edge id
    inline long insert_nn_to_edge_list(long nodeId) {
        if (isFull(nodeId)) {
            return -1;
        }
        if (isSpatial) {
            Edge e(E.size());
            e.capacity = 1;
            e.fromid = nodeId;
            e.toid = QryIt[nodeId]->second;
            double d = bg::distance(coords[nodeId], QryIt[nodeId]->first);
            e.weight = static_cast<int>(d);
            e.lower = 0;
            E.push_back(e);
            QryIt[nodeId]++;
            return E.size()-1;
        } else {
            //assuming everything is in E already
            return fullE[nodeId][QryCnt[nodeId]++];
        }
    }

    // graph generators
    void generate_full_bipartite_graph(long size, long param1, long param2 = 1000, long distr = 0);

    void print_graph();
    void save_graph(string filename);
    void load_graph(string filename, string log_filename = "", long experiment_id = 0);
    void load_lgf_graph(string filename);
    void load_adj_graph(string filename);
    void save_graph_info(string filename, long experiment_id);
    void save_graph_blossom(string filname);

//    /*
//     * Spatial data
//     */
    void load_points(string filename, string log_filename = "", long experiment_id = 0);


    // operators
    bool operator==(const Graph &other) const {
        // compare E, V, fullE
        if (V.size() != other.V.size() || E.size() != other.E.size())
            return false;
        for (long i = 0; i < n; i++) {
            if (V[i].E != other.V[i].E ||
                    V[i].ID != other.V[i].ID ||
                    V[i].supply != other.V[i].supply)
                return false;
            if (fullE[i] != other.fullE[i])
                return false;
        }
        for (long i = 0; i < m; i++) {
            if (E[i].ID != other.E[i].ID ||
                    E[i].fromid != other.E[i].fromid ||
                    E[i].toid != other.E[i].toid ||
                    E[i].capacity != other.E[i].capacity ||
                    E[i].lower != other.E[i].lower)
                return false;
        }
        return true;
    }

    /*
     * Part of Ligra to read Adj
     */
    typedef pair<long,long> longPair;
    template <class E>
    struct pairFirstCmp {
        bool operator() (pair<long,E> a, pair<long,E> b) {
            return a.first < b.first; }
    };
    struct vertex {
        long* inNeighbors, *outNeighbors;
        long outDegree;
        long inDegree;
        void del() {free(inNeighbors); free(outNeighbors);}
        vertex(long* iN, long* oN, long id, long od)
        : inNeighbors(iN), outNeighbors(oN), inDegree(id), outDegree(od) {}
        long getInNeighbor(long j) { return inNeighbors[j]; }
        long getOutNeighbor(long j) { return outNeighbors[j]; }
        void setInNeighbors(long* _i) { inNeighbors = _i; }
        void setOutNeighbors(long* _i) { outNeighbors = _i; }
        long getInDegree() { return inDegree; }
        long getOutDegree() { return outDegree; }
        void setInDegree(long _d) { inDegree = _d; }
        void setOutDegree(long _d) { outDegree = _d; }
        void flipEdges() { swap(inNeighbors,outNeighbors); swap(inDegree,outDegree); }
    };
    struct words {
        long n; // total number of characters
        char* Chars;  // array storing all strings
        long m; // number of substrings
        char** Strings; // polongers to strings (all should be null terminated)
        words() {}
        words(char* C, long nn, char** S, long mm)
                : Chars(C), n(nn), Strings(S), m(mm) {}
        void del() {free(Chars); free(Strings);}
    };
    _seq<char> readStringFromFile(char *fileName) {
        ifstream file (fileName, ios::in | ios::binary | ios::ate);
        if (!file.is_open()) {
            std::cout << "Unable to open file: " << fileName << std::endl;
            abort();
        }
        long end = file.tellg();
        file.seekg (0, ios::beg);
        long n = end - file.tellg();
        char* bytes = newA(char,n+1);
        file.read (bytes,n);
        file.close();
        return _seq<char>(bytes,n);
    }
    // parallel code for converting a string to words
    words stringToWords(char *Str, long n) {
        {parallel_for (long i=0; i < n; i++)
                if (isSpace(Str[i])) Str[i] = 0; }

        // mark start of words
        bool *FL = newA(bool,n);
        FL[0] = Str[0];
        {parallel_for (long i=1; i < n; i++) FL[i] = Str[i] && !Str[i-1];}

        // offset for each start of word
        _seq<long> Off = sequence::packIndex<long>(FL, n);
        long m = Off.n;
        long *offsets = Off.A;

        // polonger to each start of word
        char **SA = newA(char*, m);
        {parallel_for (long j=0; j < m; j++) SA[j] = Str+offsets[j];}

        free(offsets); free(FL);
        return words(Str,n,SA,m);
    }

    /*
     * Unit Tests
     */
    bool test_graph_structure();
    bool test_sorting();
    bool test_save_load();
};


#endif //NETWORKFLOW_GRAPH_H

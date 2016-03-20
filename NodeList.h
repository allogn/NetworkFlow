//
// Created by alvis on 18.03.16.
//

#ifndef NETWORKFLOW_NODELIST_H
#define NETWORKFLOW_NODELIST_H

#include <iostream>

#include "Common.h"

using namespace std;

typedef struct NodeExtended {
    uintT nodeId;
    uintT dist;
    struct NodeExtended* next;
} Node;

class NodeList {
public:
    Node** nodes;
    Node* visited;
    uintT size; // number of buckets
    uintT smallest_dist;

    NodeList(uintT _n) {
        size = _n*3+1;
        nodes = new Node*[size];
        visited = NULL;
        for (uintT i = 0; i < size; i++) {
            nodes[i] = NULL;
        }

        smallest_dist = 0;
    }

    void addToBucket(uintT dist, uintT nodeid) {
        // adding only if it is closer than 3*n
        // ignore otherwise
        //@todo add checking if no vertices too far
        if (dist < size) {
            Node* node = new Node;
            node->nodeId = nodeid;
            node->next = nodes[dist];
            nodes[dist] = node;
        }
    }

    void deleteNearest(uintT dist) {
        Node* node = nodes[dist];
        nodes[dist] = nodes[dist]->next;
        delete node;
    }

    uintT getNearest() {
        // smallest_dist is a pointer (id) of a bucket that is non-empty in current state
        // because of monotonicity every smaller bucket must be empty
        while (smallest_dist < size && nodes[smallest_dist] == NULL) {
            smallest_dist++;
        }
        assert(smallest_dist < size); // Error: nothing to pop in shortestPath
        uintT res = nodes[smallest_dist]->nodeId;
        return res;
    }

    void saveAndPopNearest() {
        Node* tmp = nodes[smallest_dist];
        nodes[smallest_dist] = nodes[smallest_dist]->next;
        tmp->next = visited;
        tmp->dist = smallest_dist;
        visited = tmp;
    }

    void popNearest() {
        nodes[smallest_dist] = nodes[smallest_dist]->next;
        //do not delete anything, it exists in visited
    }

    inline uintT getIdVisited() {
        return visited->nodeId;
    }

    inline uintT getDistVisited() {
        return visited->dist;
    }

    inline void popVisited() {
        Node* tmp = visited->next;
        delete visited;
        visited = tmp;
    }

    inline void clearBuckets(bool* _visited) {
        for (uintT i = 0; i < size; i++ ) {
            while (nodes[i] != NULL) {
                uintT tmp = nodes[i]->nodeId;
                if (!_visited[tmp]) {
                    deleteNearest(i);
                    _visited[tmp] = true;
                } else {
                    nodes[i] = nodes[i]->next;
                }
            }
        }
    }

    void printBucket() {
        cout << "Buckets: " << endl;
        for ( int i = 0; i < size; i++) {
            cout << "Bucket " << i << ": ";
            Node* node = nodes[i];
            while (node != NULL) {
                cout << node->nodeId << ",";
                node = node->next;
            }
            cout << endl;
        }
        cout << "Visited: ";
        Node* node = visited;
        while (node != NULL) {
            cout << node->nodeId << "(" << node->dist << "),";
            node = node->next;
        }
        cout << endl;
    }
};

#endif //NETWORKFLOW_NODELIST_H

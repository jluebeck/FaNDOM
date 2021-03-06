//AUTHOR: Siavash Raisi, sraeisid@ucsd.edu

#ifndef FANDOM_SEEDGRAPH_H
#define FANDOM_SEEDGRAPH_H

#include <vector>
#include <list>
#include <stack>
#include <climits>
#include <cmath>
#include <set>

#define INF INT_MAX

using namespace std;

///////GENERIC GRAPH CLASS DEFS
class AdjListNode {
    int v;
    long weight;
public:
    AdjListNode(int _v, long _w) {
        v = _v;
        weight = _w;
    }
    int getV() const { return v; }
    long getWeight() const { return weight; }
};

class Graph {
    int V;    // No. of vertices'
    // Pointer to an array containing adjacency lists
    list<AdjListNode> *adj;
    // A function used by shortestPath
    void topologicalSortUtil(int v, bool visited[], stack<int> &Stack);
public:
    explicit Graph(int V);   // Constructor
    ~Graph();   // destructor
    // function to add an edge to graph
    void addEdge(int u, int v, long weight);
    // Finds shortest paths from given source vertex
    vector<long> shortestPath(int s);
//    vector<long> shortestPath_SV(int s);
};

//////// FUNCTION DEFS

//vector<long>
//solve_graph_straight(vector<pair<int, int>> vertices, vector<double> &q_dist, vector<double> &ref_dist, int w);
//
//vector<long>
//solve_graph_reverse(vector<pair<int, int>> vertices, vector<double> &q_dist, vector<double> &ref_dist, int w);

vector<long> solve_graph(vector<pair<int, int>> vertices, vector<double> &q_dist, vector<double> &ref_dist, int w,
        int dir);


#endif //FANDOM_SEEDGRAPH_H

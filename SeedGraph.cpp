#include "SeedGraph.h"

using namespace std;

//Graph function definitions
Graph::Graph(int V) {
    this->V = V;
    adj = new list<AdjListNode>[V];
}

Graph::~Graph() {
    delete[] adj;
}
void Graph::addEdge(int u, int v, long weight) {
    AdjListNode node(v, weight);
    adj[u].push_back(node); // Add v to u's list
}

void Graph::topologicalSortUtil(int v, bool visited[], stack<int> &Stack) {
    // Mark the current node as visited
    visited[v] = true;
    // Recur for all the vertices adjacent to this vertex
    list<AdjListNode>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i) {
        AdjListNode node = *i;
        if (!visited[node.getV()])
            topologicalSortUtil(node.getV(), visited, Stack);
    }
    // Push current vertex to stack which stores topological sort
    Stack.push(v);
}
vector<long> Graph::shortestPath(int s) {//just return the score and path
    stack<int> Stack;
    vector<long> dist = vector<long>(V);
    vector<vector<int>> path = vector<vector<int>>(V);
    // Mark all the vertices as not visited
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++)
        visited[i] = false;
    // Call the recursive helper function to store Topological Sort
    // starting from all vertices one by one
    topologicalSortUtil(s, visited, Stack);
    // Initialize distances to all vertices as infinite and distance
    // to source as 0
    for (int i = 0; i < V; i++)
        dist[i] = INF;
    dist[s] = 0;
    // Process vertices in topological order
    while (! Stack.empty() ) {
        // Get the next vertex from topological order
        int u = Stack.top();
        Stack.pop();
        // Update distances of all adjacent vertices
        list<AdjListNode>::iterator i;
        if (dist[u] != INF) {
            for (i = adj[u].begin(); i != adj[u].end(); ++i)
                if (dist[i->getV()] > dist[u] + i->getWeight()) {
                    dist[i->getV()] = dist[u] + i->getWeight();
                    path[i->getV()] = path[u];
                    path[i->getV()].push_back(u);
                }
        }
    }
    delete[] visited;
    return dist;
}

/////////////calculate shortest patch in straight DAG////////////
//vector<long>
//solve_graph_straight(vector<pair<int, int>> vertices, vector<double> &q_dist, vector<double> &ref_dist, const int w) {
//    Graph *g;
//    int size_vertices = vertices.size();
//    g = new Graph(2 + w * size_vertices);
//    int source = 0;
//    int end = 2 + w * size_vertices - 1;
//    for (int ard = 0; ard < size_vertices; ard++) {
//        g[0].addEdge((ard + 1) * w - 1, (ard + 1) * w, 0);
//        g[0].addEdge((ard + 1) * w - 2, (ard + 1) * w - 1, 0);
//    }
//    for (int ard = 0; ard < size_vertices; ard++) {
//        for (int srd = 0; srd < size_vertices; srd++) {
//            if (srd != ard) {
//                int first_end = vertices[ard].first + w - 1;
//                int second_start = vertices[srd].first;
//                int first_end_ref = vertices[ard].second + w - 1;
//                int second_start_ref = vertices[srd].second;
//                if (first_end <= second_start && first_end_ref <= second_start_ref) {
//                    long X_ref = long(ref_dist[second_start_ref - 1] - ref_dist[first_end_ref - 1]);
//                    long Y_query = long(q_dist[second_start - 1] - q_dist[first_end - 1]);
//                    long weight_edge = long(sqrt(X_ref * X_ref + Y_query * Y_query));
//                    g[0].addEdge((ard + 1) * w, (srd + 1) * w - 2, weight_edge);
//                }
//                if (first_end - 1 == second_start && first_end_ref - 1 <= second_start_ref) {
//                    long X_ref = long(ref_dist[second_start_ref - 1] - ref_dist[first_end_ref - 2]);
//                    g[0].addEdge((ard + 1) * w - 1, (srd + 1) * w - 2,long(sqrt(X_ref * X_ref)));
//                }
//
//            }
//        }
//        long query_to_start_dist = q_dist[vertices[ard].first - 1] - q_dist[0];
//        long ref_length = 0;
//        g[0].addEdge(source, (ard + 1) * w - 2,
//                     long(sqrt(2 * query_to_start_dist * query_to_start_dist + ref_length * ref_length)));
//        query_to_start_dist = q_dist[q_dist.size() - 1] - q_dist[vertices[ard].first + w - 1 - 1];
//        ref_length = 0;
//        g[0].addEdge((ard + 1) * w, end,
//                     long(sqrt(2 * query_to_start_dist * query_to_start_dist + ref_length * ref_length)));
//    }
//    vector<long> score = g[0].shortestPath(source);
//    delete g;
//    return score;
//}
//
/////////////calculate shortest patch in reverse DAG////////////
//vector<long>
//solve_graph_reverse(vector<pair<int, int>> vertices, vector<double> &q_dist, vector<double> &ref_dist, const int w) {
//    Graph *g;
//    int size_vertices = vertices.size();
//    g = new Graph(2 + w * size_vertices);
//    int source = 0;
//    int end = 2 + w * size_vertices - 1;
//    for (int ard = 0; ard < size_vertices; ard++) {
//        g[0].addEdge((ard + 1) * w - 1, (ard + 1) * w, 0);
//        g[0].addEdge((ard + 1) * w - 2, (ard + 1) * w - 1, 0);
//    }
//    for (int ard = 0; ard < size_vertices; ard++) {
//        for (int srd = 0; srd < size_vertices; srd++) {
//            if (srd != ard) {
//                int first_end = vertices[ard].first - w + 1;
//                int second_start = vertices[srd].first;
//                int first_end_ref = vertices[ard].second + w - 1;
//                int second_start_ref = vertices[srd].second;
//                if (first_end >= second_start && first_end_ref <= second_start_ref) {
//                    long X_ref = long(ref_dist[second_start_ref - 1] - ref_dist[first_end_ref - 1]);
//                    long Y_query = long(q_dist[second_start - 1] - q_dist[first_end - 1]);
//                    long weight_edge = long(sqrt(X_ref * X_ref + Y_query * Y_query));
//                    g[0].addEdge((ard + 1) * w, (srd + 1) * w - 2, weight_edge);
//                }
//                if (first_end + 1 == second_start && first_end_ref - 1 <= second_start_ref) {
//                    long X_ref = long(ref_dist[second_start_ref - 1] - ref_dist[first_end_ref - 2]);
//                    long weight_edge = long(sqrt(X_ref * X_ref));
//                    g[0].addEdge((ard + 1) * w - 1, (srd + 1) * w - 2, weight_edge);
//                }
//            }
//        }
//        long query_to_start_dist = q_dist[vertices[ard].first - 1] - q_dist[q_dist.size() - 1];
//        long ref_length = 0;
//        g[0].addEdge(source, (ard + 1) * w - 2,long(sqrt(2 * query_to_start_dist * query_to_start_dist + ref_length * ref_length)));
//        query_to_start_dist = q_dist[vertices[ard].first - 1 - w + 1] - q_dist[0];
//        ref_length = 0;
//        g[0].addEdge((ard + 1) * w, end,
//                     long(sqrt(2 * query_to_start_dist * query_to_start_dist + ref_length * ref_length)));
//    }
//    vector<long> score = g[0].shortestPath(source);
//    delete g;
//    return score;
//}


/////////////calculate shortest patch in straight DAG////////////
vector<long> solve_graph(vector<pair<int, int>> vertices, vector<double> &q_dist, vector<double> &ref_dist,
        const int w, const int dir) {

    Graph *g;
    int size_vertices = vertices.size();
    g = new Graph(2 + w * size_vertices);
    int source = 0;
    int end = 2 + w * size_vertices - 1;
    for (int ard = 0; ard < size_vertices; ard++) {
        g[0].addEdge((ard + 1) * w - 1, (ard + 1) * w, 0);
        g[0].addEdge((ard + 1) * w - 2, (ard + 1) * w - 1, 0);
    }
    for (int ard = 0; ard < size_vertices; ard++) {
        for (int srd = 0; srd < size_vertices; srd++) {
            if (srd != ard) {
                int first_end = vertices[ard].first + dir*(w - 1);
                int second_start = vertices[srd].first;
                int first_end_ref = vertices[ard].second + w - 1;
                int second_start_ref = vertices[srd].second;
                if (dir*first_end <= dir*second_start && first_end_ref <= second_start_ref) {
                    long X_ref = long(ref_dist[second_start_ref - 1] - ref_dist[first_end_ref - 1]);
                    long Y_query = long(q_dist[second_start - 1] - q_dist[first_end - 1]);
                    long weight_edge = long(sqrt(X_ref * X_ref + Y_query * Y_query));
                    g[0].addEdge((ard + 1) * w, (srd + 1) * w - 2, weight_edge);
                }
                if (first_end - dir == second_start && first_end_ref - 1 <= second_start_ref) {
                    long X_ref = long(ref_dist[second_start_ref - 1] - ref_dist[first_end_ref - 2]);
                    g[0].addEdge((ard + 1) * w - 1, (srd + 1) * w - 2, long(abs(X_ref)));
                }
            }
        }
        long query_to_start_dist = q_dist[vertices[ard].first - 1];
        long query_to_end_dist;
        if (dir > 0) {
            query_to_start_dist -= q_dist[0];
            query_to_end_dist = q_dist[q_dist.size() - 1] - q_dist[vertices[ard].first + w - 2]; //(+w - 1 -1 )

        } else {
            query_to_start_dist -= q_dist[q_dist.size() - 1];
            query_to_end_dist = q_dist[vertices[ard].first - w] - q_dist[0]; //(-1 - w + 1)
        }

        g[0].addEdge(source, (ard + 1) * w - 2,long(sqrt(2 * pow(query_to_start_dist,2))));
        g[0].addEdge((ard + 1) * w, end,long(sqrt(2 * pow(query_to_end_dist, 2))));
    }
    vector<long> score = g[0].shortestPath(source);
    delete g;
    return score;
}
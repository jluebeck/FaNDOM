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


pair<long, vector<int>> Graph::shortestPath(int s) {//just return the score and path
    stack<int> Stack;
    long dist[V];
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
    return make_pair(dist[V - 1], path[V - 1]);
}

vector<long> Graph::shortestPath_SV(int s) {//return the score to all other vertices
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

///////////calculate shortest patch in straight DAG////////////
pair<long, vector<pair<int, int>>>
solve_graph_straight(vector<pair<int, int>> vertices, vector<double> &q_dist, vector<double> &ref_dist, const int w) {
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
                int payan_avali = vertices[ard].first + w - 1;
                int avale_dovomi = vertices[srd].first;
                int payan_avali_ref = vertices[ard].second + w - 1;
                int avale_dovomi_ref = vertices[srd].second;
                if (payan_avali <= avale_dovomi && payan_avali_ref <= avale_dovomi_ref) {
                    long dif_ref =
                            long(
                                    ref_dist[avale_dovomi_ref - 1] - ref_dist[payan_avali_ref - 1]);
                    long query_dif = long(
                            q_dist[avale_dovomi - 1] - q_dist[payan_avali - 1]);
                    long weight_edge = long(sqrt(dif_ref * dif_ref + query_dif * query_dif));
                    g[0].addEdge((ard + 1) * w, (srd + 1) * w - 2, weight_edge);
                }
                if (payan_avali - 1 == avale_dovomi && payan_avali_ref - 1 <= avale_dovomi_ref) {
                    long dif_ref = long(
                            ref_dist[avale_dovomi_ref - 1] - ref_dist[payan_avali_ref - 2]);
                    g[0].addEdge((ard + 1) * w - 1, (srd + 1) * w - 2,
                                 long(sqrt(dif_ref * dif_ref)));
                }

            }
        }
//////////////////////////////////
///////////query_length - first label?/////////////
////////////////////////////////////////////
        long query_length = q_dist[vertices[ard].first - 1] - q_dist[0];
        long ref_length = 0;
        g[0].addEdge(source, (ard + 1) * w - 2,
                     long(sqrt(2 * query_length * query_length + ref_length * ref_length)));
        query_length =
                q_dist[q_dist.size() - 1] - q_dist[vertices[ard].first + w - 1 - 1];
        ref_length = 0;
        g[0].addEdge((ard + 1) * w, end,
                     long(sqrt(2 * query_length * query_length + ref_length * ref_length)));
    }
    pair<long, vector<int>> score = g[0].shortestPath(source);
    set<int> score_point;
    for (int i = 1; i < score.second.size(); i++) {
        score_point.insert(ceil(float(score.second[i]) / 3));
    }
    vector<pair<int, int>> best_path;
    for (auto elem :score_point) {
        best_path.push_back(vertices[elem - 1]);
    }
    delete g;
    return make_pair(score.first, best_path);
}


///////////calculate shortest patch in reverse DAG////////////
pair<long, vector<pair<int, int>>>
solve_graph_reverse(vector<pair<int, int>> vertices, vector<double> &q_dist, vector<double> &ref_dist, const int w) {
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
                int payan_avali = vertices[ard].first - w + 1;
                int avale_dovomi = vertices[srd].first;
                int payan_avali_ref = vertices[ard].second + w - 1;
                int avale_dovomi_ref = vertices[srd].second;
                if (payan_avali >= avale_dovomi && payan_avali_ref <= avale_dovomi_ref) {
                    long dif_ref =
                            long(
                                    ref_dist[avale_dovomi_ref - 1] - ref_dist[payan_avali_ref - 1]);
                    long query_dif = long(
                            q_dist[avale_dovomi - 1] - q_dist[payan_avali - 1]);
                    long weight_edge = long(sqrt(dif_ref * dif_ref + query_dif * query_dif));
                    g[0].addEdge((ard + 1) * w, (srd + 1) * w - 2, weight_edge);
                }
                if (payan_avali + 1 == avale_dovomi && payan_avali_ref - 1 <= avale_dovomi_ref) {
                    long dif_ref = long(ref_dist[avale_dovomi_ref - 1] - ref_dist[payan_avali_ref - 2]);
                    long weight_edge = long(sqrt(dif_ref * dif_ref));
                    g[0].addEdge((ard + 1) * w - 1, (srd + 1) * w - 2, weight_edge);
                }
            }
        }
        long query_length = q_dist[vertices[ard].first - 1] - q_dist[q_dist.size() - 1];
        long ref_length = 0;
        g[0].addEdge(source, (ard + 1) * w - 2,
                     long(sqrt(2 * query_length * query_length + ref_length * ref_length)));
        query_length = q_dist[vertices[ard].first - 1 - w + 1] - q_dist[0];
        ref_length = 0;
        g[0].addEdge((ard + 1) * w, end,
                     long(sqrt(2 * query_length * query_length + ref_length * ref_length)));
    }
    pair<long, vector<int>> score = g[0].shortestPath(source);
    set<int> score_point;
    for (int i = 1; i < score.second.size(); i++) {
        score_point.insert(ceil(float(score.second[i]) / 3));
    }
    vector<pair<int, int>> best_path;
    for (auto elem :score_point) {
        best_path.push_back(vertices[elem - 1]);
    }

    delete g;
    return make_pair(score.first, best_path);
}

//////////////calculate shortest patch in straight DAG////////////////
vector<long>
solve_graph_straight_SV(vector<pair<int, int>> vertices, vector<double> &q_dist, vector<double> &ref_dist, int dir,
        const int w) {
    Graph *g;
    int size_vertices = vertices.size();
    g = new Graph(2 + w * size_vertices);
    int source = 0;
    int end = 2 + w * size_vertices - 1;
    if (dir == 1) {
        for (int ard = 0; ard < size_vertices; ard++) {
            g[0].addEdge((ard + 1) * w - 1, (ard + 1) * w, 0);
            g[0].addEdge((ard + 1) * w - 2, (ard + 1) * w - 1, 0);
        }
        for (int ard = 0; ard < size_vertices; ard++) {
            for (int srd = 0; srd < size_vertices; srd++) {
                if (srd != ard) {
                    int payan_avali = vertices[ard].first + w - 1;
                    int avale_dovomi = vertices[srd].first;
                    int payan_avali_ref = vertices[ard].second + w - 1;
                    int avale_dovomi_ref = vertices[srd].second;
                    if (payan_avali <= avale_dovomi && payan_avali_ref <= avale_dovomi_ref) {
                        long dif_ref =
                                long(
                                        ref_dist[avale_dovomi_ref - 1] - ref_dist[payan_avali_ref - 1]);
                        long query_dif = long(
                                q_dist[avale_dovomi - 1] - q_dist[payan_avali - 1]);
                        long weight_edge = long(sqrt(dif_ref * dif_ref + query_dif * query_dif));
                        g[0].addEdge((ard + 1) * w, (srd + 1) * w - 2, weight_edge);
                    }
                    if (payan_avali - 1 == avale_dovomi && payan_avali_ref - 1 <= avale_dovomi_ref) {
                        long dif_ref = long(
                                ref_dist[avale_dovomi_ref - 1] - ref_dist[payan_avali_ref - 2]);
                        g[0].addEdge((ard + 1) * w - 1, (srd + 1) * w - 2,
                                     long(sqrt(dif_ref * dif_ref)));
                    }
                }
            }
            long query_length = q_dist[vertices[ard].first - 1] - q_dist[0];
            long ref_length = 0;
            g[0].addEdge(source, (ard + 1) * w - 2,
                         long(sqrt(2 * query_length * query_length + ref_length * ref_length)));
            query_length =
                    q_dist[q_dist.size() - 1] - q_dist[vertices[ard].first + w - 1 - 1];
            ref_length = 0;
            g[0].addEdge((ard + 1) * w, end,
                         long(sqrt(2 * query_length * query_length + ref_length * ref_length)));
        }
        vector<long> score = g[0].shortestPath_SV(source);
        delete g;
        return score;
    } else {
        for (int ard = 0; ard < size_vertices; ard++) {
            g[0].addEdge((ard + 1) * w, (ard + 1) * w - 1, 0);
            g[0].addEdge((ard + 1) * w - 1, (ard + 1) * w - 2, 0);
        }
        for (int ard = 0; ard < size_vertices; ard++) {
            for (int srd = 0; srd < size_vertices; srd++) {
                if (srd != ard) {
                    int payan_avali = vertices[ard].first + w - 1;
                    int avale_dovomi = vertices[srd].first;
                    int payan_avali_ref = vertices[ard].second + w - 1;
                    int avale_dovomi_ref = vertices[srd].second;
                    if (payan_avali <= avale_dovomi && payan_avali_ref <= avale_dovomi_ref) {
                        long dif_ref =
                                long(
                                        ref_dist[avale_dovomi_ref - 1] - ref_dist[payan_avali_ref - 1]);
                        long query_dif = long(
                                q_dist[avale_dovomi - 1] - q_dist[payan_avali - 1]);
                        long weight_edge = long(sqrt(dif_ref * dif_ref + query_dif * query_dif));
                        g[0].addEdge((srd + 1) * w - 2, (ard + 1) * w, weight_edge);
                    }
                    if (payan_avali - 1 == avale_dovomi && payan_avali_ref - 1 <= avale_dovomi_ref) {
                        long dif_ref = long(
                                ref_dist[avale_dovomi_ref - 1] - ref_dist[payan_avali_ref - 2]);
                        g[0].addEdge((srd + 1) * w - 2, (ard + 1) * w - 1,
                                     long(sqrt(dif_ref * dif_ref)));
                    }
                }
            }
            long query_length = q_dist[vertices[ard].first - 1] - q_dist[0];
            long ref_length = 0;
            g[0].addEdge((ard + 1) * w - 2, source,
                         long(sqrt(2 * query_length * query_length + ref_length * ref_length)));
            query_length =
                    q_dist[q_dist.size() - 1] - q_dist[vertices[ard].first + w - 1 - 1];
            ref_length = 0;
            g[0].addEdge(end, (ard + 1) * w,
                         long(sqrt(2 * query_length * query_length + ref_length * ref_length)));
        }
        vector<long> score = g[0].shortestPath_SV(end);
        delete g;
        return score;
    }
}


/////////////////calculate shortest patch in reverse DAG/////////////////////////////
vector<long>
solve_graph_reverse_SV(vector<pair<int, int>> vertices, vector<double> &q_dist, vector<double> &ref_dist, int dir,
        const int w) {
    Graph *g;
    int size_vertices = vertices.size();
    g = new Graph(2 + w * size_vertices);
    int source = 0;
    int end = 2 + w * size_vertices - 1;
    if (dir == 1) {
        for (int ard = 0; ard < size_vertices; ard++) {
            g[0].addEdge((ard + 1) * w - 1, (ard + 1) * w, 0);
            g[0].addEdge((ard + 1) * w - 2, (ard + 1) * w - 1, 0);
        }
        for (int ard = 0; ard < size_vertices; ard++) {
            for (int srd = 0; srd < size_vertices; srd++) {
                if (srd != ard) {
                    int payan_avali = vertices[ard].first - w + 1;
                    int avale_dovomi = vertices[srd].first;
                    int payan_avali_ref = vertices[ard].second + w - 1;
                    int avale_dovomi_ref = vertices[srd].second;
                    if (payan_avali >= avale_dovomi && payan_avali_ref <= avale_dovomi_ref) {
                        long dif_ref =
                                long(
                                        ref_dist[avale_dovomi_ref - 1] - ref_dist[payan_avali_ref - 1]);
                        long query_dif = long(
                                q_dist[avale_dovomi - 1] - q_dist[payan_avali - 1]);
                        long weight_edge = long(sqrt(dif_ref * dif_ref + query_dif * query_dif));
                        g[0].addEdge((ard + 1) * w, (srd + 1) * w - 2, weight_edge);
                    }
                    if (payan_avali + 1 == avale_dovomi && payan_avali_ref - 1 <= avale_dovomi_ref) {
                        long dif_ref = long(ref_dist[avale_dovomi_ref - 1] - ref_dist[payan_avali_ref - 2]);
                        long weight_edge = long(sqrt(dif_ref * dif_ref));
                        g[0].addEdge((ard + 1) * w - 1, (srd + 1) * w - 2, weight_edge);
                    }
                }
            }
            long query_length = q_dist[vertices[ard].first - 1] - q_dist[q_dist.size() - 1];
            long ref_length = 0;
            g[0].addEdge(source, (ard + 1) * w - 2,
                         long(sqrt(2 * query_length * query_length + ref_length * ref_length)));
            query_length = q_dist[vertices[ard].first - 1 - w + 1] - q_dist[0];
            ref_length = 0;
            g[0].addEdge((ard + 1) * w, end,
                         long(sqrt(2 * query_length * query_length + ref_length * ref_length)));
        }
        vector<long> score = g[0].shortestPath_SV(source);
        delete g;
        return score;
    } else {
        for (int ard = 0; ard < size_vertices; ard++) {
            g[0].addEdge((ard + 1) * w, (ard + 1) * w - 1, 0);
            g[0].addEdge((ard + 1) * w - 1, (ard + 1) * w - 2, 0);
        }
        for (int ard = 0; ard < size_vertices; ard++) {
            for (int srd = 0; srd < size_vertices; srd++) {
                if (srd != ard) {
                    int payan_avali = vertices[ard].first - w + 1;
                    int avale_dovomi = vertices[srd].first;
                    int payan_avali_ref = vertices[ard].second + w - 1;
                    int avale_dovomi_ref = vertices[srd].second;
                    if (payan_avali >= avale_dovomi && payan_avali_ref <= avale_dovomi_ref) {
                        long dif_ref =
                                long(
                                        ref_dist[avale_dovomi_ref - 1] - ref_dist[payan_avali_ref - 1]);
                        long query_dif = long(
                                q_dist[avale_dovomi - 1] - q_dist[payan_avali - 1]);
                        long weight_edge = long(sqrt(dif_ref * dif_ref + query_dif * query_dif));
                        g[0].addEdge((srd + 1) * w - 2, (ard + 1) * w, weight_edge);
                    }
                    if (payan_avali + 1 == avale_dovomi && payan_avali_ref - 1 <= avale_dovomi_ref) {
                        long dif_ref = long(ref_dist[avale_dovomi_ref - 1] - ref_dist[payan_avali_ref - 2]);
                        long weight_edge = long(sqrt(dif_ref * dif_ref));
                        g[0].addEdge((srd + 1) * w - 2, (ard + 1) * w - 1, weight_edge);
                    }
                }
            }
            long query_length = q_dist[vertices[ard].first - 1] - q_dist[q_dist.size() - 1];
            long ref_length = 0;
            g[0].addEdge((ard + 1) * w - 2, source,
                         long(sqrt(2 * query_length * query_length + ref_length * ref_length)));
            query_length = q_dist[vertices[ard].first - 1 - w + 1] - q_dist[0];
            ref_length = 0;
            g[0].addEdge(end, (ard + 1) * w,
                         long(sqrt(2 * query_length * query_length + ref_length * ref_length)));
        }
        vector<long> score = g[0].shortestPath_SV(end);
        delete g;
        return score;
    }
}



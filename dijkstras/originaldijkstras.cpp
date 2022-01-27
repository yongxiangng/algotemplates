#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
using namespace __gnu_pbds;
#define MIN INT32_MIN;
#define MAX INT32_MAX
#define LSOne(S) ((S) & (-(S)))
typedef long long ll;
typedef long double ld;
typedef pair<int, int> ii;
typedef tuple<int, int, int> iii;
typedef vector<char> vc;
typedef vector<int> vi;
typedef vector<ii> vii;
typedef vector<vector<char>> vvc;
typedef vector<vector<int>> vvi;
typedef unordered_map<int, vi> Al;
typedef unordered_map<int, vi> AL;
typedef tree<int,null_type,less<>,rb_tree_tag,tree_order_statistics_node_update> ordered_set;
typedef tree<ii,null_type,less<>,rb_tree_tag,tree_order_statistics_node_update> ordered_set_ii;
typedef tree<int,null_type,less_equal<>,rb_tree_tag,tree_order_statistics_node_update> ordered_multiset;
const double EPS = 1e-9;

/**
 * https://github.com/stevenhalim/cpbook-code/blob/master/ch4/sssp/dijkstra.cpp
 *
 * This potentially works faster than the modified dijkstras...
 *
 * E.g. https://open.kattis.com/problems/emptyingbaltic -> modified dijkstras fails
 * (or at least from my attempt, but original dijkstras passes)
 *
 * Modified dijkstras deletes edges "lazily" which might be good since it then solves
 * the issues caused by the greedy choosing which assumes no negative edges (but still
 * doesn't work when there is negative weight cycle)
 */
int main() {
    int n;                                                                // #of vertices
    int s;                                                                // starting vertex
    unordered_map<int, vii> al;
    vi dist(n, MAX); dist[s] = 0;                                   // INF = MAX

    set<ii> pq;                                                           // balanced BST version
    for (int u = 0; u < n; ++u)                                           // dist[u] = INF
      pq.insert({dist[u], u});                                   // but dist[s] = 0
    // sort the pairs by non-decreasing distance from s
    while (!pq.empty()) {                                                  // main loop
        auto [d, u] = *pq.begin();                                // shortest unvisited u
        pq.erase(pq.begin());
        for (auto &[v, w] : al[u]) {                              // all edges from u
            if (dist[u]+w >= dist[v]) continue;                             // not improving, skip
            pq.erase(pq.find({dist[v], v}));               // erase old pair
            dist[v] = dist[u]+w;                                            // relax operation
            pq.insert({dist[v], v});                               // enqueue better pair
        }
    }
}

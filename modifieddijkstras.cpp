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
 * https://open.kattis.com/problems/shortestpath1
 *
 * Modified dijkstras over here. Pros of modified dijkstras -> correctness when negative weight
 * (not negative cycle) disadvantage -> CP4 pg 233 (edge case), https://visualgo.net/en/sssp?slide=8-5
 */
int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);

    int n, m, q, s;
    while (cin >> n >> m >> q >> s && !(n == 0 && m == 0 && q == 0 && s == 0)) {
        unordered_map<int, vii> al;

        for (int i=0;i<m;i++) {
            int u, v, w;
            cin >> u >> v >> w;
            al[u].emplace_back(v, w);
        }
        vi dist(n, MAX); dist[s] = 0;
        priority_queue<ii, vector<ii>, greater<ii>> pq;
        pq.emplace(0, s); // starting vertex

        while (!pq.empty()) {
            auto [d, u] = pq.top(); pq.pop();
            if (d > dist[u]) continue;
            for (auto &[v, w] : al[u]) {            // all edges from u
                if (dist[u] + w >= dist[v]) continue;         // not improving
                dist[v] = dist[u] + w;
                pq.emplace(dist[v], v);
            }
        }
        for (int i=0;i<q;i++) {
            int temp;
            cin >> temp;
            if (dist[temp] == MAX) {
                cout << "Impossible" << endl;
            } else {
                cout << dist[temp] << endl;
            }
        }
        cout << endl;
    }
}

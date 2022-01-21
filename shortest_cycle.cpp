#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
using namespace __gnu_pbds;
#define MIN INT32_MIN
#define MAX INT32_MAX
#define LSOne(S) ((S) & (-(S)))
using ll = long long; using ld = long double;
using ii = pair<int, int>; using iii = tuple<int, int, int>;
using vb = vector<bool>; using vc = vector<char>; using vi = vector<int>; using vl = vector<long>;
using vii = vector<ii>; using vll = vector<ll>;
using vvb = vector<vector<bool>>; using vvc = vector<vector<char>>; using vvi = vector<vector<int>>; using vvl= vector<vector<long>>;
using Al = unordered_map<int, vi>; using AL = unordered_map<int, vi>;
using ordered_set = tree<int,null_type,less<>,rb_tree_tag,tree_order_statistics_node_update>;
using ordered_set_ii = tree<ii,null_type,less<>,rb_tree_tag,tree_order_statistics_node_update>;
using ordered_multiset =  tree<int,null_type,less_equal<>,rb_tree_tag,tree_order_statistics_node_update>; // doesn't work I think
const double EPS = 1e-7;

// https://www.geeksforgeeks.org/shortest-cycle-in-an-undirected-unweighted-graph/
// https://open.kattis.com/problems/beehives2

Al gr;

int shortest_cycle(int n)
{
    // To store length of the shortest cycle
    int ans = INT_MAX;
 
    // For all vertices
    for (int i = 0; i < n; i++) {
 
        // Make distance maximum
        vector<int> dist(n, (int)(1e9));
 
        // Take a imaginary parent
        vector<int> par(n, -1);
 
        // Distance of source to source is 0
        dist[i] = 0;
        queue<int> q;
 
        // Push the source element
        q.push(i);
 
        // Continue until queue is not empty
        while (!q.empty()) {
 
            // Take the first element
            int x = q.front();
            q.pop();
 
            // Traverse for all it's childs
            for (int child : gr[x]) {
 
                // If it is not visited yet
                if (dist[child] == (int)(1e9)) {
 
                    // Increase distance by 1
                    dist[child] = 1 + dist[x];
 
                    // Change parent
                    par[child] = x;
 
                    // Push into the queue
                    q.push(child);
                }
 
                // If it is already visited
                else if (par[x] != child and par[child] != x)
                    ans = min(ans, dist[x] + dist[child] + 1);
            }
        }
    }
 
    // If graph contains no cycle
    if (ans == INT_MAX)
        return -1;
 
    // If graph contains cycle
    else
        return ans;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);
    
    int n,m; cin >> n >> m;
    
    for (int i=0;i<m;i++) {
		int u,v;
		cin >> u >> v;
		gr[u].push_back(v);
		gr[v].push_back(u);
	}
	int ans = shortest_cycle(n);
	if (ans == -1) {
		cout << "impossible" << endl;
	} else {
		cout << ans << endl;
	}
	
}

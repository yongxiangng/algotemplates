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
typedef vector<int> vi;
typedef vector<ii> vii;
typedef tree<int,null_type,less<>,rb_tree_tag,tree_order_statistics_node_update> ordered_set;
typedef tree<ii,null_type,less<>,rb_tree_tag,tree_order_statistics_node_update> ordered_set_ii;
typedef tree<int,null_type,less_equal<>,rb_tree_tag,tree_order_statistics_node_update> ordered_multiset;
const double EPS = 1e-9;

//struct pair_hash {
//    inline size_t operator()(const pair<int,int> & v) const {
//        return v.first*31+v.second;
//    }
//};

//vector<string> split(string s, string del = " ") {
//    vector<string> vec;
//    int start = 0;
//    int end = s.find(del);
//    while (end != -1) {
//        vec.push_back(s.substr(start, end - start));
//        start = end + del.size();
//        end = s.find(del, start);
//    }
//    vec.push_back(s.substr(start, end - start));
//    return vec;
//}

// what is the minimum cost if we are at vertex u and have visited vertices
// that are described by the off (0 bit in mask?) Here 1 bit means unvisted
// and the city is available.
int TSP(int u, int mask, vector<vector<int>> &dist, vector<vector<int>> &memo) {
    if (mask == 0) {                                                                    // mask = free coordinates
        return dist[u][0];                                                              // close the tour
    }
    int &ans = memo[u][mask];
    if (ans != -1) return ans;                                                          // computed before
    ans = 2000000000;
    int m = mask;
    while (m) {                                                                         // up to O(n)
        int two_pow_v = LSOne(m);                                                       // but this is fast
        int v = __builtin_ctz(two_pow_v) + 1;                                           // offset v by + 1 !!! note the possible bug here, confusing because n * 2 ^ (n - 1) runtime
        ans = min(ans, dist[u][v] + TSP(v, mask^two_pow_v, dist, memo));       // keep the min
        m -= two_pow_v;
    }
    return ans;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);
    cout.tie(NULL);

    int tc;
    cin >> tc;
    while (tc --) {
        int x, y;
        cin >> x >> y;
        ii start;
        cin >> start.first >> start.second;
        int n;
        cin >> n;
        vii arr(n);
        for (int i=0;i<n;i++) {
            cin >> arr[i].first >> arr[i].second;
        }
        vector<vector<int>> memo(n + 1, vector<int>((1 << (n))));
        for (auto &a : memo) {
            for (auto &b : a) {
                b = -1;
            }
        }

        vector<vector<int>> dist(n + 1, vector<int>(n + 1));

        for (int i=0;i<n+1;i++) {
            for (int j=0;j<n+1;j++) {
                ii a = i == 0 ? start : arr[i - 1];
                ii b = j == 0 ? start : arr[j - 1];
                dist[i][j] = abs(a.first - b.first) + abs(a.second - b.second);
            }
        }

        cout << TSP(0, (1 << (n)) - 1, dist, memo) << endl;
    }
}

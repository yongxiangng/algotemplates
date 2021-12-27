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

//int directions[4][2] = {{0,1},{0,-1},{1,0},{-1,0}};

//vi S;
//vi T;

//void Kosaraju(int u, int pass, int &size, vi &dfs_num, AL &al, AL &alt) {
//    size++;
//    dfs_num[u] = 1;
//    vi &neighbour = (pass == 1) ? al[u] : alt[u];
//    for (auto &v : neighbour) {
//        if (dfs_num[v] == -1) {
//            Kosaraju(v, pass, size, dfs_num, al, alt);
//        }
//    }
//    S.push_back(u);
//    T.push_back(u);
//}

void process(int u,
             unordered_map<int, vector<pair<int, double>>> &al, // u to [v, w]
             priority_queue<pair<double, int>> &pq,
             vi &taken) {
    taken[u] = 1;
    for (auto &[v, w] : al[u]) {
        if (!taken[v]) {
            pq.emplace(-w, -v);
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);

    cout << fixed << setprecision(3);

    int tc;
    cin >> tc;
    while (tc--) {
        unordered_map<int, vector<pair<int, double>>> al;
        int n;
        cin >> n;
        vector<pair<double, double>> vertices;
        for (int i=0;i<n;i++) {
            double x, y;
            cin >> x >> y;
            vertices.emplace_back(x, y);
            for (int j=0;j<i;j++) {
                double vx = vertices[j].first;
                double vy = vertices[j].second;
                double dist = sqrt((vx - x) * (vx - x) + (vy - y) * (vy - y));
                al[i].push_back({j, dist});
                al[j].push_back({i, dist});
            }
        }
        vi taken;
        priority_queue<pair<double, int>> pq;

        taken.assign(n, 0);
        process(0, al, pq, taken);
        double mst_cost = 0.0;
        int numTaken = 0;

        while (!pq.empty()) {
            auto [w, u] = pq.top(); pq.pop();
            w = -w; u = -u;
            if (taken[u]) continue;
            mst_cost += w;
            process(u, al, pq, taken);
            numTaken++;
            if (numTaken == n - 1) break;
        }

        cout << mst_cost << endl;
    }
}

#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
using namespace __gnu_pbds;
#define MIN INT32_MIN
#define MAX INT32_MAX
#define LSOne(S) ((S) & (-(S)))
typedef long long ll; typedef long double ld;
typedef pair<int, int> ii; typedef tuple<int, int, int> iii;
typedef vector<bool> vb; typedef vector<char> vc; typedef vector<int> vi; typedef vector<long> vl;
typedef vector<ii> vii; typedef vector<ll> vll;
typedef vector<vector<bool>> vvb; typedef vector<vector<char>> vvc; typedef vector<vector<int>> vvi; typedef vector<vector<long>> vvl;
typedef unordered_map<int, vi> Al; typedef unordered_map<int, vi> AL;
typedef tree<int,null_type,less<>,rb_tree_tag,tree_order_statistics_node_update> ordered_set;
typedef tree<ii,null_type,less<>,rb_tree_tag,tree_order_statistics_node_update> ordered_set_ii;
typedef tree<int,null_type,less_equal<>,rb_tree_tag,tree_order_statistics_node_update> ordered_multiset; // doesn't work I think
const double EPS = 1e-7;

struct pair_hash {
    inline size_t operator()(const pair<int,int> & v) const {
        return v.first*31+v.second;
    }
};

//struct triple_hash {
//    inline size_t operator()(const tuple<int,int,int> & v) const {
//        auto &[a,b,c]=v;
//        return a*31*31+b*31+c;
//    }
//};

//struct quad_hash {
//    inline size_t operator()(const tuple<int,int,int,int> & v) const {
//        auto &[a,b,c,d]=v;
//        return a*31*31*31+b*31*31+c*31+d;
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

//int directions[4][2] = {{1,0},{0,1},{-1,0},{0,-1}};
//int directions[8][2] = {{-1,0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};

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

//class UnionFind {                                // OOP style
//private:
//    vi p, rank, setSize;                           // vi p is the key part
//    int numSets;
//public:
//    UnionFind(int N) {
//        p.assign(N, 0); for (int i = 0; i < N; ++i) p[i] = i;
//        rank.assign(N, 0);                           // optional speedup
//        setSize.assign(N, 1);                        // optional feature
//        numSets = N;                                 // optional feature
//    }
//
//    int findSet(int i) { return (p[i] == i) ? i : (p[i] = findSet(p[i])); }
//    bool isSameSet(int i, int j) { return findSet(i) == findSet(j); }
//
//    int numDisjointSets() { return numSets; }      // optional
//    int sizeOfSet(int i) { return setSize[findSet(i)]; } // optional
//
//    void unionSet(int i, int j) {
//        if (isSameSet(i, j)) return;                 // i and j are in same set
//        int x = findSet(i), y = findSet(j);          // find both rep items
//        if (rank[x] > rank[y]) swap(x, y);           // keep x 'shorter' than y
//        p[x] = y;                                    // set x under y
//        if (rank[x] == rank[y]) ++rank[y];           // optional speedup
//        setSize[y] += setSize[x];                    // combine set sizes at y
//        --numSets;                                   // a union reduces numSets
//    }
//};

typedef tuple<int, ll, ll> edge;

const ll INF = 1e18;                             // large enough

class max_flow {
private:
    int V;
    vector<edge> EL;
    vector<vi> AL;
    vi d, last;
    vector<ii> p;

    bool BFS(int s, int t) {                       // find augmenting path
        d.assign(V, -1); d[s] = 0;
        queue<int> q({s});
        p.assign(V, {-1, -1});                       // record BFS sp tree
        while (!q.empty()) {
            int u = q.front(); q.pop();
            if (u == t) break;                         // stop as sink t reached
            for (auto &idx : AL[u]) {                  // explore neighbors of u
                auto &[v, cap, flow] = EL[idx];          // stored in EL[idx]
                if ((cap-flow > 0) && (d[v] == -1))      // positive residual edge
                    d[v] = d[u]+1, q.push(v), p[v] = {u, idx}; // 3 lines in one!
            }
        }
        return d[t] != -1;                           // has an augmenting path
    }

    ll send_one_flow(int s, int t, ll f = INF) {   // send one flow from s->t
        if (s == t) return f;                        // bottleneck edge f found
        auto &[u, idx] = p[t];
        auto &cap = get<1>(EL[idx]), &flow = get<2>(EL[idx]);
        ll pushed = send_one_flow(s, u, min(f, cap-flow));
        flow += pushed;
        auto &rflow = get<2>(EL[idx^1]);             // back edge
        rflow -= pushed;                             // back flow
        return pushed;
    }

    ll DFS(int u, int t, ll f = INF) {             // traverse from s->t
        if ((u == t) || (f == 0)) return f;
        for (int &i = last[u]; i < (int)AL[u].size(); ++i) { // from last edge
            auto &[v, cap, flow] = EL[AL[u][i]];
            if (d[v] != d[u]+1) continue;              // not part of layer graph
            if (ll pushed = DFS(v, t, min(f, cap-flow))) {
                flow += pushed;
                auto &rflow = get<2>(EL[AL[u][i]^1]);     // back edge
                rflow -= pushed;
                return pushed;
            }
        }
        return 0;
    }

public:
    max_flow(int initialV) : V(initialV) {
        EL.clear();
        AL.assign(V, vi());
    }

    // if you are adding a bidirectional edge u<->v with weight w into your
    // flow graph, set directed = false (default value is directed = true)
    void add_edge(int u, int v, ll w, bool directed = true) {
        if (u == v) return;                          // safeguard: no self loop
        EL.emplace_back(v, w, 0);                    // u->v, cap w, flow 0
        AL[u].push_back(EL.size()-1);                // remember this index
        EL.emplace_back(u, directed ? 0 : w, 0);     // back edge
        AL[v].push_back(EL.size()-1);                // remember this index
    }

    ll edmonds_karp(int s, int t) {
        ll mf = 0;                                   // mf stands for max_flow
        while (BFS(s, t)) {                          // an O(V*E^2) algorithm
            ll f = send_one_flow(s, t);                // find and send 1 flow f
            if (f == 0) break;                         // if f == 0, stop
            mf += f;                                   // if f > 0, add to mf
        }
        return mf;
    }

    ll dinic(int s, int t) {
        ll mf = 0;                                   // mf stands for max_flow
        while (BFS(s, t)) {                          // an O(V^2*E) algorithm
            last.assign(V, 0);                         // important speedup
            while (ll f = DFS(s, t))                   // exhaust blocking flow
                mf += f;
        }
        return mf;
    }
};

/**
 * Max flow modelling, baseball elimination variant.
 * https://open.kattis.com/problems/chesscompetition
 * https://open.kattis.com/problems/unfairplay
 */
int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);

    int tc;
    cin >> tc;
    while (tc--) {
        int n;
        cin >> n;
        vvc arr(n, vc(n));
        for (int i=0;i<n;i++) {
            for (int j=0;j<n;j++) {
                cin >> arr[i][j];
            }
        }
        for (int contestant = 0; contestant < n; contestant++) {
            // we are going to do maxflow on this contestant

            // calculate his score if he wins everything
            int bestContestantScore = 0;
            for (auto &c : arr[contestant]) {
                if (c == '1' || c == '.') bestContestantScore+=2;
                else if (c == 'd') bestContestantScore++;
            }
            int matchesCnt = 0;
            for (int i=0;i<n;i++) {
                for (int j=i+1;j<n;j++) {
                    if (i == contestant || j == contestant) continue;
                    if (arr[i][j] != '.') continue;
                    // contestant i vs contestant j
                    // unplayed game
                    matchesCnt++;
                }
            }

            int s = 0; // this is our source node
            int t = 2 + n + matchesCnt - 1;
            max_flow mf(2 + n + matchesCnt);

            // match vertex cnt
            int matchVertex = 1 + n;
            for (int i=0;i<n;i++) {
                for (int j=i+1;j<n;j++) {
                    if (i == contestant || j == contestant) continue;
                    if (arr[i][j] != '.') continue;
                    mf.add_edge(s, matchVertex,2);
                    mf.add_edge(matchVertex, i + 1, 2);
                    mf.add_edge(matchVertex, j + 1, 2);
                    matchVertex++;
                }
            }

            bool flag = true;
            // add players to sink
            for (int i=0;i<n;i++) {
                int currScore = 0;
                for (auto &c : arr[i]) {
                    if (c == '1') {
                        currScore += 2;
                    } else if (c == 'd') {
                        currScore += 1;
                    }
                }
                int w = bestContestantScore - currScore; // points to match the current contestant
                if (w < 0) {
                    flag = false;
                    break;
                }
                mf.add_edge(i + 1, t, w);
            }

            if (flag) {
                ll maxFlowVal = mf.dinic(s, t);
                if (maxFlowVal == matchesCnt * 2) {
                    printf("%d ", contestant + 1);
                }
            }
        }
        printf("\n");
    }
}

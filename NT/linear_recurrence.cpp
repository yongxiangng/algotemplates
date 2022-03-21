#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <ctime>
#pragma GCC optimize("O3")
using namespace std;
using namespace __gnu_pbds;
#define MIN INT32_MIN
#define MAX INT32_MAX
#define LSOne(S) ((S) & (-(S)))
using ll = long long; using ld = long double;
using ii = pair<int, int>; using iii = tuple<int, int, int>;
using vb = vector<bool>; using vc = vector<char>; using vi = vector<int>; using vl = vector<long long>;
using vii = vector<ii>; using vll = vector<ll>;
using vvb = vector<vector<bool>>; using vvc = vector<vector<char>>; using vvi = vector<vector<int>>; using vvl= vector<vector<long long>>;
using Al = unordered_map<int, vi>; using AL = unordered_map<int, vi>;
using ordered_set = tree<int,null_type,less<>,rb_tree_tag,tree_order_statistics_node_update>;
using ordered_set_ii = tree<ii,null_type,less<>,rb_tree_tag,tree_order_statistics_node_update>;
using ordered_multiset =  tree<int,null_type,less_equal<>,rb_tree_tag,tree_order_statistics_node_update>; // doesn't work I think
//const double EPS = 1e-7;

struct pair_hash {
    inline size_t operator()(const pair<int,int> & v) const {
        return v.first*31+v.second;
    }
};

//auto hashF = [](const pair<int,int> & v) {
//    return 1000000000L*(v.first+1000000L)+v.second+1000000L;
//};
//
//auto eq = [](auto &a, auto &b) {return a == b;};
//
struct triple_hash {
    inline size_t operator()(const tuple<int,int,int> & v) const {
        auto &[a,b,c]=v;
        return a*31*31+b*31+c;
    }
};

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

int directions[4][2] = {{1,0},{0,1},{-1,0},{0,-1}};
//int directions[6][2] = {{-1, 0},{-1,-1},{0,-1},{1,0},{1,1},{0,1}};
//int directions[8][2] = {{-1,0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};

//vi S;

//void Kosaraju(int u, int pass, vi &dfs_num, AL &al, AL &alt, vi &temp) {
//    dfs_num[u] = 1;
//    vi &neighbour = (pass == 1) ? al[u] : alt[u];
//    for (auto &v : neighbour) {
//        if (dfs_num[v] == -1) {
//            Kosaraju(v, pass, dfs_num, al, alt, temp);
//        }
//    }
//    S.push_back(u);
//    if (pass == 2) temp.push_back(u);
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

// --- Original fenwick ---

class FenwickTree {                              // index 0 is not used
private:
    vll ft;                                        // internal FT is an array
public:
    FenwickTree(int m) { ft.assign(m+1, 0); }      // create an empty FT

    void build(const vll &f) {
        int m = (int)f.size()-1;                     // note f[0] is always 0
        ft.assign(m+1, 0);
        for (int i = 1; i <= m; ++i) {               // O(m)
            ft[i] += f[i];                             // add this value
            if (i+LSOne(i) <= m)                       // i has parent
                ft[i+LSOne(i)] += ft[i];                 // add to that parent
        }
    }

    FenwickTree(const vll &f) { build(f); }        // create FT based on f

    FenwickTree(int m, const vi &s) {              // create FT based on s
        vll f(m+1, 0);
        for (int i = 0; i < (int)s.size(); ++i)      // do the conversion first
            ++f[s[i]];                                 // in O(n)
        build(f);                                    // in O(m)
    }

    ll rsq(int j) {                                // returns RSQ(1, j)
        ll sum = 0;
        for (; j; j -= LSOne(j))
            sum += ft[j];
        return sum;
    }

    ll rsq(int i, int j) { return rsq(j) - rsq(i-1); } // inc/exclusion

    // updates value of the i-th element by v (v can be +ve/inc or -ve/dec)
    void update(int i, ll v) {
        for (; i < (int)ft.size(); i += LSOne(i))
            ft[i] += v;
    }

    int select(ll k) {                             // O(log m)
        int p = 1;
        while (p*2 < (int)ft.size()) p *= 2;
        int i = 0;
        while (p) {
            if (k > ft[i+p]) {
                k -= ft[i+p];
                i += p;
            }
            p /= 2;
        }
        return i+1;
    }
};

class RUPQ {                                     // RUPQ variant
private:
    FenwickTree ft;                                // internally use PURQ FT
public:
    RUPQ(int m) : ft(FenwickTree(m)) {}
    void range_update(int ui, int uj, ll v) {
        ft.update(ui, v);                            // [ui, ui+1, .., m] +v
        ft.update(uj+1, -v);                         // [uj+1, uj+2, .., m] -v
    }                                              // [ui, ui+1, .., uj] +v
    ll point_query(int i) { return ft.rsq(i); }    // rsq(i) is sufficient
};

class RURQ  {                                    // RURQ variant
private:                                         // needs two helper FTs
    RUPQ rupq;                                     // one RUPQ and
    FenwickTree purq;                              // one PURQ
public:
    RURQ(int m) : rupq(RUPQ(m)), purq(FenwickTree(m)) {} // initialization
    void range_update(int ui, int uj, ll v) {
        rupq.range_update(ui, uj, v);                // [ui, ui+1, .., uj] +v
        purq.update(ui, v*(ui-1));                   // -(ui-1)*v before ui
        purq.update(uj+1, -v*uj);                    // +(uj-ui+1)*v after uj
    }
    ll rsq(int j) {
        return rupq.point_query(j)*j -               // optimistic calculation
               purq.rsq(j);                          // cancelation factor
    }
    ll rsq(int i, int j) { return rsq(j) - rsq(i-1); } // standard
};

// --- Modified fenwick ---

//class FenwickTree {                              // index 0 is not used
//private:
//    unordered_map<long, long> ft;                                        // internal FT is an array
//    long maxVal = 0;                                    // maxVal in the tree
//public:
//    FenwickTree(long m) { maxVal = m; }      // create an empty FT
//
//    void build(const vll &f) {
//        long m = (long)f.size()-1;                     // note f[0] is always 0
//        maxVal = m;
//        for (long i = 1; i <= m; ++i) {               // O(m)
//            ft[i] += f[i];                             // add this value
//            if (i+LSOne(i) <= maxVal)                       // i has parent
//                ft[i+LSOne(i)] += ft[i];                 // add to that parent
//        }
//    }
//
//    FenwickTree(const vll &f) { build(f); }        // create FT based on f
//
//    FenwickTree(long m, const vl &s) {              // create FT based on s
//        vll f(m+1, 0);
//        maxVal = m;
//        for (long i = 0; i < (long)s.size(); ++i)      // do the conversion first
//            ++f[s[i]];                                 // in O(n)
//        build(f);                                    // in O(m)
//    }
//
//    ll rsq(long j) {                                // returns RSQ(1, j)
//        ll sum = 0;
//        for (; j; j -= LSOne(j))
//            sum += ft[j];
//        return sum;
//    }
//
//    ll rsq(long i, long j) { return rsq(j) - rsq(i-1); } // inc/exclusion
//
//    // updates value of the i-th element by v (v can be +ve/inc or -ve/dec)
//    void update(long i, ll v) {
//        for (; i <= maxVal; i += LSOne(i))
//            ft[i] += v;
//    }
//
//    long select(ll k) {                             // O(log m)
//        long p = 1;
//        while (p*2 <= maxVal) p *= 2;
//        long i = 0;
//        while (p) {
//            if (k > ft[i+p]) {
//                k -= ft[i+p];
//                i += p;
//            }
//            p /= 2;
//        }
//        return i+1;
//    }
//};
//
//class RUPQ {                                     // RUPQ variant
//private:
//    FenwickTree ft;                                // internally use PURQ FT
//public:
//    RUPQ(long m) : ft(FenwickTree(m)) {}
//    void range_update(long ui, long uj, ll v) {
//        ft.update(ui, v);                            // [ui, ui+1, .., m] +v
//        ft.update(uj+1, -v);                         // [uj+1, uj+2, .., m] -v
//    }                                              // [ui, ui+1, .., uj] +v
//    ll point_query(long i) { return ft.rsq(i); }    // rsq(i) is sufficient
//};
//
//class RURQ  {                                    // RURQ variant
//private:                                         // needs two helper FTs
//    RUPQ rupq;                                     // one RUPQ and
//    FenwickTree purq;                              // one PURQ
//public:
//    RURQ(long m) : rupq(RUPQ(m)), purq(FenwickTree(m)) {} // initialization
//    void range_update(long ui, long uj, ll v) {
//        rupq.range_update(ui, uj, v);                // [ui, ui+1, .., uj] +v
//        purq.update(ui, v*(ui-1));                   // -(ui-1)*v before ui
//        purq.update(uj+1, -v*uj);                    // +(uj-ui+1)*v after uj
//    }
//    ll rsq(long j) {
//        return rupq.point_query(j)*j -               // optimistic calculation
//               purq.rsq(j);                          // cancelation factor
//    }
//    ll rsq(long i, long j) { return rsq(j) - rsq(i-1); } // standard
//};

//class SegmentTree {                              // OOP style
//private:
//    int n;                                         // n = (int)A.size()
//    vi A, st, lazy;                                // the arrays
//
//    int l(int p) { return  p<<1; }                 // go to left child
//    int r(int p) { return (p<<1)+1; }              // go to right child
//
//    int conquer(int a, int b) {
//        if (a == -1) return b;                       // corner case
//        if (b == -1) return a;
//        return min(a, b);                            // RMQ
//    }
//
//    void build(int p, int L, int R) {              // O(n)
//        if (L == R)
//            st[p] = A[L];                              // base case
//        else {
//            int m = (L+R)/2;
//            build(l(p), L  , m);
//            build(r(p), m+1, R);
//            st[p] = conquer(st[l(p)], st[r(p)]);
//        }
//    }
//
//    void propagate(int p, int L, int R) {
//        if (lazy[p] != -1) {                         // has a lazy flag
//            st[p] = lazy[p];                           // [L..R] has same value
//            if (L != R)                                // not a leaf
//                lazy[l(p)] = lazy[r(p)] = lazy[p];       // propagate downwards
//            else                                       // L == R, a single index
//                A[L] = lazy[p];                          // time to update this
//            lazy[p] = -1;                              // erase lazy flag
//        }
//    }
//
//    int RMQ(int p, int L, int R, int i, int j) {   // O(log n)
//        propagate(p, L, R);                          // lazy propagation
//        if (i > j) return -1;                        // infeasible
//        if ((L >= i) && (R <= j)) return st[p];      // found the segment
//        int m = (L+R)/2;
//        return conquer(RMQ(l(p), L  , m, i          , min(m, j)),
//                       RMQ(r(p), m+1, R, max(i, m+1), j        ));
//    }
//
//    void update(int p, int L, int R, int i, int j, int val) { // O(log n)
//        propagate(p, L, R);                          // lazy propagation
//        if (i > j) return;
//        if ((L >= i) && (R <= j)) {                  // found the segment
//            lazy[p] = val;                             // update this
//            propagate(p, L, R);                        // lazy propagation
//        }
//        else {
//            int m = (L+R)/2;
//            update(l(p), L  , m, i          , min(m, j), val);
//            update(r(p), m+1, R, max(i, m+1), j        , val);
//            int lsubtree = (lazy[l(p)] != -1) ? lazy[l(p)] : st[l(p)];
//            int rsubtree = (lazy[r(p)] != -1) ? lazy[r(p)] : st[r(p)];
//            st[p] = (lsubtree <= rsubtree) ? st[l(p)] : st[r(p)];
//        }
//    }
//
//public:
//    SegmentTree(int sz) : n(sz), st(4*n), lazy(4*n, -1) {}
//
//    SegmentTree(const vi &initialA) : SegmentTree((int)initialA.size()) {
//        A = initialA;
//        build(1, 0, n-1);
//    }
//
//    void update(int i, int j, int val) { update(1, 0, n-1, i, j, val); }
//
//    int RMQ(int i, int j) { return RMQ(1, 0, n-1, i, j); }
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

//typedef tuple<int, ll, ll, ll> edge;
//
//const ll INF = 1e18; // INF = 1e18, not 2^63-1 to avoid overflow
//
//class min_cost_max_flow {
//private:
//  int V;
//  ll total_cost;
//  vector<edge> EL;
//  vector<vi> AL;
//  vll d;
//  vi last, vis;
//
//  bool SPFA(int s, int t) { // SPFA to find augmenting path in residual graph
//    d.assign(V, INF); d[s] = 0; vis[s] = 1;
//    queue<int> q({s});
//    while (!q.empty()) {
//      int u = q.front(); q.pop(); vis[u] = 0;
//      for (auto &idx : AL[u]) {                  // explore neighbors of u
//        auto &[v, cap, flow, cost] = EL[idx];          // stored in EL[idx]
//        if ((cap-flow > 0) && (d[v] > d[u] + cost)) {      // positive residual edge
//          d[v] = d[u]+cost;
//          if(!vis[v]) q.push(v), vis[v] = 1;
//        }
//      }
//    }
//    return d[t] != INF;                           // has an augmenting path
//  }
//
//  ll DFS(int u, int t, ll f = INF) {             // traverse from s->t
//    if ((u == t) || (f == 0)) return f;
//    vis[u] = 1;
//    for (int &i = last[u]; i < (int)AL[u].size(); ++i) { // from last edge
//      auto &[v, cap, flow, cost] = EL[AL[u][i]];
//      if (!vis[v] && d[v] == d[u]+cost) {                      // in current layer graph
//        if (ll pushed = DFS(v, t, min(f, cap-flow))) {
//      total_cost += pushed * cost;
//          flow += pushed;
//          auto &[rv, rcap, rflow, rcost] = EL[AL[u][i]^1]; // back edge
//          rflow -= pushed;
//          vis[u] = 0;
//          return pushed;
//        }
//      }
//    }
//    vis[u] = 0;
//    return 0;
//  }
//
//public:
//  min_cost_max_flow(int initialV) : V(initialV), total_cost(0) {
//    EL.clear();
//    AL.assign(V, vi());
//    vis.assign(V, 0);
//  }
//
//  // if you are adding a bidirectional edge u<->v with weight w into your
//  // flow graph, set directed = false (default value is directed = true)
//  void add_edge(int u, int v, ll w, ll c, bool directed = true) { // VERY IMPORTANT: w is capacity, c is cost
//    if (u == v) return;                          // safeguard: no self loop
//    EL.emplace_back(v, w, 0, c);                    // u->v, cap w, flow 0
//    AL[u].push_back(EL.size()-1);                // remember this index
//    EL.emplace_back(u, 0, 0, -c);     // back edge
//    AL[v].push_back(EL.size()-1);                // remember this index
//    if (!directed) add_edge(v, u, w, c);
//  }
//
//  pair<ll, ll> mcmf(int s, int t) {
//    ll mf = 0;                                   // mf stands for max_flow
//    while (SPFA(s, t)) {                          // an O(V^2*E) algorithm
//      last.assign(V, 0);                         // important speedup
//      while (ll f = DFS(s, t))                   // exhaust blocking flow
//        mf += f;
//    }
//    return {mf, total_cost};
//  }
//
//  vector<edge> getEdges() {
//    return EL;
//  }
//};

long mod(long a, long m) {                          // returns a (mod m)
    return ((a%m) + m) % m;                        // ensure positive answer
}

long globalAns[42][42];

void matPow(vvl &mat, long n, long MOD) {
    if (n == 1) {
        n = mat.size();
        for (int i=0;i<n;i++) {
            for (int j=0;j<n;j++) {
                globalAns[i][j] = mat[i][j];
            }
        }
        return;
    }
    if (n % 2 == 0) {
        matPow(mat, n / 2, MOD);
        n = mat.size();
        vvl ans = vvl(n, vl(n));
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < n; ++k) {
                for (int j = 0; j < n; ++j) {
                    ans[i][j] += globalAns[i][k] * globalAns[k][j];
                    ans[i][j] = mod(ans[i][j], MOD);
                }
            }
        }
        for (int i=0;i<n;i++) {
            for (int j=0;j<n;j++) {
                globalAns[i][j] = ans[i][j];
            }
        }
    } else {
        matPow(mat, n - 1, MOD);
        n = mat.size();
        vvl ans = vvl(n, vl(n));
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < n; ++k) {
                for (int j = 0; j < n; ++j) {
                    ans[i][j] += globalAns[i][k] * mat[k][j];
                    ans[i][j] = mod(ans[i][j], MOD);
                }
            }
        }
        for (int i=0;i<n;i++) {
            for (int j=0;j<n;j++) {
                globalAns[i][j] = ans[i][j];
            }
        }
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);

    int n;
    cin >> n;
    vl coeff(n + 1), init(n);
    for (auto &a : coeff) cin >> a;
    for (auto &a : init) cin >> a;
    int m;
    cin >> m;
    for (int i=0;i<m;i++) {
        long q, MOD; cin >> q >> MOD;

        if (q < n) {
            long ans = init[q];
            if (ans < 0) {
                ans += MOD * ((abs(ans) / MOD) + 1L);
            }
            ans %= MOD;
            cout << ans << endl;
            continue;
        }
        vl colVector(n + 1);
        colVector[0] = 1;
        for (int i=0;i<n;i++) colVector[n-i] = init[i];
        vvl mat(n + 1, vl(n + 1, 0));

        for (int i=0;i<=n;i++) {
            mat[1][i] = coeff[i];
        }

        mat[0][0] = 1;
        for (int i = 2;i<=n;i++) {
            mat[i][i - 1] = 1;
        }

        matPow(mat, q - n + 1, MOD);
        ll ans = 0;
        for (int i=0;i<=n;i++) {
            ans += globalAns[1][i] * colVector[i];
            ans = mod(ans, MOD);
        }

        cout << ans << endl;
    }
}



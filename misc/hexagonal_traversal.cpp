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
//int directions[6][2] = {{-1, 0},{-1,-1},{0,-1},{1,0},{1,1},{0,1}};
//int directions[8][2] = {{-1,0},{-1,1},{0,1},{1,1},{1,0},{1,-1},{0,-1},{-1,-1}};

//vi S;
//vi T;

//void Kosaraju(int u, int pass, int &size, vi &dfs_num, unordered_map<int, unordered_set<int>> &al, unordered_map<int, unordered_set<int>> &alt, vi& temp) {
//    size++;
//    dfs_num[u] = 1;
//    unordered_set<int> &neighbour = (pass == 1) ? al[u] : alt[u];
//    for (auto &v : neighbour) {
//        if (dfs_num[v] == -1) {
//            Kosaraju(v, pass, size, dfs_num, al, alt, temp);
//        }
//    }
//    S.push_back(u);
//    T.push_back(u);
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

//class FenwickTree {                              // index 0 is not used
//private:
//    vll ft;                                        // internal FT is an array
//public:
//    FenwickTree(int m) { ft.assign(m+1, 0); }      // create an empty FT
//
//    void build(const vll &f) {
//        int m = (int)f.size()-1;                     // note f[0] is always 0
//        ft.assign(m+1, 0);
//        for (int i = 1; i <= m; ++i) {               // O(m)
//            ft[i] += f[i];                             // add this value
//            if (i+LSOne(i) <= m)                       // i has parent
//                ft[i+LSOne(i)] += ft[i];                 // add to that parent
//        }
//    }
//
//    FenwickTree(const vll &f) { build(f); }        // create FT based on f
//
//    FenwickTree(int m, const vi &s) {              // create FT based on s
//        vll f(m+1, 0);
//        for (int i = 0; i < (int)s.size(); ++i)      // do the conversion first
//            ++f[s[i]];                                 // in O(n)
//        build(f);                                    // in O(m)
//    }
//
//    ll rsq(int j) {                                // returns RSQ(1, j)
//        ll sum = 0;
//        for (; j; j -= LSOne(j))
//            sum += ft[j];
//        return sum;
//    }
//
//    ll rsq(int i, int j) { return rsq(j) - rsq(i-1); } // inc/exclusion
//
//    // updates value of the i-th element by v (v can be +ve/inc or -ve/dec)
//    void update(int i, ll v) {
//        for (; i < (int)ft.size(); i += LSOne(i))
//            ft[i] += v;
//    }
//
//    int select(ll k) {                             // O(log m)
//        int p = 1;
//        while (p*2 < (int)ft.size()) p *= 2;
//        int i = 0;
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
//    RUPQ(int m) : ft(FenwickTree(m)) {}
//    void range_update(int ui, int uj, ll v) {
//        ft.update(ui, v);                            // [ui, ui+1, .., m] +v
//        ft.update(uj+1, -v);                         // [uj+1, uj+2, .., m] -v
//    }                                              // [ui, ui+1, .., uj] +v
//    ll point_query(int i) { return ft.rsq(i); }    // rsq(i) is sufficient
//};
//
//class RURQ  {                                    // RURQ variant
//private:                                         // needs two helper FTs
//    RUPQ rupq;                                     // one RUPQ and
//    FenwickTree purq;                              // one PURQ
//public:
//    RURQ(int m) : rupq(RUPQ(m)), purq(FenwickTree(m)) {} // initialization
//    void range_update(int ui, int uj, ll v) {
//        rupq.range_update(ui, uj, v);                // [ui, ui+1, .., uj] +v
//        purq.update(ui, v*(ui-1));                   // -(ui-1)*v before ui
//        purq.update(uj+1, -v*uj);                    // +(uj-ui+1)*v after uj
//    }
//    ll rsq(int j) {
//        return rupq.point_query(j)*j -               // optimistic calculation
//               purq.rsq(j);                          // cancelation factor
//    }
//    ll rsq(int i, int j) { return rsq(j) - rsq(i-1); } // standard
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

//typedef tuple<int, ll, ll> edge;
//
//const ll INF = 1e18;                             // large enough
//
//class max_flow {
//private:
//    int V;
//    vector<edge> EL;
//    vector<vi> AL;
//    vi d, last;
//    vector<ii> p;
//
//    bool BFS(int s, int t) {                       // find augmenting path
//        d.assign(V, -1); d[s] = 0;
//        queue<int> q({s});
//        p.assign(V, {-1, -1});                       // record BFS sp tree
//        while (!q.empty()) {
//            int u = q.front(); q.pop();
//            if (u == t) break;                         // stop as sink t reached
//            for (auto &idx : AL[u]) {                  // explore neighbors of u
//                auto &[v, cap, flow] = EL[idx];          // stored in EL[idx]
//                if ((cap-flow > 0) && (d[v] == -1))      // positive residual edge
//                    d[v] = d[u]+1, q.push(v), p[v] = {u, idx}; // 3 lines in one!
//            }
//        }
//        return d[t] != -1;                           // has an augmenting path
//    }
//
//    ll send_one_flow(int s, int t, ll f = INF) {   // send one flow from s->t
//        if (s == t) return f;                        // bottleneck edge f found
//        auto &[u, idx] = p[t];
//        auto &cap = get<1>(EL[idx]), &flow = get<2>(EL[idx]);
//        ll pushed = send_one_flow(s, u, min(f, cap-flow));
//        flow += pushed;
//        auto &rflow = get<2>(EL[idx^1]);             // back edge
//        rflow -= pushed;                             // back flow
//        return pushed;
//    }
//
//    ll DFS(int u, int t, ll f = INF) {             // traverse from s->t
//        if ((u == t) || (f == 0)) return f;
//        for (int &i = last[u]; i < (int)AL[u].size(); ++i) { // from last edge
//            auto &[v, cap, flow] = EL[AL[u][i]];
//            if (d[v] != d[u]+1) continue;              // not part of layer graph
//            if (ll pushed = DFS(v, t, min(f, cap-flow))) {
//                flow += pushed;
//                auto &rflow = get<2>(EL[AL[u][i]^1]);     // back edge
//                rflow -= pushed;
//                return pushed;
//            }
//        }
//        return 0;
//    }
//
//public:
//    max_flow(int initialV) : V(initialV) {
//        EL.clear();
//        AL.assign(V, vi());
//    }
//
//    // if you are adding a bidirectional edge u<->v with weight w into your
//    // flow graph, set directed = false (default value is directed = true)
//    void add_edge(int u, int v, ll w, bool directed = true) {
//        if (u == v) return;                          // safeguard: no self loop
//        EL.emplace_back(v, w, 0);                    // u->v, cap w, flow 0
//        AL[u].push_back(EL.size()-1);                // remember this index
//        EL.emplace_back(u, directed ? 0 : w, 0);     // back edge
//        AL[v].push_back(EL.size()-1);                // remember this index
//    }
//
//    ll edmonds_karp(int s, int t) {
//        ll mf = 0;                                   // mf stands for max_flow
//        while (BFS(s, t)) {                          // an O(V*E^2) algorithm
//            ll f = send_one_flow(s, t);                // find and send 1 flow f
//            if (f == 0) break;                         // if f == 0, stop
//            mf += f;                                   // if f > 0, add to mf
//        }
//        return mf;
//    }
//
//    ll dinic(int s, int t) {
//        ll mf = 0;                                   // mf stands for max_flow
//        while (BFS(s, t)) {                          // an O(V^2*E) algorithm
//            last.assign(V, 0);                         // important speedup
//            while (ll f = DFS(s, t))                   // exhaust blocking flow
//                mf += f;
//        }
//        return mf;
//    }
//};

/**
 * https://open.kattis.com/problems/settlers2
 */
int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);


	int size = 20; // change this to 20 to see better
	vvi arr(size, vi(size, -1));
	
	vi ans;
	unordered_map<int, int> freq;
	
	int currI = size / 2;
	int currJ = size / 2;
	
	arr[currI][currJ] = 1;
	freq[1] = 1;
	ans.push_back(1);
	
	int currDirection = 0;
	int directions[6][2] = {{-1, 0},{-1,-1},{0,-1},{1,0},{1,1},{0,1}};
	int i = 2;
	while (i < size * size / 4) {
		auto &direction = directions[currDirection];
		currI += direction[0]; currJ += direction[1];
		unordered_set<int> neighbours;
		for (auto &direction : directions) {
			int newI = currI + direction[0];
			int newJ = currJ + direction[1];
			if (arr[newI][newJ] != -1) {
				neighbours.insert(arr[newI][newJ]);
			}
		}
		
		int bestVal;
		int bestFreq = MAX;
		for (int val = 1; val <= 5; val++) {
			if (neighbours.count(val) == 0) {
				if (freq[val] < bestFreq) {
					bestFreq = freq[val];
					bestVal = val;
				}
			}
		}
		
		arr[currI][currJ] = bestVal;
		ans.push_back(bestVal);
		freq[bestVal]++;
		
		
		while (arr[currI + directions[(currDirection + 1) % 6][0]][currJ + directions[(currDirection + 1) % 6][1]] == -1) {
			currDirection = (currDirection + 1) % 6;
		}
		
		i++;
	}
	
	
	for (auto &a : arr) {
		for (auto &b : a) {
			string temp;
			temp = to_string(b);
			while (temp.size() < 3) {
				temp = " " + temp;
			}
			cout << temp;
		}
		cout << endl;
	}
	
	//int tc;
	//cin >> tc;
	//while (tc--) {
		//int q;
		//cin >> q;
		//cout << ans[q - 1] << endl;
	//}
}

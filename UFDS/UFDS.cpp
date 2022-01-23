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
typedef vector<ll> vll;
typedef vector<vector<char>> vvc;
typedef vector<vector<int>> vvi;
typedef unordered_map<int, vi> Al;
typedef unordered_map<int, vi> AL;
typedef tree<int,null_type,less<>,rb_tree_tag,tree_order_statistics_node_update> ordered_set;
typedef tree<ii,null_type,less<>,rb_tree_tag,tree_order_statistics_node_update> ordered_set_ii;
typedef tree<int,null_type,less_equal<>,rb_tree_tag,tree_order_statistics_node_update> ordered_multiset;
const double EPS = 1e-7;

/**
 * Stolen, as usual. https://github.com/stevenhalim/cpbook-code/blob/master/ch2/ourown/unionfind_ds.cpp
 */
class UnionFind {                                // OOP style
private:
    vi p, rank, setSize;                           // vi p is the key part
    int numSets;
public:
    UnionFind(int N) {
        p.assign(N, 0); for (int i = 0; i < N; ++i) p[i] = i;
        rank.assign(N, 0);                           // optional speedup
        setSize.assign(N, 1);                        // optional feature
        numSets = N;                                 // optional feature
    }

    int findSet(int i) { return (p[i] == i) ? i : (p[i] = findSet(p[i])); }
    bool isSameSet(int i, int j) { return findSet(i) == findSet(j); }

    int numDisjointSets() { return numSets; }      // optional
    int sizeOfSet(int i) { return setSize[findSet(i)]; } // optional

    void unionSet(int i, int j) {
        if (isSameSet(i, j)) return;                 // i and j are in same set
        int x = findSet(i), y = findSet(j);          // find both rep items
        if (rank[x] > rank[y]) swap(x, y);           // keep x 'shorter' than y
        p[x] = y;                                    // set x under y
        if (rank[x] == rank[y]) ++rank[y];           // optional speedup
        setSize[y] += setSize[x];                    // combine set sizes at y
        --numSets;                                   // a union reduces numSets
    }
};


/**
 * Sample usage
 *
 * https://open.kattis.com/problems/swaptosort
 */


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
//
//int main() {
//    ios::sync_with_stdio(false);
//    cin.tie(NULL); cout.tie(NULL);
//
//    // ufds
//    Al al;
//    int n, k;
//    cin >> n >> k;
//    UnionFind UF(n);
//    for (int i=0;i<k;i++) {
//        int u, v;
//        cin >> u >> v;
//        u--;v--;
//        UF.unionSet(u, v);
//    }
//    for (int i=0;i<n;i++) {
//        if (!UF.isSameSet(i, n - 1 - i)) {
//            cout << "No" << endl;
//            return 0;
//        }
//    }
//    cout << "Yes" << endl;
//}

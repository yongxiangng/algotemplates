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
const double EPS = 1e-9;

/**
 * https://open.kattis.com/problems/powerstrings
 *
 * Part of KMP. Will update with search code if I encounter such questions.
 */

string P;
vi b;
int m;                                        // n = |T|, m = |P|
void kmpPreprocess() {                           // call this first
    int i = 0, j = -1; b[0] = -1;                  // starting values
    while (i < m) {                                // pre-process P
        while ((j >= 0) && (P[i] != P[j])) j = b[j]; // different, reset j
        ++i; ++j;                                    // same, advance both
        b[i] = j;
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);

    while (true) {
        getline(cin, P);
        if (P == ".") break;
        m = P.size();
        b.resize(m + 1);
        kmpPreprocess();
        cout << m / (m - b[m]) << endl;
    }
}

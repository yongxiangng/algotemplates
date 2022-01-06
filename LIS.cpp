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

void print_LIS(int i, vi &A, vi &p) {                             // backtracking routine
    if (p[i] == -1) { printf("%d", A[i]); return; } // base case
    print_LIS(p[i], A, p);                                  // backtrack
    printf(" %d", A[i]);
}

int LIS(vi &arr) {
    int n = arr.size();
    vi p;

    int k = 0, lis_end = 0;
    vi L(n, 0), L_id(n, 0);
    p.assign(n, -1);

    for (int i = 0; i < n; ++i) {                    // O(n)
        // use lower_bound if you want strictly increasing
        // use upper_bound if you want non strictly increasing
        int pos = lower_bound(L.begin(), L.begin()+k, arr[i]) - L.begin();
        L[pos] = arr[i];                             // greedily overwrite this
        L_id[pos] = i;                               // remember the index too
        p[i] = pos ? L_id[pos-1] : -1;               // predecessor info
        if (pos == k) {                              // can extend LIS?
            k = pos+1;                               // k = longer LIS by +1
            lis_end = i;                             // keep best ending i
        }
    }

    print_LIS(lis_end, arr, p);

    return k;
}

int main() {
    vi arr({-7, 10, 9, 2, 3, 8, 8, 1, 2, 3, 4, 99});
    int size = LIS(arr);
    cout << endl << "The size is " << size << endl;
}

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

/**
 * Square Root Decomposition Example. CP4 Book 2, pg 489.
 * 
 * https://www.youtube.com/watch?v=BJhzd_VG61k
 */

int Arr[200010];
int bucket[450][450]; // If bucket[b][a] = c it means that we're adding c to all k = a mod b

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);
    
    int n, q;
    cin >> n >> q;
    
    memset(Arr, 0, sizeof Arr);
    memset(bucket, 0, sizeof bucket);
        
	int bound = (int) sqrt(n) + 1;
	
    for (int i=0;i<q;i++) {
		int t;
		cin >> t;
		if (t == 1) {
			int a,b,c;
			cin >> a >> b >> c;
			if (b > bound) {
				// then we just manually add it into the array
				while (a <= n) {
					Arr[a] += c;
					a += b;
				}
			} else {
				bucket[b][a] += c;
			}
		} else {
			// t == 2
			int d;
			cin >> d;
			long ans = 0;
			ans += Arr[d];
			for (int i=1;i<=bound;i++) { // iterate thru each bucket
				ans += bucket[i][d % i]; // d % i (mod i) -> d is incremented by this bucket
			}
			cout << ans << endl;
		}
		
	}
}

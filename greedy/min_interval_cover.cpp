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
 * https://www.geeksforgeeks.org/minimum-number-of-intervals-to-cover-the-target-interval/
 * 
 * https://open.kattis.com/problems/intervalcover
 */
int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);
    double a, b;
    while (cin >> a >> b) {
		int n;
		cin >> n;
		vector<pair<double, double>> arr(n);
		map<pair<double, double>, int> toIdx;
		for (int i=0;i<n;i++) {
			cin >> arr[i].first >> arr[i].second;
			toIdx[{arr[i].first, arr[i].second}] = i;
		}
		sort(arr.begin(), arr.end());
		arr.push_back({1e10, 1e10});
		vi ans;
		double bestStart;
		double bestEnd;
		double start = a;
		double end = a - 1;
		long cnt = 0;
		for (int i=0;i<n+1;) {
			if (arr[i].first <= start) {
				end = max(end, arr[i++].second);
				if (arr[i - 1].second == end) {
					bestStart = arr[i - 1].first;
					bestEnd = arr[i - 1].second;
				}
			} else {
				start = end;
				++cnt;
				ans.push_back(toIdx[{bestStart, bestEnd}]);
				if (arr[i].first > end || end >= b) {
					break;
				}
			}
		}
		
		if (end < b) {
			printf("impossible\n");
		} else {
			printf("%ld\n", cnt);
			for (auto &i : ans) {
				printf("%d ", i);
			}
			printf("\n");
		}
		
	}
}


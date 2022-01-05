#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
typedef pair<int, int> ii;

/**
* Min max variant
*/
void SlidingWindow(int A[], int n, int K) {
    deque<ii> window;                                       // we sort window in ascending order
    for (int i=0;i<n;++i) {
        while (!window.empty() && window.back().first >= A[i]) {
            window.pop_back();                              // keep window sorted
        }

        window.push_back({A[i], i});               // store value by index

        while (window.front().second <= i - K) {
            window.pop_front();
        }
        if (i+1 >= K) {
            printf("%d\n", window.front().first);     // min for this window
        }
    }
}

#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace std;
using namespace __gnu_pbds;

int main() {
    int mask = 18;
    for (int subset = mask; subset; subset = (mask & (subset - 1))) {
        cout << subset << endl;
    }
}

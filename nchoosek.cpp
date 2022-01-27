// nCk implementation by Rezwan
#include <bits/stdc++.h>
using namespace std;

int main() {
  int n = 6;
  int k = 3; // C(6, 3)
  for (int mask = (1<<k)-1; mask < (1<<n); ) {
    cout << mask << '\n'; // only masks with k bits on here
    int x = mask & -mask; // LSOne(mask)
    int y = mask + x;
    // mask = ((mask & ~y) / x >> 1) | y; // same # of on bits, involve /
    mask = ((mask & ~y) >> (1 + __builtin_ctz(x))) | y; // same # of on bits
  }
  return 0;
}


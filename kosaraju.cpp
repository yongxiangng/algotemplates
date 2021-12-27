// finds SCC

#include <bits/stdc++.h>
using namespace std;
typedef vector<int> vi;
typedef unordered_map<int, vi> AL;

vi S;

void Kosaraju(int u, int pass, int &size, vi &dfs_num, AL &al, AL &alt) {
    size++;
    dfs_num[u] = 1;
    vi &neighbour = (pass == 1) ? al[u] : alt[u];
    for (auto &v : neighbour) {
        if (dfs_num[v] == -1) {
            Kosaraju(v, pass, size, dfs_num, al, alt);
        }
    }
    S.push_back(u);
}

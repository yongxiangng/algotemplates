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

int main() {
    int n;
    AL al, alt;

    S.clear();
    vi dfs_num;
    dfs_num.assign(n, -1);
    for (int i=0;i<n;i++) {
        if (dfs_num[i] == -1) {
            int size = 0;
            Kosaraju(i, 1, size, dfs_num, al, alt);
        }
    }

    int maxSize = 0;
    dfs_num.assign(n, -1);
    for (int i=n-1;i>=0;i--) {
        if (dfs_num[S[i]] == -1) {
            int size = 0;
            Kosaraju(S[i], 2, size, dfs_num, al, alt);
            maxSize = max(maxSize, size);
        }
    }
}

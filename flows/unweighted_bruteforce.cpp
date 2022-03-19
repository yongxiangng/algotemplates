// run with MCM(0) I believe

int MCM(int bitmask) {
    if (bitmask == (1 << N) - 1) return 0;
    int &ans = memo[bitmask];
    if (ans != -1) return ans;

    int p1, p2;
    for (p1 = 0; p1 < N; ++p1) {
        if (!(bitmask & (1<<p1))) {
            break;
        }
    }

    ans = MCM(bitmask | (1<<p1));

    for (p2 = 0; p2 < N; ++p2) {
        if (AM[p1][p2] && (p2 != p1) && !(bitmask & (1<<p2))) {
            ans = max(ans, 1 + MCM(bitmask | (1<<p1) | (1<<p2)));
        }
    }
    return ans;
}


// bruteforce wMCM steven

ii wMCM(int mask) {
    if (mask == (1<<N)-1) return ii(0, 0);
    if (memo[mask] != ii(-1, -1)) return memo[mask];
    int p1, p2;
    for (p1 = 0; p1 < N; p1++) if (!(mask & (1<<p1))) break;
    ii ans = wMCM(mask | (1<<p1)); // p1 unmatched
    for (p2 = p1+1; p2 < N; p2++)
        if (!(mask & (1<<p2)) && cost[p1][p2]) {
            ii nxt = wMCM(mask | (1<<p1) | (1<<p2)); // match p1-p2
            nxt.first += 2; nxt.second += cost[p1][p2];
            if ((nxt.first > ans.first) || // more # matching
                    ((nxt.first == ans.first) && // or equal # matching
                     (nxt.second < ans.second))) // with smaller cost
                ans = nxt;
        }
    return memo[mask] = ans;
}

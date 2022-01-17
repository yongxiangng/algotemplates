#include <bits/stdc++.h>
using namespace std;

#define LSOne(S) ((S) & -(S))                    // the key operation

typedef long long ll;                            // for extra flexibility
typedef vector<ll> vll;
typedef vector<int> vi;

/**
 * Avoid if possible, this is 4 times slower than standard implementation
 */
class FenwickTree {                              // index 0 is not used
private:
    unordered_map<long, long> ft;                                        // internal FT is an array
    long maxVal = 0;                                    // maxVal in the tree
public:
    FenwickTree(long m) { maxVal = m; }      // create an empty FT

    void build(const vll &f) {
        long m = (long)f.size()-1;                     // note f[0] is always 0
        maxVal = m;
        for (long i = 1; i <= m; ++i) {               // O(m)
            ft[i] += f[i];                             // add this value
            if (i+LSOne(i) <= maxVal)                       // i has parent
                ft[i+LSOne(i)] += ft[i];                 // add to that parent
        }
    }

    FenwickTree(const vll &f) { build(f); }        // create FT based on f

    FenwickTree(long m, const vl &s) {              // create FT based on s
        vll f(m+1, 0);
        maxVal = m;
        for (long i = 0; i < (long)s.size(); ++i)      // do the conversion first
            ++f[s[i]];                                 // in O(n)
        build(f);                                    // in O(m)
    }

    ll rsq(long j) {                                // returns RSQ(1, j)
        ll sum = 0;
        for (; j; j -= LSOne(j))
            sum += ft[j];
        return sum;
    }

    ll rsq(long i, long j) { return rsq(j) - rsq(i-1); } // inc/exclusion

    // updates value of the i-th element by v (v can be +ve/inc or -ve/dec)
    void update(long i, ll v) {
        for (; i <= maxVal; i += LSOne(i))
            ft[i] += v;
    }

    long select(ll k) {                             // O(log m)
        long p = 1;
        while (p*2 <= maxVal) p *= 2;
        long i = 0;
        while (p) {
            if (k > ft[i+p]) {
                k -= ft[i+p];
                i += p;
            }
            p /= 2;
        }
        return i+1;
    }
};

class RUPQ {                                     // RUPQ variant
private:
    FenwickTree ft;                                // internally use PURQ FT
public:
    RUPQ(long m) : ft(FenwickTree(m)) {}
    void range_update(long ui, long uj, ll v) {
        ft.update(ui, v);                            // [ui, ui+1, .., m] +v
        ft.update(uj+1, -v);                         // [uj+1, uj+2, .., m] -v
    }                                              // [ui, ui+1, .., uj] +v
    ll point_query(long i) { return ft.rsq(i); }    // rsq(i) is sufficient
};

class RURQ  {                                    // RURQ variant
private:                                         // needs two helper FTs
    RUPQ rupq;                                     // one RUPQ and
    FenwickTree purq;                              // one PURQ
public:
    RURQ(long m) : rupq(RUPQ(m)), purq(FenwickTree(m)) {} // initialization
    void range_update(long ui, long uj, ll v) {
        rupq.range_update(ui, uj, v);                // [ui, ui+1, .., uj] +v
        purq.update(ui, v*(ui-1));                   // -(ui-1)*v before ui
        purq.update(uj+1, -v*uj);                    // +(uj-ui+1)*v after uj
    }
    ll rsq(long j) {
        return rupq.point_query(j)*j -               // optimistic calculation
               purq.rsq(j);                          // cancelation factor
    }
    ll rsq(long i, long j) { return rsq(j) - rsq(i-1); } // standard
};

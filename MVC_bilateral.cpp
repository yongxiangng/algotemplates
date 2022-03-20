// https://open.kattis.com/problems/bilateral

vi match, vis;
Al al;

int Aug(int L) {
	if (vis[L]) return 0;
	vis[L] = 1;
	for (auto &R : al[L]) {
		if ((match[R] == -1) || Aug(match[R])) {
			match[R] = L;
			return 1;
		}
	}
	return 0;
}

int rig(int L) {
	if (vis[L]) return 0;
	vis[L] = 1;
	for (auto &R : al[L]) {
		if (match[R] == L) return 1;
	}
	
	for (auto &R : al[L]) {
		if (rig(match[R])) {
			match[R] = L;
			return 1;
		}
	}
	return 0;
}

void getZ(int L, unordered_set<int> &z) {
	z.insert(L);
	if (vis[L]) return;
	vis[L] = 1;
	for (auto &R : al[L]) {
		if ((match[R] != -1) && match[R] != L) {
			z.insert(R);
			getZ(match[R], z);
		}
	}
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);

	unordered_set<int> l;
	unordered_set<int> r;

	// friend has id of 9
	int m; cin >> m;
	for (int i=0;i<m;i++) {
		int a, b; cin >> a >> b;
		a -= 1000;
		b -= 1000;
		al[a].push_back(b);
		l.insert(a);
		r.insert(b);
	}

	int VLeft = 1000;
	int V = 2000;
	
	unordered_set<int> freeV;
	for (int L = 0; L < VLeft; ++L) {
		freeV.insert(L);
	}
	
	match.assign(V, -1);
	
	int MCBM = 0;
	
	for (int L = 0; L < VLeft; ++L) {
		vi candidates;
		for (auto &R : al[L]) {
			if (match[R] == -1) {
				candidates.push_back(R);
			}
		}
		if (candidates.size() > 0) {
			++MCBM;
			freeV.erase(L);
			int a = rand() % (int) candidates.size();
			match[candidates[a]] = L;
		}
	}
	
	for (auto &f : freeV) {
		vis.assign(VLeft, 0);
		MCBM += Aug(f);
	}
	
	// try to get your friend in
	vis.assign(VLeft, 0);
	rig(9);

	cout << MCBM << endl;
	
	unordered_set<int> u;
	unordered_set<int> z;
	unordered_set<int> mvc;
	unordered_set<int> matched;
	
	for (int i=0;i<V;i++) {
		if (match[i] != -1) {
			matched.insert(i);
			matched.insert(match[i]);
		}
	}
	
	// find u
	for (auto &a : l) {
		if (matched.count(a) == 0) {
			u.insert(a);
		}
	}
	
	// find z
	for (auto &L : u) {
		vis.assign(VLeft, 0);
		getZ(L, z);
	}
	
	// take union
	// l not in z
	for (auto &a : l) {
		if (z.count(a) == 0) {
			mvc.insert(a);
		}
	}
	
	// r in z
	for (auto &a : r) {
		if (z.count(a) != 0) {
			mvc.insert(a);
		}
	}

	for (auto &a : mvc) {
		cout << a + 1000 << endl;
	}
}


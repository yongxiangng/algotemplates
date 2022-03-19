// https://open.kattis.com/problems/halloweenloot

int offset = 10010;
int mem[101][20050];

// returns abs(actual - target)
int findClosest(int i, vi &arr, int target) {
    if (i == arr.size()) return abs(target);

    int &ans = mem[i][target + offset];
    if (ans != -1) return ans;

    int a = findClosest(i + 1, arr, target - arr[i]); // take current
    int b = findClosest(i + 1, arr, target); // don't take current

    return ans = min(a, b);
}

void printClosest(int i, vi &arr, int target) {
    if (i == arr.size()) {
        cout << endl;
        return;
    }

    if (findClosest(i + 1, arr, target - arr[i]) == mem[i][target+offset]) {
        cout << "A";
        printClosest(i + 1, arr, target - arr[i]);
    } else {
        cout << "B";
        printClosest(i + 1, arr, target);
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL); cout.tie(NULL);

    memset(mem, -1, sizeof mem);
	int n;
	cin >> n;
	vi arr1(n), arr2(n);
	for (auto &a : arr1) cin >> a;
	for (auto &a : arr2) cin >> a;
    int sum = 0;
	vi arr(n);
    for (int i=0;i<n;i++) arr[i] = arr1[i] + arr2[i];
    for (int i=0;i<n;i++) sum += arr2[i];

    findClosest(0, arr,sum);
	printClosest(0, arr, sum);
}


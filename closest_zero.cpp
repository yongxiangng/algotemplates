#include <bits/stdc++.h>
using namespace std;

#define arrSize 51
#define maxSum 201
#define MAX 100
#define inf 999999

// Variable to store states of dp
int dp[arrSize][maxSum];
bool visit[arrSize][maxSum];

// Function to return the number closer to integer s
int RetClose(int a, int b, int s)
{
	if (abs(a - s) < abs(b - s))
		return a;
	else
		return b;
}

// To find the sum closest to zero
// Since sum can be negative, we will add MAX
// to it to make it positive
int MinDiff(int i, int sum, int arr[], int n)
{

	// Base cases
	if (i == n)
		return 0;
	// Checks if a state is already solved
	if (visit[i][sum + MAX])
		return dp[i][sum + MAX];
	visit[i][sum + MAX] = 1;

	// Recurrence relation
	dp[i][sum + MAX] = RetClose(arr[i] +
						MinDiff(i + 1, sum + arr[i], arr, n),
						MinDiff(i + 1, sum, arr, n), -1 * sum);

	// Returning the value
	return dp[i][sum + MAX];
}

// Function to calculate the closest sum value
void FindClose(int arr[],int n)
{
	int ans=inf;

	// Calculate the Closest value for every
	// subarray arr[i-1:n]
	for (int i = 1; i <= n; i++)
		ans = RetClose(arr[i - 1] +
				MinDiff(i, arr[i - 1], arr, n), ans, 0);

	cout<<ans<<endl;
}

// Driver function
int main()
{
	// Input array
	int arr[] = { 25, -9, -10, -4, -7, -33 };
	int n = sizeof(arr) / sizeof(int);
	
	FindClose(arr,n);
	return 0;
}


﻿#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <map>
using namespace std;

struct ListNode 
{
	int val;
	ListNode *next;
	ListNode() : val(0), next(nullptr) {}
	ListNode(int x) : val(x), next(nullptr) {}
	ListNode(int x, ListNode *next) : val(x), next(next) {}
};

struct TreeNode {
	int val;
	TreeNode *left;
	TreeNode *right;
	TreeNode() : val(0), left(nullptr), right(nullptr) {}
	TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}
	TreeNode(int x, TreeNode *left, TreeNode *right) : val(x), left(left), right(right) {}
	
};

static int max(int a, int b)
{
	return a >= b ? a : b;
}

static int min(int a, int b)
{
	return a <= b ? a : b;
}

void printVectorString(vector<string>&& S)
{
	for (auto s : S)
	{
		cout << s << endl;
	}

	cout << '\n' << endl;
}

void printDoubleDPTable(vector<int> T1, vector<int> T2)
{
	if (T1.size() != T2.size()) return;

	int n = T1.size();

	for (int i = 0; i < n; ++i)
	{
		cout << T1[i] << " * " << T2[i] << endl;
	}
}

void printDoubleDPTable(vector<vector<int>> T1, vector<vector<int>> T2)
{
	if (T1.empty() || T2.empty()) return;
	if (T1.size() != T2.size() || T1[0].size() != T2[0].size()) return;

	int n = T1.size(), m = T1[0].size();

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			cout << T1[i][j] << " * " << T2[i][j] << '\t';
		}

		cout << '\n';
	}

	cout << endl;
}

void printDPTable(vector<int>& T)
{
	int n = T.size();

	for (int i = 0; i < n; ++i)
		cout << i << "\t";
	cout << "\n" << endl;

	for (auto t : T)
		cout << t << "\t";

	cout << '\n' << endl;
	cout << '\n' << endl;
}

void printDPTable(vector<bool>& T)
{
	int n = T.size();

	for (int i = 0; i < n; ++i)
		cout << i << "\t";
	cout << "\n" << endl;

	for (auto t : T)
	{
		string s;
		if (t)
			s = "true";
		else
			s = "false";

		cout << s << "\t";
	}

	cout << '\n' << endl;
	cout << '\n' << endl;
}

void printDPTable(vector<vector<bool>>& T)
{
	int len = T[0].size();

	cout << "\t";
	for (int i = 0; i < len; ++i)
		cout << i << "\t";
	cout << '\n' << endl;

	int col = 0;
	for (auto a : T)
	{
		cout << col << "\t";
		++col;

		for (auto b : a)
		{
			cout << b << "\t";
		}
		cout << '\n' << endl;
	}
	cout << '\n';
}

void printDPTable(vector<vector<int>>& T)
{
	int len = T[0].size();

	cout << "\t";
	for (int i = 0; i < len; ++i)
		cout << i << "\t";
	cout << '\n' << endl;

	int col = 0;
	for (auto a : T)
	{
		cout << col << "\t";
		++col;

		for (auto b : a)
		{
			cout << b << "\t";
		}
		cout << '\n' << endl;
	}
	cout << '\n';
}

int binarySearch(vector<int>& nums, int target)
{
	int left = 0, right = nums.size() - 1, length = nums.size();
	int mid;

	while (left <= right)
	{
		mid = (left + right) / 2;
		if (nums[mid] == target) return mid;
		else if (nums[mid] < target) left = mid + 1;
		else if (nums[mid] > target) right = mid - 1;
	}

	return -1;
}

bool canDevide3(int n)
{
	int sum = 0;
	while (n != 0)
	{
		sum += (n % 10);
		n /= 10;
	}

	return sum % 3 == 0;
}

#pragma region Sort

// 快排
class QuickSortSolution
{
public:
	vector<int>& Sort(vector<int>& nums)
	{
		//QuickSort(0, nums.size() - 1, nums);
		QuickSortTailRecursion(0, nums.size() - 1, nums);

		return nums;
	}

	void QuickSort(int left, int right, vector<int>& nums)
	{
		if (left < right)
		{
			int p = Partition(left, right, nums);
			QuickSort(left, p - 1, nums);
			QuickSort(p + 1, right, nums);
		}
	}

	// Space: O(n) -> O(lgn)
	void QuickSortTailRecursion(int left, int right, vector<int>& nums)
	{
		while (left < right)
		{
			int p = Partition(left, right, nums);
			if (p - left < right - p)
			{
				QuickSortTailRecursion(left, p - 1, nums);
				left = p + 1;
			}
			else
			{
				QuickSortTailRecursion(p + 1, right, nums);
				right = p - 1;
			}
		}
	}

	int Partition(int left, int right, vector<int>& nums)
	{
		int pivot = nums[right], index = left - 1;

		for (int i = left; i < right; ++i)
		{
			if (nums[i] <= pivot)
			{
				++index;
				int temp = nums[i];
				nums[i] = nums[index];
				nums[index] = temp;
			}
		}

		nums[right] = nums[index + 1];
		nums[index + 1] = pivot;

		return index + 1;
	}
};

// 归并
class MergeSortSolution
{
public:
	vector<int> Sort(vector<int>& nums)
	{

	}
};

// 堆排序
class HeapSortSolution
{
public:
	vector<int> Sort(vector<int>& nums)
	{

	}

	// MAX-HEAPIFY(A, i)
	// BUILD-HEAP(A)
	// HEAP-SORT(A)
};

// 中位数

// LRU

#pragma endregion

#pragma region DP

#pragma region Stock Trade

// 121
int maxProfitForKEquals1(vector<int>& prices)
{
	if (prices.size() == 0) return 0;

	int n = prices.size();
	int T_i0 = 0, T_i1 = INTPTR_MIN;
	for (int i = 0; i < n; ++i)
	{
		T_i0 = max(T_i0, T_i1 + prices[i]);
		T_i1 = max(T_i1, -prices[i]);
	}

	return T_i0;
}

// 123
int maxProfitForKEquals2(vector<int>& prices)
{
	if (prices.size() == 0) return 0;

	int n = prices.size();
	int T_i20 = 0, T_i21 = INT32_MIN, T_i10 = 0, T_i11 = INT32_MIN;
	for (int i = 0; i < n; ++i)
	{
		T_i20 = max(T_i20, T_i21 + prices[i]);
		T_i21 = max(T_i21, T_i10 - prices[i]);
		T_i10 = max(T_i10, T_i11 + prices[i]);
		T_i11 = max(T_i11, -prices[i]);

		//cout << "Day " << i+1 << "   Price " << prices[i] << endl;
		//cout << "T_i20:	" << T_i20 << "    T_i21: " << T_i21 << "    T_i10: " << T_i10 << "    T_i11: " << T_i11 << '\n' << endl;
	}

	return T_i20;
}

int maxProfitForKInfinity(vector<int>& prices)
{
	if (prices.size() == 0) return 0;

	int n = prices.size();
	int T_ik0 = 0, T_ik1 = INT32_MIN;
	for (int i = 0; i < n; ++i)
	{
		T_ik0 = max(T_ik0, T_ik1 + prices[i]);
		T_ik1 = max(T_ik1, T_ik0 - prices[i]);
	}

	return T_ik0;
}

// 188
int maxProfitForKArbitrary(int k, vector<int>& prices)
{
	if (k >= prices.size() >> 1)
	{
		return maxProfitForKInfinity(prices);
	}

	vector<int> T_ik0(k + 1, 0);
	vector<int> T_ik1(k + 1, INT32_MIN);

	for (auto price : prices)
	{
		for (int j = k; j > 0; --j)
		{
			T_ik0[j] = max(T_ik0[j], T_ik1[j] + price);
			T_ik1[j] = max(T_ik1[j], T_ik0[j - 1] - price);
		}
	}

	return T_ik0[k];
}

// 309
int maxProfitForCooldown(vector<int>& prices)
{
	if (prices.size() == 0) return 0;

	int n = prices.size();
	int T_ik0 = 0, T_ik0_pre = 0, T_ik1 = INT32_MIN, T_ik0_old = 0;
	for (int i = 0; i < n; ++i)
	{
		T_ik0_pre = T_ik0_old;
		T_ik0_old = T_ik0;
		T_ik0 = max(T_ik0, T_ik1 + prices[i]);
		T_ik1 = max(T_ik1, T_ik0_pre - prices[i]);

		//cout << "Day " << i+1 << "   Price " << prices[i] << endl;
		//cout << "T_ik0_pre: " << T_ik0_pre << "		T_ik0_old: " << T_ik0_old << "	 T_ik0: " << T_ik0 << "		T_ik1: " << T_ik1 << '\n' << endl;
	}

	return T_ik0;
}

// 714
int maxProfitWithFee(vector<int>& prices, int fee)
{
	if (prices.size() == 0) return 0;

	int n = prices.size();
	int T_ik0 = 0, T_ik1 = INT32_MIN;
	for (int i = 0; i < n; ++i)
	{
		T_ik0 = max(T_ik0, T_ik1 + prices[i]);
		T_ik1 = max(T_ik1, T_ik0 - prices[i] - fee);
	}

	return T_ik0;
}

#pragma endregion

#pragma region Single

#pragma region Rely on O(1) Sub Problem

// 746
int minCostClimbingStairs(vector<int>& cost)
{
	int n = cost.size();
	vector<int> T(n + 1, 0);
	T[1] = cost[0];

	for (int i = 2; i <= n; ++i)
	{
		T[i] = min(T[i - 1], T[i - 2]) + cost[i - 1];
	}

	return min(T[n], T[n - 1]);
}

// 1137 
int tribonacci(int n)
{
	if (n == 0) return 0;
	else if (n == 1 || n == 2) return 1;

	vector<int> T(n + 1, 0);
	T[1] = T[2] = 1;

	for (int i = 3; i <= n; ++i)
	{
		T[i] = T[i - 2] + T[i - 1] + T[i - 3];
	}

	return T[n];
}

// 392
bool isSubsequence(string s, string t)
{
	int index = 0;

	for (auto c : t)
	{
		if (c == s[index])
			++index;
	}

	return index == s.size() ? true : false;
}

// 70
int climbStairs(int n)
{
	if (n == 0) return 0;
	if (n == 1) return 1;
	if (n == 2) return 2;

	int target = 2, target_before1 = 1, target_before2 = 0;

	for (int i = 2; i < n; ++i)
	{
		target_before2 = target_before1;
		target_before1 = target;
		target = target_before1 + target_before2;
	}

	return target;
}

// 53
// O(n^2) 优化到 O(n)
int maxSubArray(vector<int>& nums)
{
	int length = nums.size();
	vector<int> T(length, 0);
	int m = T[0] = nums[0];

	for (int i = 1; i < length; ++i)
	{
		T[i] = T[i - 1] > 0 ? T[i - 1] + nums[i] : nums[i];
		m = max(m, T[i]);
	}

	return m;
}

// 62
int uniquePaths1(int m, int n)
{
	vector<vector<int>> T(m, vector<int>(n, 1));

	for (int i = 1; i < m; ++i)
	{
		for (int j = 1; j < n; ++j)
		{
			T[i][j] = T[i - 1][j] + T[i][j - 1];
		}
	}

	return T[m - 1][n - 1];
}

// 63
int uniquePaths2(vector<vector<int>>& obstacleGrid)
{
	int m = obstacleGrid.size(), n = obstacleGrid[0].size();
	vector<vector<int>> T(m + 1, vector<int>(n + 1, 0));
	T[1][0] = 1;

	for (int i = 1; i <= m; ++i)
	{
		for (int j = 1; j <= n; ++j)
		{
			if (obstacleGrid[i - 1][j - 1] == 0)
				T[i][j] = T[i - 1][j] + T[i][j - 1];
		}
	}

	return T[m][n];
}

// 64
int minPathSum(vector<vector<int>>& grid)
{
	int m = grid.size(), n = grid[0].size();
	vector<vector<int>> T(m + 1, vector<int>(n + 1, INT16_MAX));
	T[1][0] = 0;

	for (int i = 1; i <= m; ++i)
	{
		for (int j = 1; j <= n; ++j)
		{
			T[i][j] = min(T[i - 1][j], T[i][j - 1]) + grid[i - 1][j - 1];
		}
	}

	return T[m][n];
}

// 120
int minimumTotal(vector<vector<int>>& triangle)
{
	vector<vector<int>> T(triangle);

	for (int i = 1; i < triangle.size(); ++i)
	{
		for (int j = 0; j < triangle[i].size(); ++j)
		{
			if (j == 0)
			{
				T[i][j] = T[i - 1][j] + triangle[i][j];
			}
			else if (j >= i)
			{
				T[i][j] = T[i - 1][j - 1] + triangle[i][j];
			}
			else
			{
				T[i][j] = min(T[i - 1][j - 1], T[i - 1][j]) + triangle[i][j];
			}
		}
	}

	int min = INT32_MAX;
	for (int c : T[triangle.size() - 1])
	{
		if (c < min)
			min = c;
	}

	return min;
}

// 264
int nthUglyNumber(int n)
{
	if (n == 0) return 0;
	if (n == 1) return 1;

	vector<int> T(n, 0);
	int p2 = 0, p3 = 0, p5 = 0;
	T[0] = 1;

	for (int i = 1; i < n; ++i)
	{
		//cout << "n: " << i << endl;
		//cout << "p2: " << p2 << "\tp3: " << p3 << "\tp5: " << p5 << endl;

		T[i] = min(T[p2] * 2, min(T[p3] * 3, T[p5] * 5));
		//cout << "T[i]: " << T[i] << '\n' << endl;
		if (T[i] == T[p2] * 2) ++p2;
		if (T[i] == T[p3] * 3) ++p3;
		if (T[i] == T[p5] * 5) ++p5;
	}

	return T[n - 1];
}

// 1326 **
// O(n^3) 优化到 O(n^2)
int minTaps(int n, vector<int>& ranges)
{
	//vector<vector<int>> T(n + 1, vector<int>(n + 1, 0));

	//for (int right = 0; right <= n; ++right)
	//{
	//	for (int left = right; left >= 0; --left)
	//	{
	//		for (int i = left; i <= right; ++i)
	//		{

	//		}
	//	}
	//}

	vector<int> T(n + 1, n + 2);
	T[0] = 0;

	for (int i = 0; i <= n; ++i)
	{
		for (int j = max(0, i - ranges[i]); j <= min(n, i + ranges[i]); ++j)
		{
			T[j] = min(T[j], T[max(0, i - ranges[i])] + 1);
		}
	}

	return T[n] == n + 2 ? -1 : T[n];
}

// 174 ***
int calculateMinimumHP(vector<vector<int>>& dungeon)
{
#pragma region Wrong Answer
	// XXX: Current answer is wrong. Min is lower, but Cur is higher.
// 	vector<vector<int>> b
//{
//	{1, -3, 3},
//	{ 0, -2, 0 },
//	{ -3, -3, -3 },
//};
	//int n = dungeon.size(), m = dungeon[0].size();
	//vector<vector<int>> T1(n + 1, vector<int>(m + 1, 0));	// minimun HP needed
	//vector<vector<int>> T2(n + 1, vector<int>(m + 1, 0));	// current HP

	//T1[1][1] = T2[1][1] = dungeon[0][0];
	//for (int i = 1; i <= n; ++i)
	//	T1[i][0] = T2[i][0] = INT_MIN >> 1;
	//for (int i = 1; i <= m; ++i)
	//	T1[0][i] = T2[0][i] = INT_MIN >> 1;

	//for (int i = 1; i <= n; ++i)
	//{
	//	for (int j = 1; j <= m; ++j)
	//	{
	//		if (i == 1 && j == 1) continue;

	//		int up_cur = T2[i - 1][j] + dungeon[i - 1][j - 1];
	//		int up_min = min(up_cur, T1[i - 1][j]);

	//		int left_cur = T2[i][j - 1] + dungeon[i - 1][j - 1];
	//		int left_min = min(left_cur, T1[i][j - 1]);

	//		if (left_min > up_min)
	//		{
	//			T1[i][j] = left_min;
	//			T2[i][j] = left_cur;
	//		}
	//		else if (left_min < up_min)
	//		{
	//			T1[i][j] = up_min;
	//			T2[i][j] = up_cur;
	//		}
	//		else
	//		{
	//			T1[i][j] = up_min;
	//			T2[i][j] = max(up_cur, left_cur);
	//		}
	//	}
	//}

	//return T1[n][m] > 0 ? 1 : -T1[n][m] + 1;

#pragma endregion

	int n = dungeon.size(), m = dungeon[0].size();
	vector<vector<int>> T(n + 1, vector<int>(m + 1, INT_MAX));

	T[n][m - 1] = T[n - 1][m] = 1;

	for (int i = n - 1; i >= 0; i--) 
	{
		for (int j = m - 1; j >= 0; j--) 
		{
			int need = min(T[i + 1][j], T[i][j + 1]) - dungeon[i][j];
			T[i][j] = need <= 0 ? 1 : need;

			printDPTable(T);
		}
	}

	return T[0][0];
}

// 198  *
int rob(vector<int>& nums)
{
	int len = nums.size();
	if (len == 0) return 0;
	if (len == 1) return nums[0];

	vector<int> T(len, 0);
	T[0] = nums[0];
	T[1] = nums[1] > nums[0] ? nums[1] : nums[0];

	for (int i = 2; i < len; ++i)
	{
		T[i] = max(T[i - 1], T[i - 2] + nums[i]);
	}

	return T[len - 1];
}

// 740 **
int deleteAndEarn(vector<int>& nums)
{
	vector<int> T(10001, 0);
	for (auto n : nums)
	{
		T[n] += n;
	}

	return rob(T);
}

// 221 **
int maximalSquare(vector<vector<int>>& matrix)
{
	if (matrix.empty()) return 0;
	int m = matrix.size(), n = matrix[0].size();
	vector<vector<int>> T(m, vector<int>(n, 0));

	int size = 0;
	for (int i = 0; i < m; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (i == 0 || j == 0 || matrix[i][j] == 0)
				T[i][j] = matrix[i][j] - 0;
			else
				T[i][j] = min(T[i - 1][j - 1], min(T[i - 1][j], T[i][j - 1])) + 1;

			size = max(size, T[i][j]);
		}
	}

	return size * size;
}

// 5 ***
string longestPalindrome(string s)
{
	int len = s.size();
	if (len <= 1) return s;

	vector<vector<bool>> T(len, vector<bool>(len, false));
	for (int i = 0; i < len; ++i)
		T[i][i] = true;

	int start = 0, maxDis = 0;
	for (int i = len; i >= 0; --i)
	{
		for (int dis = 1; dis < len - i; ++dis)
		{
			if (dis == 1 && s[i] == s[i + dis])
			{
				T[i][i + dis] = true;
			}
			else if (s[i] == s[i + dis])
			{
				T[i][i + dis] = T[i + 1][i + dis - 1];
			}

			if (T[i][i + dis] && dis > maxDis)
			{
				start = i;
				maxDis = dis;
			}

			cout << i << "   " << i + dis << endl;
			printDPTable(T);
		}
	}

	cout << "The last: " << endl;
	printDPTable(T);

	return s.substr(start, maxDis + 1);
}

string longestPalindrome2(string s)
{
	int len = s.size();
	if (len <= 1) return s;

	vector<vector<bool>> T(len, vector<bool>(len, false));
	for (int i = 0; i < len; ++i)
		T[i][i] = true;

	int start = 0, maxLen = 0;
	for (int l = 1; l < len; ++l)
	{
		for (int i = 0; i < len - l; ++i)
		{
			if (l == 1)
				T[i][i + l] = s[i] == s[i + l];
			else
				T[i][i + l] = (s[i] == s[i + l] && T[i + 1][i + l - 1]);

			if (T[i][i + l])
			{
				start = i;
				maxLen = l;
				cout << s.substr(start, maxLen + 1) << endl;
			}

			cout << i << "   " << i + l << endl;
			printDPTable(T);
		}
	}

	printDPTable(T);
	return s.substr(start, maxLen + 1);
}

// 139
bool wordBreak(string s, vector<string>& wordDict)
{
	if (wordDict.empty()) return false;

	vector<bool> T(s.size() + 1, false);
	T[0] = true;

	for (int i = 1; i <= s.size(); ++i)
	{
		for (int j = i - 1; j >= 0; --j)
		{
			if (T[j] && find(wordDict.cbegin(), wordDict.cend(), s.substr(j, i - j)) != wordDict.end())
			{
				T[i] = true;
				break;
			}
		}
	}

	return T[s.size()];
}

// 1553
class NOrangesSolution {
public:
	unordered_map<int, int> dp;
	int minDays(int n) {
		if (n <= 1)
			return n;

		if (dp.count(n) == 0)
			dp[n] = 1 + min(n % 2 + minDays(n / 2), n % 3 + minDays(n / 3));

		return dp[n];
	}
};

#pragma endregion

#pragma region Rely on O(n) Sub Problems

// 322
int coinChange(vector<int>& coins, int amount)
{
	vector<int> T(amount + 1, amount + 1);
	T[0] = 0;
	for (auto c : coins)
	{
		for (int i = c; i <= amount; ++i)
		{
			T[i] = min(T[i], T[i - c] + 1);
		}
	}

	return T[amount] > amount ? -1 : T[amount];
}

// 377 **
// 数组中不必连续地选取n个数，使和为目标值。每个数都可以重复使用
int combinationSum4(vector<int>& nums, int target)
{
	int n = nums.size();
	vector<int> T(target + 1, 0);
	T[0] = 1;

	for (int i = 1; i <= target; ++i)
	{
		for (auto n : nums)
		{
			if (n <= i)
				T[i] += T[i - n];
		}

		printDPTable(T);
	}

	
	// 这已经默认有序排列了：
	// 当计算目标值为3，2与1组合一定是先用1在用2.
	// 因为遍历nums是有序的。
	//for (auto num : nums)
	//{
	//	for (int i = num; i <= target; ++i)
	//		T[i] += T[i - num];

	//	cout << num << "  ->  " << target << endl;
	//	printDPTable(T);
	//}

	return T.back();
}

// 494 **
// 数组中不必连续地选取n个数，使和为目标值。每个数只能使用一次
// 输出选取方法总和数
int findTargetSumWays(vector<int>& nums, int S)
{
	int sum = accumulate(nums.begin(), nums.end(), 0);
	if (S > sum || S < -sum || ((sum + S) & 1) == 1 ) return 0;
	int target = (sum + S) / 2;

	vector<int> T(target + 1, 0);
	T[0] = 1;

	for (auto num : nums)
	{
		for (int i = target; i >= num; --i)
		//for (int i = num; i <= target; ++i)  // wrong
		{
			T[i] += T[i - num];
		}

		cout << target << "  ->  " << num << endl;
		printDPTable(T);
	}

	// Wrong
	// 这样默认数字可以交换位置，但题目表述来看是有序的
	//for (int i = 1; i <= target; ++i)
	//{
	//	for (int num : nums)
	//	{
	//		if (num <= i)
	//			T[i] += T[i - num];

	//		//cout << target << "  ->  " << num << endl;
	//		printDPTable(T);
	//	}
	//}


	return T.back();
}

// 416 **
// 数组中不必连续地选取n个数，使和为目标值。每个数只能使用一次
// 输出是否有选取方法
bool canPartition(vector<int>& nums)
{
	int sum = accumulate(nums.begin(), nums.end(), 0);
	if ((sum & 1) == 1) return false;
	sum >>= 1;

	//int n = nums.size();
	//vector<vector<bool>> T(n + 1, vector<bool>(sum + 1, false));
	//for (int i = 0; i <= n; ++i)
	//	T[i][0] = true;
	//for (int i = 1; i <= n; ++i)
	//{
	//	for (int j = 1; j <= sum; ++j)
	//	{
	//		T[i][j] = T[i - 1][j];
	//		if (j >= nums[i - 1])
	//			T[i][j] = T[i][j] || T[i - 1][j - nums[i - 1]];
	//	}
	//}
	//return T[n][sum];

	// Space: O(n*m) -> O(m)
	vector<bool> T(sum + 1, false);
	T[0] = true;

	for (int num : nums)
	{
		for (int i = sum; i >= num; --i)
		{
			T[i] = T[i] || T[i - num];
		}

		cout << sum << "  ->  " << num << endl;
		printDPTable(T);
	}

	return T.back();
}

// 279 *
int numSquares(int n)
{
	if (n <= 1) return n;

	vector<int> T(n + 1, INT_MAX >> 1);
	T[0] = 0;

	for (int i = 2; i <= n; ++i)
		for (int j = 1; j * j <= i; ++j)
			T[i] = min(T[i], T[i - j * j] + 1);

	return T[n];
}

// 312 ***
// 戳气球
int maxCoins(vector<int>& nums)
{
	int n = nums.size();
	nums.insert(nums.begin(), 1);
	nums.push_back(1);

	vector<vector<int>> T(n + 2, vector<int>(n + 2, 0));

	for (int r = 1; r <= n; ++r)
	{
		for (int l = r; l > 0; --l)
		{
			for (int i = l; i <= r; ++i)
			{
				//cout << "left: " << l << "    right: " << r << "    i: " << i << "    res: " << T[l][i - 1] + nums[l - 1] * nums[i] * nums[r + 1] + T[i + 1][r] << endl;
				T[l][r] = max(T[l][r], T[l][i - 1] + nums[l - 1] * nums[i] * nums[r + 1] + T[i + 1][r]);
			}

			//cout << "left: " << l << "    right: " << r << "    max: " << T[l][r] << '\n' << endl;
		}
	}

	return T[1][n];
}

// 300 ***
int lengthOfLIS(vector<int>& nums)
{
	if (nums.empty()) return 0;
	vector<int> T(nums.size(), 1);

	for (int i = 1; i < nums.size(); ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			if (nums[j] < nums[i])
				T[i] = max(T[j] + 1, T[i]);
		}
	}

	return *max_element(T.begin(), T.end());
}

// 85 ***
int maximalRectangle(vector<vector<char> > &matrix) {
	if (matrix.empty()) return 0;
	int n = matrix.size(), m = matrix[0].size(), res = 0;
	vector<vector<int>> T(n, vector<int>(m, 0));

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			if (matrix[i][j] == '1')
				T[i][j] = j == 0 ? 1 : (T[i][j - 1] + 1);
		}
	}

	printDPTable(T);

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < m; ++j)
		{
			int len = INT_MAX;

			for (int k = i; k < n; ++k)
			{
				len = min(len, T[k][j]);
				res = max(res, len * (k - i + 1));
			}
		}
	}

	return res;
}

#pragma endregion

// TODO
// 898
int subarrayBitwiseORs(vector<int>& arr) 
{
	return 0;
}

// 39
vector<vector<int>> combinationSum(vector<int>& candidates, int target)
{
	return vector<vector<int>>();
}

#pragma endregion

#pragma region Double

// 72
int MinDistance(string word1, string word2)
{
	int m = word1.size(), n = word2.size();
	vector<vector<int>> T(m + 1, vector<int>(n + 1, 0));
	for (int i = 1; i <= m; ++i) T[i][0] = i;
	for (int i = 1; i <= n; ++i) T[0][i] = i;

	for (int i = 1; i <= m; ++i)
	{
		for (int j = 1; j <= n; ++j)
		{
			if (word1[i - 1] == word2[j - 1])
				T[i][j] = T[i - 1][j - 1];
			else
				T[i][j] = min(min(T[i - 1][j - 1], T[i - 1][j]), T[i][j - 1]) + 1;
		}
	}

	return T[m][n];
}

// 1143
int longestCommonSubsequence(string text1, string text2)
{
	int n = text1.size(), m = text2.size();
	vector<vector<int>> T(n + 1, vector<int>(m + 1, 0));

	for (int i = 1; i <= n; ++i)
	{
		for (int j = 1; j <= m; ++j)
		{
			if (text1[i - 1] == text2[j - 1])
				T[i][j] = T[i - 1][j - 1] + 1;
			else
				T[i][j] = max(T[i - 1][j], T[i][j - 1]);
		}
	}

	return T[n][m];
}

// 718
int findLength(vector<int>& A, vector<int>& B)
{
	vector<vector<int>> T(A.size() + 1, vector<int>(B.size() + 1, 0));
	int res = 0;

	for (int i = 1; i <= A.size(); ++i)
	{
		for (int j = 1; j <= B.size(); ++j)
		{
			res = max(res, T[i][j] = A[i - 1] == B[j - 1] ? T[i - 1][j - 1] + 1 : 0);
		}
	}

	printDPTable(T);

	return res;
}

// TODO
// 10
bool isMatch(string s, string p)
{
	return false;
}

#pragma endregion

#pragma endregion

#pragma region Array and List and Hash Table

// 2
ListNode* addTwoNumbers(ListNode* l1, ListNode* l2)
{
	ListNode target(0), *p = &target;
	int extra = 0;

	while (l1 || l2 || extra)
	{
		if (l1) { extra += l1->val; l1 = l1->next; }
		if (l2) { extra += l2->val; l2 = l2->next; }

		p->next = new ListNode(extra % 10);
		p = p->next;

		extra /= 10;
	}

	return target.next;
}

// 4
double findMedianSortedArrays(vector<int>& nums1, vector<int>& nums2)
{
	int const len1 = nums1.size(), len2 = nums2.size();
	int n1 = 0, n2 = 0, n = 0;
	vector<int> newNums(len1 + len2, 0);

	while (n1 < len1 || n2 < len2)
	{
		if (n1 >= len1) { newNums[n++] = (nums2[n2++]); continue; }
		if (n2 >= len2) { newNums[n++] = (nums1[n1++]); continue; }

		if (nums2[n2] >= nums1[n1])
			newNums[n++] = (nums1[n1++]);
		else
			newNums[n++] = (nums2[n2++]);
	}

	if ((len1 + len2) % 2 == 1)
		return newNums[(len1 + len2) / 2];
	else
		return (double)(newNums[(len1 + len2) / 2 - 1] + newNums[(len1 + len2) / 2]) / 2;
}

// 141
bool hasCycle(ListNode *head) 
{
	ListNode *walker, *runner;
	walker = runner = head;

	while (walker != nullptr && runner != nullptr)
	{
		walker = walker->next;
		runner = runner->next;
		if (runner == nullptr) return false;
		runner = runner->next;

		if (walker == runner) return true;
	}

	return false;
}

// 142
ListNode *detectCycle(ListNode *head) 
{
	if (head == nullptr) return nullptr;

	ListNode *walker, *runner;
	walker = runner = head;
	bool hasCycle = false;

	while (runner->next != nullptr && runner->next->next != nullptr)
	{
		walker = walker->next;
		runner = runner->next->next;
		if (walker == runner)
		{
			hasCycle = true;
			break;
		}
	}

	if (hasCycle == false) return nullptr;

	walker = head;
	while (walker != runner)
	{
		walker = walker->next;
		runner = runner->next;
	}

	return walker;
}

// 206
ListNode* reverseList(ListNode* head) 
{
	ListNode* prev = nullptr;
	ListNode* curr = head;
	ListNode* next;

	while (curr) {
		next = curr->next;
		curr->next = prev;
		prev = curr;
		curr = next;
	}

	return prev;
}

// 560 **
int subarraySum(vector<int>& nums, int k)
{
	// DP
	//int n = nums.size();
	//vector<int> sums(n, 0);
	//sums[0] = nums[0];
	//for (int i = 1; i < n; ++i)
	//	sums[i] = sums[i - 1] + nums[i];

	//int length = 0;
	//for (int i = 0; i < n; ++i)
	//{
	//	for (int j = i; j >= 0; --j)
	//	{
	//		
	//		if (sums[i] - sums[j] == k)
	//			length = min(length, i - j + 1);
	//	}
	//}
	//return length;

	int sum = 0, result = 0;
	unordered_map<int, int> preSum;
	preSum[0] = 1;

	for (int i = 0; i < nums.size(); ++i)
	{
		sum += nums[i];
		auto it = preSum.find(sum - k);
		if (it != preSum.end())
			result += it->second;

		it = preSum.find(sum);
		preSum[sum] = (it == preSum.end() ? 0 : it->second) + 1;
	}

	return result;
}

#pragma endregion

#pragma region LinkedList

// 23
ListNode* mergeKLists(vector<ListNode*>& lists)
{
	if (lists.empty()) return nullptr;
	if (lists.size() == 1) return lists[0];
}
ListNode* mergeTwoLists(ListNode* l1, ListNode* l2)
{
	return nullptr;
}

#pragma endregion

#pragma region Two Pointers

// 3
int lengthOfLongestSubstring(string s)
{
	map<char, int> map;
	int start = -1, maxLen = 0;

	for (int i = 0; i < s.size(); ++i)
	{
		if (map.count(s[i]) != 0)
			start = max(start, map[s[i]]);

		map[s[i]] = i;
		maxLen = max(maxLen, i - start);
	}

	return maxLen;
}

// 11
int maxArea(vector<int>& height)
{
	int l = 0, r = height.size() - 1;
	int res = 0;

	while (l < r)
	{
		res = max(res, min(height[l], height[r]) * (r - l));

		if (height[l] > height[r])
			--r;
		else
			++l;
	}

	return res;
}

// 42
int trap(vector<int>& height)
{
	if (height.size() < 3) return 0;

	auto l = height.begin(), r = height.end() - 1;
	int lower = 0, current = 0, result = 0;

	while (l != r + 1)
	{
		current = *l > *r ? *r-- : *l++;
		lower = max(lower, current);
		result += lower - current;
	}

	return result;
}

#pragma endregion

#pragma region Tree

// 94
//vector<int> inorderTraversal(TreeNode* root)
//{
//	if (root)
//	{
//		inorderTraversal(root->left);
//		cout << root->val << endl;
//		inorderTraversal(root->right);
//	}
//}

// 96
int numTrees(int n)
{
	vector<int> T(n + 1, 0);
	T[0] = T[1] = 1;

	for (int i = 2; i <= n; ++i)
	{
		for (int j = 1; j <= i; ++j)
		{
			T[i] += T[j - 1] * T[i - j];
		}
	}

	return T[n];
}

// 102
vector<vector<int>> levelOrder(TreeNode* root)
{
	vector<vector<int>> res(0);
	queue<TreeNode*> q;
	q.push(root);

	while (!q.empty())
	{
		int n = q.size();
		vector<int> curLevel;
		curLevel.reserve(n);

		for (int i = 0; i < n; ++i)
		{
			TreeNode* p = q.front();
			q.pop();

			if (p != nullptr)
			{
				curLevel.push_back(p->val);
				q.push(p->left);
				q.push(p->right);
			}
		}

		if (curLevel.size() > 0)
			res.push_back(move(curLevel));
	}

	return res;
}

#pragma endregion

#pragma region DFS and Backtracking

// 22
class GenerateParenthesesSolution {
public:
	vector<string> generateParenthesis(int n) 
	{
		vector<string> res;
		helper(res, "", n, n);
		return res;
	}

	void helper(vector<string>& res, string&& s, int n, int m)
	{
		if (n == 0 && m == 0)
		{
			res.push_back(s);
			return;
		}

		if (n > 0)
			helper(res, s + '(', n - 1, m);

		if (m > n)
			helper(res, s + ')', n, m - 1);
	}
};

// TODO
// 37
class SudokuSolverSolution {
public:
	void solveSudoku(vector<vector<char>>& board) 
	{
		_row = vector<vector<bool>>(9, vector<bool>(10, false));
		_col = vector<vector<bool>>(9, vector<bool>(10, false));
		_box = vector<vector<bool>>(9, vector<bool>(10, false));

		for (int i = 0; i < 9; ++i)
		{
			for (int j = 0; j < 9; ++j)
			{
				if (board[i][j] != '.')
				{
					_col[i][board[i][j] - '0'] = true;
					_row[j][board[i][j] - '0'] = true;
					_box[i / 3 * 3 + j / 3][board[i][j] - '0'] = true;
				}
			}
		}

		fill(board, 0, 0);
	}

	void fill(vector<vector<char>>& board, int n, int m)
	{

	}

private:
	vector<vector<bool>> _row, _col, _box;
};

// 46
class PermutationsSolution {
public:
	vector<vector<int>> permute(vector<int>& nums) 
	{
		for (auto i : nums)
			map[i] = false;

		fill(nums);

		return res;
	}

	void fill(vector<int>& nums)
	{
		if (temp.size() == nums.size())
			res.push_back(temp);

		for (auto i : nums)
		{
			if (map[i]) continue;

			map[i] = true;
			temp.push_back(i);
			fill(nums);
			map[i] = false;
			temp.pop_back();
		}
	}

private:
	vector<vector<int>> res;
	vector<int> temp;
	unordered_map<int, bool> map;
};

// TODO
// 51
class NQueensSolution {
public:
	vector<vector<string>> solveNQueens(int n) {

	}
};

// 79
class WordSearchSolution {
public:
	bool exist(vector<vector<char>>& board, string word) 
	{
		if (board.empty()) return false;

		n = board.size(), m = board[0].size();

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				if (find(i, j, 0, board, word)) return true;
			}
		}

		return false;
	}

	bool find(int x, int y, int index, vector<vector<char>>& board, string& word)
	{
		if (x < 0 || x >= n || y < 0 || y >= m || board[x][y] != word[index])
			return false;

		if (index == word.length() - 1) return true;

		char cur = word[index];
		board[x][y] = 0;

		bool b = find(x - 1, y, index + 1, board, word)
			|| find(x + 1, y, index + 1, board, word)
			|| find(x, y - 1, index + 1, board, word)
			|| find(x, y + 1, index + 1, board, word);

		board[x][y] = cur;

		return b;
	}

private:
	int n, m;
};

// 200
class IslandsNumberSolution {
public:
	int numIslands(vector<vector<char>>& grid) 
	{

	}

	void find(vector<vector<char>>& grid, int i, int j)
	{

	}
};

// 526
class BeautifulArrangementSolution {
public:
	int countArrangement(int n) 
	{
		for (int i = 1; i <= n; ++i)
			map[i] = false;

		fill(n);

		return res;
	}

	void fill(int n)
	{
		if (temp.size() == n)
		{
			++res;
			return;
		}

		for (int i = 1; i <= n; ++i)
		{
			if (map[i])
				continue;

			map[i] = true;
			temp.push_back(i);
			if (check()) fill(n);
			map[i] = false;
			temp.pop_back();
		}
	}

	bool check()
	{
		for (int i = 1; i <= temp.size(); ++i)
		{
			if (i % temp[i - 1] != 0 && temp[i - 1] % i != 0)
				return false;
		}

		return true;
	}

private:
	unordered_map<int, bool> map;
	vector<int> temp;
	int res = 0;
};

// 698
class PartitionKSubsetsSolution {
public:
	// use global variables to avoid long parameter list
	int target; // of each bucket
	vector< int > ns;
	vector< int > bucket;

	//bool canPartitionKSubsets(vector<int>& nums, int k) {
	//	int sum = 0;
	//	for (int &n : nums) sum += n;
	//	if (sum % k) return false; // not divisible
	//	target = sum / k;
	//	ns = vector< int >(nums);
	//	bucket = vector< int >(k, 0);
	//	// starting with bigger ones makes it faster
	//	sort(ns.begin(), ns.end(), greater<int>());
	//	//reverse(ns.begin(), ns.end());
	//	return put(0);
	//}

	//// put #n item of ns into some bucket to meet target
	//bool put(int n) {
	//	for (int i = 0; i < bucket.size(); ++i) {
	//		if (bucket[i] + ns[n] > target) continue; // try next bucket
	//		bucket[i] += ns[n]; // put it in!
	//		if (n == ns.size() - 1) return true; // all items in bucket, no overflow
	//		if (put(n + 1)) return true; // move on to next item
	//		else { // no solution = wrong bucket
	//			bucket[i] -= ns[n]; // take it out
	//			if (bucket[i] == 0) return false; // no need to try other empty bucket
	//		}
	//	}
	//	return false; // no bucket fits
	//}

	bool canPartitionKSubsets(vector<int>& nums, int k) {
		const int sum = accumulate(nums.begin(), nums.end(), 0);
		if (sum%k != 0) return false;
		sort(nums.rbegin(), nums.rend());
		return dfs(nums, sum / k, 0, k, 0);
	}

	bool dfs(const vector<int>&nums, int target, int cur, int k, int used) {
		if (k == 0) return used == (1 << nums.size()) - 1;
		for (int i = 0; i < nums.size(); i++) {
			if (used&(1 << i)) continue;
			int t = cur + nums[i];
			if (t > target) break;
			int new_used = used | (1 << i);
			if (t == target && dfs(nums, target, 0, k - 1, new_used)) return true;
			else if (dfs(nums, target, t, k, new_used)) return true;
		}
		return false;
	}
};

#pragma endregion

#pragma region Math

// 7
int reverse(int x)
{
	int prevRec = 0, rev = 0;

	while (x != 0)
	{
		rev = rev * 10 + x % 10;
		if ((rev - x % 10) / 10 != prevRec)
			return 0;
		prevRec = rev;
		x /= 10;
	}

	return rev;
}

// 48
void rotate(vector<vector<int>>& matrix)
{
	reverse(matrix.begin(), matrix.end());

	for (int i = 0; i < matrix.size(); ++i)
	{
		for (int j = i + 1; j < matrix[i].size(); ++j)
		{
			swap(matrix[i][j], matrix[j][i]);
		}
	}
}

#pragma endregion

vector<int> CurrentProblem(vector<int>& nums)
{
	priority_queue<int> q;
	for (auto num : nums)
		q.push(num);

	vector<int> res;

	for (int i = 1; i <= 100; ++i)
	{
		res.push_back(q.top());
		q.pop();
	}

	return res;
}

void main()
{
	vector<int> a = { 1,1,1 };
	vector<int> a1 = { 10,20,30 };
	vector<int> a2 = { 1, 100, 1, 1, 1, 100, 1, 1, 100, 1 };
	//vector<vector<char>> T
	//{
	//	{ '0','0','0','0','0','0','1' },
	//	{ '0','0','0','0','1','1','1'},
	//	{ '1','1','1','1','1','1','1' },
	//	{ '0','0','0','1','1','1','1' }
	//};
	vector<vector<char>> T
	{
		{ '1','0','1','0','0' },
		{ '1','0','1','1','1'},
		{ '1','1','1','1','1' },
		{ '1','0','0','1','0' }
	};
	vector<vector<int>> b
	{ 
		{1, -3, 3},
		{0, -2, 0},
		{-3, -3, -3},
	};
	vector<string> c{ "Leet", "Code" };
	string s1 = "asdsadag";
	string s2 = "b";

	TreeNode *root = &TreeNode(1);
	root->left = &TreeNode(2);
	root->right = &TreeNode(3);
	root->left->left = &TreeNode(4);
	root->left->right = &TreeNode(5);
	root->right->left = &TreeNode(6);
	root->right->right = &TreeNode(7);

	ListNode *node = &ListNode(1);
	node->next = &ListNode(2);
	node->next->next = &ListNode(3);
	node->next->next->next = &ListNode(4);
	node->next->next->next->next = &ListNode(5);
	node->next->next->next->next->next = node->next->next;

	cout << findTargetSumWays(a, 1) << endl;

	// 网易一面 4.6 60min
	//1. map和unordered_map；哪种占用内存更大；unordered_map扩容与删除。
	//2. vector扩容；push_back平均时间复杂度。
	//3. 有序的1到1000寻找缺失的一个数；无序情况怎么做。
	//4. TopK；内存不够的情况；有满分的情况。
	//5. 说说A*；cost函数怎么定义；NavMesh怎么寻路；其他寻路算法。

	// 紫龙一面 4.6 30min
	//1. vector；怎么插入；怎么扩容。
	//2. map与unordered_map；什么时候用哪个。
	//3. static关键字；修饰函数中一临时变量时，多久生成、多久销毁。
	//4. 四个cast；dynamic_cast转换失败返回什么。
	//5. c#值类型与引用类型；装箱与拆箱。
	//6. c# GC。
	//7. 二叉树遍历；递归怎么改成循环。
	//8. 向量点乘与叉乘；举例使用时机。
	//9. 变换矩阵为什么是4x4；一般是（0, 0, 0, 1）的最后一列\行什么作用。
	//10. 一个顶点经历哪些空间与变换才会到屏幕上。
	//11. BlinnPhong光照模型；n点乘l和n点乘h意义。
	//12. shadowmap过程；shadowmap问题与解决方案；照向四面八方的点光源怎么生成shadowmap。
	//13. 骨骼与蒙皮；蒙皮算法。
	//14. 寻路算法；它们异同；说说A*。

	// 紫龙二面 4.14 80min
	// 1. i++线程安全吗。（local static线程安全，shared_ptr线程安全）
	// 2. const修饰函数什么作用。
	// 3. 字节对齐；数据传输时对齐到8个字节的数据怎么压缩回5个字节。
	// 4. float占用几个字节；64位中呢；float类型数据在while循环中一直加一，会溢出吗。
	// 5. 虚继承使用意义是什么。
	// 6. static修饰全局变量有什么意义。
	// 7. 什么情况使用map而不是unordered_map；hash函数定义；hash后的均匀分布。
	// 8. STL中sort是稳定的吗；什么样的排序算法称之为稳定的。
	// 9. dynamic_cast什么时候抛出异常。
	// 10. 树的广搜。
	// 11. Dijkstra寻路算法；A*什么时候退化为Dijkstra。
	// （Dijkstra属于Greedy算法，负权边要用DP的Bellman Ford算法）
	// 12. Mipmap作用；详细的纹理映射流程，软光栅与硬件的不同。
	// 13. 由四个点围成一个矩形，如何判断第五个是否在矩形内。
	// 14. 欧拉角和Quaternion；欧拉角旋转。
	// 15. 变换矩阵为什么是4x4；在一个矩阵同时表述旋转、缩放和平移时，三个动作的先后顺序。
	// 紧接就HR面，差不多半小时。

	// “高塔不立于浮沙之上”，“知识就像树一样，有离根很近的，有离根很远的”

	// 腾讯一面 4.8 70min
	//1. 3C部分（3C是Character, Control和Camera）。扩展聊了聊IK，更真实的动作系统，基于物理的动作，相机策略。
	//2. 项目中的跑图算法。基础是Flocking群聚算法，项目中引入了视野和道路感知解决某些实际问题。
	//3. 游戏软体。
	//4. 运行时资源管理。
	//5. ECS和帧同步。
	//6. GC。
	//7. 图形学和渲染管线。
	//8. 碰撞检测。
	//9. 骨骼与蒙皮。

	// 腾讯二面 4.15 30min
	// 1. 更完善的动画。
	// 2. ECS；为什么会有误差（float）；Data Oriented和Object Oriented在项目上有什么区别。
	// 3. 跑图算法细节，随机性怎么产生的。

	// 腾讯三面 4.21 30min
	// 1. 相机运动策略；引入弹簧系统怎么运作；弹簧简谐运动怎么停下来。
	// 2. 质点弹簧系统；timestep大小有要求吗；timestep为什么是恒定的；基于力和基于位置哪种更适合游戏。
	// 3. AABB和多面体与AABB和三角形哪种碰撞检测更复杂；描述下AABB与三角形碰撞。
	// 4. 物理引擎；碰撞检测流程；BVH和KDTree优劣；KDTree和Octree比较。
	// 5. 物理引擎怎么优化；多线程相关；ECS将相同类型数据放在连续内存中，为什么对性能有提升。
	// 6. SIMD；向量运算优化；
	// 7. 反射了解吗；反射内存偏移怎么计算；Unreal怎么做的反射。
	// 8. std::move了解吗；move一个const变量会发生什么；将其传给构造函数，会调用const构造还是非const移动构造。

	// 莉莉丝一二面 4.9 60min
	// 1. 对于空间中多个点，最小AABB怎么生成；最小包围球怎么生成。
	// 2. 欧拉角和Quaternion；对于一个正向物体，欧拉角怎么偏转到(x,y,z)固定角度；欧拉角怎么转化到矩阵。
	// 3. string s; vector<string> v; v.emplace_back(std::move(s));
	// 4. 跟随的人物转过墙角时，相机怎么运动。
	// 5. 纹理映射。
	// 6. 软体模拟。
	// 7. 物理引擎了解吗，碰撞相关。

	// 米哈游一面 4.13 50min
	//1. C++与C#区别。
	//2. C++的struct与class；C++4个cast；C#的ref与out。
	//3. 模板。
	//4. quaternion。
	//5. 渲染管线；深度测试。
	//6. 行为树与状态机。

	// 米哈游二面 4.15 60min
	// 1. 任务编辑器。如果策划想要在编辑完立刻验证逻辑，怎么处理单机与联网状态下的任务数据。
	// 2. 角色。设计一个良好的启、停动作逻辑。
	// 3. 跑图算法。这里面的避让和碰撞相关，如果场景很复杂你怎么重新设计算法中的碰撞检测和处理。
	// 4. 软体；物理引擎。（这部分是我发起的讨论，因为比较熟）
	// 5. 一个游戏需要长时间运行，怎么设计整个游戏系统与各个模块，避免内存碎片化。
	// 6. 设计Graphics Driver，光栅流程兼容OpenGL和DX。
	// 紧接着就是HR面，不到半小时。
}

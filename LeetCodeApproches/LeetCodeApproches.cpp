#include <iostream>
#include <string>
#include <vector>
#include <queue>
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

void printDPTable(vector<vector<int>>&& T)
{
	for (auto a : T)
	{
		for (auto b : a)
		{
			cout << b << "\t";
		}
		cout << '\n' << endl;
	}
	cout << '\n';
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

#pragma endregion

#pragma region Rely on O(n) Sub Problems

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

// 312 *
int maxCoins(vector<int>& nums)
{
	int n = nums.size();
	nums.insert(nums.begin(), 1);
	nums.push_back(1);

	vector<vector<int>> T(n + 2, vector<int>(n + 2, 0));

	for (int len = 1; len <= n; ++len)
	{
		for (int left = 1; left <= n - len + 1; ++left)
		{
			int right = left + len - 1;
			for (int k = left; k <= right; ++k)
			{
				T[left][right] = max(T[left][right], T[left][k - 1] + (nums[left - 1] * nums[k] * nums[right + 1]) + T[k + 1][right]);
			}
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

#pragma endregion

#pragma region Double

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

	return res;
}

#pragma endregion

#pragma endregion

#pragma region Array and List

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

//for (int i = 0; i < n; ++i)
//{
//	if (!available()) continue;	// 剪枝
//	put();	// 放置
//	fill(); // 进入
//	back(); // 回退
//}

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

#pragma region LinkedLists

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

int CurrentProblem(vector<int>& nums)
{
	int n = nums.size();
	if (n == 0) return 0;
	vector<int> T(n + 1, 1);

	for (int i = 1; i <= n; ++i)
	{
		for (int j = 1; j < i; ++j)
		{
			if (nums[j - 1] < nums[i - 1])
			{
				T[i] = max(T[i], T[j] + 1);
			}
		}
	}

	return *max_element(T.begin(), T.end());
}

void main()
{
	vector<int> a = { 0,1,0,2,1,0,1,3,2,1,2,1 };
	vector<int> a1 = { 7, 5, 6, 8, 1, 4, 2, 9, 3 };
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
	vector<vector<int>> b{ {1} };
	vector<string> c{ "Leet", "Code" };
	string s1 = "a";
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

	
	cout << CurrentProblem(a1) << endl;
}

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <numeric>
#include <algorithm>
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

void printDPTable(vector<vector<int>>& T)
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

void printDPTable(vector<vector<bool>>& T)
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

// 5
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
		}
	}

	printDPTable(T);

	return s.substr(start, maxDis + 1);
}

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

// 53
int maxSubArray(vector<int>& nums)
{
	int length = nums.size();
	vector<int> T(length, 0);
	int m = T[0];

	for (int i = 1; i < length; ++i)
	{
		T[i] = T[i - 1] > 0 ? T[i - 1] + nums[i] : nums[i];
		m = max(m, T[i]);
	}

	return m;
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

// 198
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

// 221
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
			if (T[j])
			{
				string word = s.substr(j, i - j);
				if (find(wordDict.begin(), wordDict.end(), word) != wordDict.end())
				{
					T[i] = true;
					break;
				}
			}
		}
	}

	return T[s.size()];
}

// 1143
int longestCommonSubsequence(string text1, string text2)
{
	vector<vector<int>> T(text1.size(), vector<int>(text2.size(), 0));
	T[0][0] = text1[0] == text2[0];

	for (int i = 0; i < text1.size(); ++i)
	{
		for (int j = 0; j < text2.size(); ++j)
		{
			if (i == 0 && j == 0) continue;

			if (i == 0)
				T[i][j] = T[i][j - 1] | text1[i] == text2[j];
			else if (j == 0)
				T[i][j] = T[i - 1][j] | text1[i] == text2[j];
			else if (text1[i] != text2[j])
				T[i][j] = max(T[i - 1][j], T[i][j - 1]);
			else
				T[i][j] = T[i - 1][j - 1] + 1;
		}
	}

	return T[text1.size() - 1][text2.size() - 1];
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

// 300
int lengthOfLIS(vector<int>& nums)
{
	if (nums.empty()) return 0;
	vector<int> T(nums.size(), 1);
	int res = 1;

	for (int i = 1; i < nums.size(); ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			if (nums[j] < nums[i])
				T[i] = max(T[j] + 1, T[i]);
		}
	}

	for (auto i : T)
		res = max(res, i);

	return res;
}

// 312
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

// 898
class Solution {
public:
	int subarrayBitwiseORs(vector<int>& arr) {

	}
};

// 698
class Solution {
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

#pragma region Array

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
	if (n == 1) return 1;
	if (n == 2) return 2;

	vector<int> T(n + 1, 0);
	T[0] = T[1] = 1;
	T[2] = 2;

	for (int i = 3; i <= n; ++i)
	{
		for (int j = 0; j < i; ++j)
		{
			T[i] += T[j] * T[i - j - 1];
		}
	}

	return T[n];
}

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

void main()
{
	//vector<int> a = { 1,2 ,3,4,5};
	//vector<int> a1 = { 3,4 };
	//vector<vector<int>> b{ {1} };
	//vector<string> c{ "Leet", "Code" };
}

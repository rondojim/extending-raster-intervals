#include "../../include/mbr_algos/axis_intersection.h"
#include <map>
#include <set>
#include <vector>

static void update(std::vector<std::vector<int>> &segt, int root, int lo,
                   int hi, int x, int y) {
  segt[root].push_back(y);
  if (lo == hi) {
    return;
  } else {
    int mid = (lo + hi) / 2;
    int left = 2 * root;
    int right = left + 1;

    if (x <= mid) {
      update(segt, left, lo, mid, x, y);
    } else {
      update(segt, right, mid + 1, hi, x, y);
    }
  }
}

// Function to query the segment tree for points within a given range
static bool query(std::vector<std::vector<int>> &segt, int root, int lo, int hi,
                  int x1, int x2, int y1, int y2) {
  if (x1 <= lo && hi <= x2) {
    return ((upper_bound(segt[root].begin(), segt[root].end(), y2) -
             lower_bound(segt[root].begin(), segt[root].end(), y1)));
  } else if (hi < x1 || lo > x2) {
    return false;
  } else {
    int mid = (lo + hi) / 2;
    int left = 2 * root;
    int right = left + 1;

    return (query(segt, left, lo, mid, x1, x2, y1, y2) ||
            query(segt, right, mid + 1, hi, x1, x2, y1, y2));
  }
}

long long calculate_axis_intersection(std::vector<mbr> &lhs,
                                      std::vector<mbr> &rhs, bool Xaxis,
                                      int MAXC) {
  std::vector<std::pair<int, int>> lsegm;
  std::vector<std::pair<int, int>> rsegm;

  for (auto &it : lhs) {
    if (Xaxis) {
      lsegm.push_back({it.minc.first, it.maxc.first});
    } else {
      lsegm.push_back({it.minc.second, it.maxc.second});
    }
  }

  for (auto &it : rhs) {
    if (Xaxis) {
      rsegm.push_back({it.minc.first, it.maxc.first});
    } else {
      rsegm.push_back({it.minc.second, it.maxc.second});
    }
  }

  std::vector<std::vector<int>> segt;
  for (int i = 0; i <= 4 * MAXC; ++i) {
    std::vector<int> vec;
    segt.push_back(vec);
  }

  for (auto &it : lsegm) {
    update(segt, 1, 1, MAXC, it.first, it.second);
  }

  long long result = 0;

  for (auto &it : rsegm) {
    int x = it.first;
    int y = it.second;

    result += query(segt, 1, 1, MAXC, 1, x, x, MAXC);
    result += query(segt, 1, 1, MAXC, x + 1, MAXC, 1, x);
  }

  return result;
}

#include "../../include/join_algos/node.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/progressbar.hpp"
#include <algorithm>
#include <future>
#include <map>
#include <set>
#include <vector>

struct mbr {
  std::pair<int, int> minc;
  std::pair<int, int> maxc;
  int id;
};

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

static bool query(std::vector<std::vector<int>> &segt, int root, int lo, int hi,
                  int x1, int x2, int y1, int y2) {
  if (x1 <= lo && hi <= x2) {
    return ((upper_bound(segt[root].begin(), segt[root].end(), y2) -
             lower_bound(segt[root].begin(), segt[root].end(), y1)) > 0);
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

void mbr_no_outside(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                    std::map<std::pair<int, int>, int> &compress, int MAXX,
                    int MAXY, std::set<int> &result) {
  std::vector<std::vector<int>> segt;
  for (int i = 0; i <= 4 * MAXX; ++i) {
    std::vector<int> vec;
    segt.push_back(vec);
  }
  std::vector<mbr> mlhs;
  std::vector<mbr> mrhs;

  for (auto &it : lhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];
    mbr cur;
    cur.minc = {min_x, min_y};
    cur.maxc = {max_x, max_y};
    cur.id = it.polygon_id;
    mlhs.push_back(cur);
  }

  for (auto &it : rhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];
    mbr cur;
    cur.minc = {min_x, min_y};
    cur.maxc = {max_x, max_y};
    cur.id = it.polygon_id;
    mrhs.push_back(cur);
  }

  int total = 0;

  for (auto &it : mlhs) {
    update(segt, 1, 1, MAXX, it.minc.first, it.minc.second);
    update(segt, 1, 1, MAXX, it.maxc.first, it.maxc.second);
    update(segt, 1, 1, MAXX, it.minc.first, it.maxc.second);
    update(segt, 1, 1, MAXX, it.maxc.first, it.minc.second);
  }

  for (auto &it : segt) {
    if (it.size() > 1) {
      std::sort(it.begin(), it.end());
    }
  }

  for (auto &it : mrhs) {
    int min_x = it.minc.first;
    int max_x = it.maxc.first;
    int min_y = it.minc.second;
    int max_y = it.maxc.second;

    if (query(segt, 1, 1, MAXX, min_x, max_x, min_y, max_y) > 0) {
      result.insert(it.id);
    }
  }

  for (int i = 0; i <= 4 * MAXX; ++i) {
    segt[i].clear();
  }

  for (auto &it : mrhs) {
    update(segt, 1, 1, MAXX, it.minc.first, it.minc.second);
    update(segt, 1, 1, MAXX, it.maxc.first, it.maxc.second);
    update(segt, 1, 1, MAXX, it.minc.first, it.maxc.second);
    update(segt, 1, 1, MAXX, it.maxc.first, it.minc.second);
  }

  for (auto &it : segt) {
    if (it.size() > 1) {
      std::sort(it.begin(), it.end());
    }
  }

  for (auto &it : mlhs) {
    int min_x = it.minc.first;
    int max_x = it.maxc.first;
    int min_y = it.minc.second;
    int max_y = it.maxc.second;

    if (query(segt, 1, 1, MAXX, min_x, max_x, min_y, max_y) > 0) {
      result.insert(it.id);
    }
  }
}
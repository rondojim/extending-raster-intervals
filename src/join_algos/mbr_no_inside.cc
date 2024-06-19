#include "../../include/join_algos/mbr_no_inside.h"
#include "../../include/join_algos/interval.h"
#include "../../include/utils/geometry_types.h"
#include <algorithm>
#include <map>
#include <set>
#include <vector>

struct mbr {
  std::pair<int, int> minc;
  std::pair<int, int> maxc;
  int id;
};

static void add(std::vector<int> &B, int x, int val) {
  for (; x < B.size(); x += x & -x) {
    B[x] += val;
  }
}

static void range_add(std::vector<int> &B1, std::vector<int> &B2, int l, int r,
                      int val) {
  add(B1, l, val);
  add(B1, r + 1, -val);
  add(B2, l, val * (l - 1));
  add(B2, r + 1, -val * r);
}

static int sum(std::vector<int> &B, int x) {
  int res = 0;
  for (; x > 0; x -= x & -x) {
    res += B[x];
  }
  return res;
}

static int prefix_sum(std::vector<int> &B1, std::vector<int> &B2, int x) {
  return sum(B1, x) * x - sum(B2, x);
}

static int range_sum(std::vector<int> &B1, std::vector<int> &B2, int l, int r) {
  return prefix_sum(B1, B2, r) - prefix_sum(B1, B2, l - 1);
}

void sweep_on_X(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                std::map<std::pair<int, int>, int> &compress, int MAXX,
                int MAXY, std::set<int> &result) {

  std::vector<Interval> intervals;

  for (auto &it : lhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];

    Interval left_interval = Interval(min_x, min_y, max_y, it.polygon_id, true);
    Interval right_interval =
        Interval(max_x, min_y, max_y, it.polygon_id, false);

    intervals.push_back(left_interval);
    intervals.push_back(right_interval);
  }

  for (auto &it : rhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];

    Interval left_interval = Interval(min_x, min_y, max_y, it.polygon_id, true);
    Interval right_interval =
        Interval(max_x, min_y, max_y, it.polygon_id, false);

    intervals.push_back(left_interval);
    intervals.push_back(right_interval);
  }

  std::vector<int> B1;
  std::vector<int> B2;

  B1.assign(MAXY + 1, 0);
  B2.assign(MAXY + 1, 0);

  std::sort(intervals.begin(), intervals.end(),
            IntervalComparePositiveIdPriority());

  for (auto &it : intervals) {
    if (it.id > 0) {
      if (range_sum(B1, B2, it.min_y, it.max_y) > 0) {
        result.insert(it.id);
      }
    } else {
      if (it.is_start) {
        range_add(B1, B2, it.min_y, it.max_y, 1);
      } else {
        range_add(B1, B2, it.min_y, it.max_y, -1);
      }
    }
  }

  B1.assign(MAXY + 1, 0);
  B2.assign(MAXY + 1, 0);

  std::sort(intervals.begin(), intervals.end(),
            IntervalCompareNegativeIdPriority());

  for (auto &it : intervals) {
    if (it.id < 0) {
      if (range_sum(B1, B2, it.min_y, it.max_y) > 0) {
        result.insert(it.id);
      }
    } else {
      if (it.is_start) {
        range_add(B1, B2, it.min_y, it.max_y, 1);
      } else {
        range_add(B1, B2, it.min_y, it.max_y, -1);
      }
    }
  }
}

void sweep_on_Y(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                std::map<std::pair<int, int>, int> &compress, int MAXX,
                int MAXY, std::set<int> &result) {

  std::vector<Interval> intervals;

  for (auto &it : lhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];

    Interval left_interval = Interval(min_y, min_x, max_x, it.polygon_id, true);
    Interval right_interval =
        Interval(max_y, min_x, max_x, it.polygon_id, false);

    intervals.push_back(left_interval);
    intervals.push_back(right_interval);
  }

  for (auto &it : rhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];

    Interval left_interval = Interval(min_y, min_x, max_x, it.polygon_id, true);
    Interval right_interval =
        Interval(max_y, min_x, max_x, it.polygon_id, false);

    intervals.push_back(left_interval);
    intervals.push_back(right_interval);
  }

  std::vector<int> B1;
  std::vector<int> B2;

  B1.assign(MAXX + 1, 0);
  B2.assign(MAXX + 1, 0);

  std::sort(intervals.begin(), intervals.end(),
            IntervalComparePositiveIdPriority());

  for (auto &it : intervals) {
    if (it.id > 0) {
      if (range_sum(B1, B2, it.min_y, it.max_y) > 0) {
        result.insert(it.id);
      }
    } else {
      if (it.is_start) {
        range_add(B1, B2, it.min_y, it.max_y, 1);
      } else {
        range_add(B1, B2, it.min_y, it.max_y, -1);
      }
    }
  }

  B1.assign(MAXX + 1, 0);
  B2.assign(MAXX + 1, 0);

  std::sort(intervals.begin(), intervals.end(),
            IntervalCompareNegativeIdPriority());

  for (auto &it : intervals) {
    if (it.id < 0) {
      if (range_sum(B1, B2, it.min_y, it.max_y) > 0) {
        result.insert(it.id);
      }
    } else {
      if (it.is_start) {
        range_add(B1, B2, it.min_y, it.max_y, 1);
      } else {
        range_add(B1, B2, it.min_y, it.max_y, -1);
      }
    }
  }
}

void mbr_no_inside(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                   std::map<std::pair<int, int>, int> &compress, int MAXX,
                   int MAXY, std::set<int> &result) {
  sweep_on_X(lhs, rhs, compress, MAXX, MAXY, result);
  sweep_on_Y(lhs, rhs, compress, MAXX, MAXY, result);
}
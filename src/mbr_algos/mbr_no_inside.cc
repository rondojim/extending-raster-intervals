#include "../../include/mbr_algos/mbr_no_inside.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/interval.h"
#include <algorithm>
#include <map>
#include <set>
#include <vector>

// Function to update a Binary Indexed Tree (BIT) at position x with value val
static void add(std::vector<int> &B, int x, int val) {
  for (; x < B.size(); x += x & -x) {
    B[x] += val;
  }
}

// Function to perform a range update on two BITs (B1 and B2) from l to r with
// value val
static void range_add(std::vector<int> &B1, std::vector<int> &B2, int l, int r,
                      int val) {
  add(B1, l, val);
  add(B1, r + 1, -val);
  add(B2, l, val * (l - 1));
  add(B2, r + 1, -val * r);
}

// Function to get the prefix sum from BIT B up to position x
static int sum(std::vector<int> &B, int x) {
  int res = 0;
  for (; x > 0; x -= x & -x) {
    res += B[x];
  }
  return res;
}

// Function to get the prefix sum up to position x using two BITs (B1 and B2)
static int prefix_sum(std::vector<int> &B1, std::vector<int> &B2, int x) {
  return sum(B1, x) * x - sum(B2, x);
}

// Function to get the range sum from l to r using two BITs (B1 and B2)
static int range_sum(std::vector<int> &B1, std::vector<int> &B2, int l, int r) {
  return prefix_sum(B1, B2, r) - prefix_sum(B1, B2, l - 1);
}

// Function to perform a sweep line algorithm along the X axis
static void sweep_on_X(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                       std::map<std::pair<int, int>, int> &compress, int MAXX,
                       int MAXY, std::set<int> &result) {

  std::vector<Interval> intervals; // Vector to store the intervals

  // Create intervals for lhs polygons
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

  // Create intervals for rhs polygons
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

  std::vector<int> B1; // BIT for range addition
  std::vector<int> B2; // BIT for range sum

  // Initialize BITs with size MAXY + 1
  B1.assign(MAXY + 1, 0);
  B2.assign(MAXY + 1, 0);

  // Sort intervals by their x-coordinate with positive id priority
  std::sort(intervals.begin(), intervals.end(),
            IntervalComparePositiveIdPriority());

  // Process the intervals
  for (auto &it : intervals) {
    if (it.id > 0) {
      // Query the BITs for range sum to check if there is an intersection
      if (range_sum(B1, B2, it.min_y, it.max_y) > 0) {
        result.insert(it.id);
      }
    } else {
      // Update the BITs for an active interval
      if (it.is_start) {
        range_add(B1, B2, it.min_y, it.max_y, 1);
      } else {
        range_add(B1, B2, it.min_y, it.max_y, -1);
      }
    }
  }

  // Reset the BITs
  B1.assign(MAXY + 1, 0);
  B2.assign(MAXY + 1, 0);

  // Sort intervals by their x-coordinate with negative id priority
  std::sort(intervals.begin(), intervals.end(),
            IntervalCompareNegativeIdPriority());

  // Process the intervals again
  for (auto &it : intervals) {
    if (it.id < 0) {
      // Query the BITs for range sum to check if there is an intersection
      if (range_sum(B1, B2, it.min_y, it.max_y) > 0) {
        result.insert(it.id);
      }
    } else {
      // Update the BITs for an active interval
      if (it.is_start) {
        range_add(B1, B2, it.min_y, it.max_y, 1);
      } else {
        range_add(B1, B2, it.min_y, it.max_y, -1);
      }
    }
  }
}

// Function to perform a sweep line algorithm along the Y axis
static void sweep_on_Y(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                       std::map<std::pair<int, int>, int> &compress, int MAXX,
                       int MAXY, std::set<int> &result) {

  std::vector<Interval> intervals; // Vector to store the intervals

  // Create intervals for lhs polygons
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

  // Create intervals for rhs polygons
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

  std::vector<int> B1; // BIT for range addition
  std::vector<int> B2; // BIT for range sum

  // Initialize BITs with size MAXX + 1
  B1.assign(MAXX + 1, 0);
  B2.assign(MAXX + 1, 0);

  // Sort intervals by their y-coordinate with positive id priority
  std::sort(intervals.begin(), intervals.end(),
            IntervalComparePositiveIdPriority());

  // Process the intervals
  for (auto &it : intervals) {
    if (it.id > 0) {
      // Query the BITs for range sum to check if there is an intersection
      if (range_sum(B1, B2, it.min_y, it.max_y) > 0) {
        result.insert(it.id);
      }
    } else {
      // Update the BITs for an active interval
      if (it.is_start) {
        range_add(B1, B2, it.min_y, it.max_y, 1);
      } else {
        range_add(B1, B2, it.min_y, it.max_y, -1);
      }
    }
  }

  // Reset the BITs
  B1.assign(MAXX + 1, 0);
  B2.assign(MAXX + 1, 0);

  // Sort intervals by their y-coordinate with negative id priority
  std::sort(intervals.begin(), intervals.end(),
            IntervalCompareNegativeIdPriority());

  // Process the intervals again
  for (auto &it : intervals) {
    if (it.id < 0) {
      // Query the BITs for range sum to check if there is an intersection
      if (range_sum(B1, B2, it.min_y, it.max_y) > 0) {
        result.insert(it.id);
      }
    } else {
      // Update the BITs for an active interval
      if (it.is_start) {
        range_add(B1, B2, it.min_y, it.max_y, 1);
      } else {
        range_add(B1, B2, it.min_y, it.max_y, -1);
      }
    }
  }
}

// Function to find MBR intersections without considering inner intersections
void mbr_no_inside(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                   std::map<std::pair<int, int>, int> &compress, int MAXX,
                   int MAXY, std::set<int> &result) {
  // Perform sweep line algorithm along X axis
  sweep_on_X(lhs, rhs, compress, MAXX, MAXY, result);
  // Perform sweep line algorithm along Y axis
  sweep_on_Y(lhs, rhs, compress, MAXX, MAXY, result);
}

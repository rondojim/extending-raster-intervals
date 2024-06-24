#include "../../include/mbr_algos/axis_intersection.h"
#include "../../include/utils/compressing.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/minimum_bounding_rectangle.h"
#include "../../include/utils/node.h"
#include "../../include/utils/progressbar.hpp"
#include <algorithm>
#include <chrono>
#include <future>
#include <limits>
#include <map>
#include <set>
#include <vector>

// Function to compare two MBRs
bool compare_mbr(const mbr &a, const mbr &b) {
  if (a.minc != b.minc)
    return a.minc < b.minc;

  if (a.maxc != b.maxc)
    return a.maxc < b.maxc;

  return a.id < b.id;
}

// Function to filter partitions and find intersections
static void filter_partition(std::vector<mbr> &mlhs, std::vector<mbr> &mrhs,
                             std::vector<std::pair<int, int>> &result_) {
  int lh = 0;
  int rh = 0;

  // Loop through both left-hand side (lhs) and right-hand side (rhs) MBRs that
  // are sorted by coordinate X
  while (lh < mlhs.size() && rh < mrhs.size()) {
    if (mlhs[lh].minc.first <= mrhs[rh].minc.first) {
      int ly1 = mlhs[lh].minc.second;
      int ly2 = mlhs[lh].maxc.second;

      int rh_ = rh;
      // Check if the MBRs intersect
      while (rh_ < mrhs.size() && mlhs[lh].maxc.first >= mrhs[rh_].minc.first) {
        int ry1 = mrhs[rh_].minc.second;
        int ry2 = mrhs[rh_].maxc.second;
        if (!(ly2 < ry1 || ly1 > ry2)) {
          result_.push_back({mlhs[lh].id, mrhs[rh_].id});
        }
        rh_++;
      }
      lh++;
    } else {
      int ry1 = mrhs[rh].minc.second;
      int ry2 = mrhs[rh].maxc.second;

      int lh_ = lh;
      // Check if the MBRs intersect
      while (lh_ < mlhs.size() && mrhs[rh].maxc.first >= mlhs[lh_].minc.first) {
        int ly1 = mlhs[lh_].minc.second;
        int ly2 = mlhs[lh_].maxc.second;
        if (!(ry2 < ly1 || ry1 > ly2)) {
          result_.push_back({mlhs[lh_].id, mrhs[rh].id});
        }
        lh_++;
      }
      rh++;
    }
  }
}

// Function to filter polygons using partitioning
static void filter_polygons(std::vector<mbr> &mlhs, std::vector<mbr> &mrhs,
                            int MAXX, int MAXY,
                            std::vector<std::pair<int, int>> &result,
                            int partitions = 3) {
  // Initialize partitions
  std::vector<int> Xpartitions;
  std::vector<int> Ypartitions;

  // make partitions for X axis
  for (int i = 0; i <= partitions; ++i) {
    Xpartitions.push_back(i * MAXX / partitions);
  }

  // make partitions for Y axis
  for (int i = 0; i <= partitions; ++i) {
    Ypartitions.push_back(i * MAXY / partitions);
  }

  // Partition the MBRs
  std::vector<mbr> lhs_parts[10][10];
  std::vector<mbr> rhs_parts[10][10];

  for (auto &it : mlhs) {
    int min_x = it.minc.first;
    int max_x = it.maxc.first;
    int min_y = it.minc.second;
    int max_y = it.maxc.second;

    std::vector<int> xids, yids;

    // a polygon might span through multiple partitions of the X axis
    for (int i = 0; i < partitions; ++i) {
      if (!(max_x < Xpartitions[i] || min_x > Xpartitions[i + 1])) {
        xids.push_back(i);
      }
    }

    // a polygon might span through multiple partitions of the Y axis
    for (int i = 0; i < partitions; ++i) {
      if (!(max_y < Ypartitions[i] || min_y > Ypartitions[i + 1])) {
        yids.push_back(i);
      }
    }

    // take all different combinations of X-partitions and Y-partitions in order
    // to find all the partition-cells the mbr belongs to
    for (auto xid : xids) {
      for (auto yid : yids) {
        lhs_parts[xid][yid].push_back(it);
      }
    }
  }

  // do the same for rhs mbr polygons
  for (auto &it : mrhs) {
    int min_x = it.minc.first;
    int max_x = it.maxc.first;
    int min_y = it.minc.second;
    int max_y = it.maxc.second;

    std::vector<int> xids, yids;

    for (int i = 0; i < partitions; ++i) {
      if (!(max_x < Xpartitions[i] || min_x > Xpartitions[i + 1])) {
        xids.push_back(i);
      }
    }

    for (int i = 0; i < partitions; ++i) {
      if (!(max_y < Ypartitions[i] || min_y > Ypartitions[i + 1])) {
        yids.push_back(i);
      }
    }

    for (auto xid : xids) {
      for (auto yid : yids) {
        rhs_parts[xid][yid].push_back(it);
      }
    }
  }

  // Sets to store IDs of intersecting polygons

  std::vector<std::pair<int, int>> result_[10][10];

  std::vector<std::future<void>> futures;

  // Launch asynchronous tasks for filtering partitions
  for (int i = 0; i < partitions; ++i) {
    for (int j = 0; j < partitions; ++j) {
      futures.emplace_back(std::async(std::launch::async, [&, i, j]() {
        filter_partition(lhs_parts[i][j], rhs_parts[i][j], result_[i][j]);
      }));
    }
  }

  // Wait for all tasks to complete
  for (auto &future : futures) {
    future.get();
  }

  // Filter partitions for the entire set
  // filter_partition(mlhs, mrhs, lhs_ids_, rhs_ids_);

  std::set<int> lhs_ids;
  std::set<int> rhs_ids;

  // Combine results from all partitions
  for (int i = 0; i < partitions; ++i) {
    for (int j = 0; j < partitions; ++j) {
      // push back in result the result_[i][j] vector
      for (auto &pair : result_[i][j]) {
        lhs_ids.insert(pair.first);
        rhs_ids.insert(pair.second);
        result.push_back(pair);
      }
    }
  }
  printf("Filtered: %d\n", lhs_ids.size() + rhs_ids.size());
  // Store results in the result pair
}

static void swap_axis(std::vector<mbr> &mlhs, std::vector<mbr> &mrhs) {

  for (auto &it : mlhs) {
    std::swap(it.minc.first, it.minc.second);
    std::swap(it.maxc.first, it.maxc.second);
  }
  for (auto &it : mrhs) {
    std::swap(it.minc.first, it.minc.second);
    std::swap(it.maxc.first, it.maxc.second);
  }
}

// Function to perform a forward scan for MBR intersection using a brute force
// approach
void forward_scan(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                  std::vector<std::pair<int, int>> &result, int partitions) {

  std::map<std::pair<int, int>, int> compress; // Map for coordinate compression

  // Perform coordinate compression and get the maximum compressed coordinates
  std::pair<int, int> coords = coordinate_compression(lhs, rhs, compress);
  int MAXX = coords.first, MAXY = coords.second;

  std::vector<mbr> mlhs;
  std::vector<mbr> mrhs;

  create_mbr_vectors(lhs, rhs, mlhs, mrhs, compress);

  long long x_inter = calculate_axis_intersection(mlhs, mrhs, true, MAXX);
  long long y_inter = calculate_axis_intersection(mlhs, mrhs, false, MAXY);

  if (x_inter > y_inter) {
    swap_axis(mlhs, mrhs);
    std::swap(MAXX, MAXY);
  }

  // Sort MBRs
  std::sort(mlhs.begin(), mlhs.end(), compare_mbr);
  std::sort(mrhs.begin(), mrhs.end(), compare_mbr);

  // Filter polygons using partitioning
  filter_polygons(mlhs, mrhs, MAXX, MAXY, result, partitions);

  // Store the resulting polygons in the result pai
}

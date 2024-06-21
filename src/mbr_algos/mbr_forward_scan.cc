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
                             std::set<int> &lhs_ids, std::set<int> &rhs_ids) {
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
          lhs_ids.insert(mlhs[lh].id);
          rhs_ids.insert(mrhs[rh_].id);
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
          lhs_ids.insert(mlhs[lh_].id);
          rhs_ids.insert(mrhs[rh].id);
        }
        lh_++;
      }
      rh++;
    }
  }
}

// Function to filter polygons using partitioning
static void
filter_polygons(std::vector<mbr> &mlhs, std::vector<mbr> &mrhs, int MAXX,
                int MAXY, std::pair<std::vector<int>, std::vector<int>> &result,
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
  std::set<int> lhs_ids_;
  std::set<int> rhs_ids_;

  std::set<int> lhs_ids[10][10];
  std::set<int> rhs_ids[10][10];

  std::vector<std::future<void>> futures;

  // Launch asynchronous tasks for filtering partitions
  for (int i = 0; i < partitions; ++i) {
    for (int j = 0; j < partitions; ++j) {
      futures.emplace_back(std::async(std::launch::async, [&, i, j]() {
        filter_partition(lhs_parts[i][j], rhs_parts[i][j], lhs_ids[i][j],
                         rhs_ids[i][j]);
      }));
    }
  }

  // Wait for all tasks to complete
  for (auto &future : futures) {
    future.get();
  }

  // Filter partitions for the entire set
  filter_partition(mlhs, mrhs, lhs_ids_, rhs_ids_);

  // Combine results from all partitions
  for (int i = 0; i < partitions; ++i) {
    for (int j = 0; j < partitions; ++j) {
      lhs_ids_.insert(lhs_ids[i][j].begin(), lhs_ids[i][j].end());
      rhs_ids_.insert(rhs_ids[i][j].begin(), rhs_ids[i][j].end());
    }
  }

  // Store results in the result pair
  for (auto &it : lhs_ids_) {
    result.first.push_back(it);
  }
  for (auto &it : rhs_ids_) {
    result.second.push_back(it);
  }
}

// Function to perform a forward scan for MBR intersection using a brute force
// approach
void forward_scan(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                  std::pair<std::vector<Polygon>, std::vector<Polygon>> &result,
                  int partitions) {

  std::map<std::pair<int, int>, int> compress; // Map for coordinate compression

  // Perform coordinate compression and get the maximum compressed coordinates
  std::pair<int, int> coords = coordinate_compression(lhs, rhs, compress);
  int MAXX = coords.first, MAXY = coords.second;

  std::vector<mbr> mlhs;
  std::vector<mbr> mrhs;

  create_mbr_vectors(lhs, rhs, mlhs, mrhs, compress);

  std::pair<std::vector<int>, std::vector<int>> result_;

  // Sort MBRs
  std::sort(mlhs.begin(), mlhs.end(), compare_mbr);
  std::sort(mrhs.begin(), mrhs.end(), compare_mbr);

  // Filter polygons using partitioning
  filter_polygons(mlhs, mrhs, MAXX, MAXY, result_, partitions);

  // Store the resulting polygons in the result pair
  for (auto &it : result_.first) {
    result.first.push_back(lhs[it - 1]);
  }
  for (auto &it : result_.second) {
    result.second.push_back(rhs[-it - 1]);
  }
}

#include "../../include/join_algos/node.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/progressbar.hpp"
#include <algorithm>
#include <chrono>
#include <future>
#include <limits>
#include <map>
#include <set>
#include <vector>

// #include <ext/pb_ds/assoc_container.hpp>
// #include <ext/pb_ds/tree_policy.hpp>

struct mbr {
  std::pair<int, int> minc;
  std::pair<int, int> maxc;
  int id;
};

bool compare_mbr(const mbr &a, const mbr &b) {
  if (a.minc != b.minc)
    return a.minc < b.minc;

  if (a.maxc != b.maxc)
    return a.maxc < b.maxc;

  return a.id < b.id;
}

static void filter_partition(std::vector<mbr> &mlhs, std::vector<mbr> &mrhs,
                             std::set<int> &lhs_ids, std::set<int> &rhs_ids) {
  int lh = 0;
  int rh = 0;

  while (lh < mlhs.size() && rh < mrhs.size()) {
    if (mlhs[lh].minc.first <= mrhs[rh].minc.first) {
      int ly1 = mlhs[lh].minc.second;
      int ly2 = mlhs[lh].maxc.second;

      int rh_ = rh;
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

static void
filter_polygons(std::vector<mbr> &mlhs, std::vector<mbr> &mrhs, int MAXX,
                int MAXY, std::pair<std::vector<int>, std::vector<int>> &result,
                int partitions = 3) {
  // std::vector<int> Xpartitions;
  // std::vector<int> Ypartitions;

  // for (int i = 0; i <= partitions; ++i)
  // {
  //     Xpartitions.push_back(i * MAXX / partitions);
  // }

  // for (int i = 0; i <= partitions; ++i)
  // {
  //     Ypartitions.push_back(i * MAXY / partitions);
  // }

  // std::vector<mbr> lhs_parts[10][10];
  // std::vector<mbr> rhs_parts[10][10];

  // for (auto &it : mlhs)
  // {
  //     int min_x = it.minc.first;
  //     int max_x = it.maxc.first;
  //     int min_y = it.minc.second;
  //     int max_y = it.maxc.second;

  //     std::vector<int> xids, yids;

  //     for (int i = 0; i < partitions; ++i)
  //     {
  //         if (!(max_x < Xpartitions[i] || min_x > Xpartitions[i + 1]))
  //         {
  //             xids.push_back(i);
  //         }
  //     }

  //     for (int i = 0; i < partitions; ++i)
  //     {
  //         if (!(max_y < Ypartitions[i] || min_y > Xpartitions[i + 1]))
  //         {
  //             yids.push_back(i);
  //         }
  //     }

  //     for (auto xid : xids)
  //     {
  //         for (auto yid : yids)
  //         {
  //             lhs_parts[xid][yid].push_back(it);
  //         }
  //     }
  // }

  // for (auto &it : mrhs)
  // {
  //     int min_x = it.minc.first;
  //     int max_x = it.maxc.first;
  //     int min_y = it.minc.second;
  //     int max_y = it.maxc.second;

  //     std::vector<int> xids, yids;

  //     for (int i = 0; i < partitions; ++i)
  //     {
  //         if (!(max_x < Xpartitions[i] || min_x > Xpartitions[i + 1]))
  //         {
  //             xids.push_back(i);
  //         }
  //     }

  //     for (int i = 0; i < partitions; ++i)
  //     {
  //         if (!(max_y < Ypartitions[i] || min_y > Xpartitions[i + 1]))
  //         {
  //             yids.push_back(i);
  //         }
  //     }

  //     for (auto xid : xids)
  //     {
  //         for (auto yid : yids)
  //         {
  //             rhs_parts[xid][yid].push_back(it);
  //         }
  //     }
  // }

  // for (int i = 0; i < partitions; ++i)
  // {
  //     for (int j = 0; j < 3; ++j)
  //     {
  //         printf("%zu ", lhs_parts[i][j].size());
  //     }
  //     printf("\n");
  // }

  // for (int i = 0; i < 3; ++i)
  // {
  //     for (int j = 0; j < 3; ++j)
  //     {
  //         printf("%zu ", rhs_parts[i][j].size());
  //     }
  //     printf("\n");
  // }

  std::set<int> lhs_ids_;
  std::set<int> rhs_ids_;

  // std::set<int> lhs_ids[10][10];
  // std::set<int> rhs_ids[10][10];

  // std::vector<std::future<void>> futures;

  // for (int i = 0; i < partitions; ++i)
  // {
  //     for (int j = 0; j < partitions; ++j)
  //     {
  //         futures.emplace_back(std::async(std::launch::async, [&, i, j]()
  //                                         { filter_partition(lhs_parts[i][j],
  //                                         rhs_parts[i][j], lhs_ids[i][j],
  //                                         rhs_ids[i][j]); }));
  //     }
  // }

  // Wait for all tasks to complete
  // for (auto &future : futures)
  // {
  //     future.get();
  // }

  filter_partition(mlhs, mrhs, lhs_ids_, rhs_ids_);

  // for (int i = 0; i < partitions; ++i)
  // {
  //     for (int j = 0; j < partitions; ++j)
  //     {
  //         lhs_ids_.insert(lhs_ids[i][j].begin(), lhs_ids[i][j].end());
  //         rhs_ids_.insert(rhs_ids[i][j].begin(), rhs_ids[i][j].end());
  //     }
  // }

  for (auto &it : lhs_ids_) {
    result.first.push_back(it);
  }
  for (auto &it : rhs_ids_) {
    result.second.push_back(it);
  }
}

void mbr_brute(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
               std::map<std::pair<int, int>, int> &compress, int MAXX, int MAXY,
               std::pair<std::vector<Polygon>, std::vector<Polygon>> &result) {
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

  std::pair<std::vector<int>, std::vector<int>> result_;

  std::sort(mlhs.begin(), mlhs.end(), compare_mbr);
  std::sort(mrhs.begin(), mrhs.end(), compare_mbr);

  filter_polygons(mlhs, mrhs, MAXX, MAXY, result_, 1);

  for (auto &it : result_.first) {
    result.first.push_back(lhs[it - 1]);
  }
  for (auto &it : result_.second) {
    result.second.push_back(rhs[-it - 1]);
  }
}
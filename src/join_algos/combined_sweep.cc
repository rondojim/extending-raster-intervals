#include "../../include/join_algos/combined_sweep.h"
#include "../../include/join_algos/interval.h"
#include "../../include/join_algos/mbr_no_inside.h"
#include "../../include/join_algos/mbr_no_outside.h"
#include "../../include/utils/geometry_types.h"
#include <algorithm>
#include <map>
#include <set>
#include <vector>

static void find_all_points(
    std::vector<Polygon> &polygons,
    std::vector<std::pair<double, std::pair<int, int>>> &all_points) {
  for (auto &it : polygons) {

    all_points.push_back({it.minCorner.x, {0, it.polygon_id}});
    all_points.push_back({it.minCorner.y, {1, it.polygon_id}});
    all_points.push_back({it.maxCorner.x, {2, it.polygon_id}});
    all_points.push_back({it.maxCorner.y, {3, it.polygon_id}});
  }
}

static std::pair<int, int>
coordinate_compression(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                       std::map<std::pair<int, int>, int> &compress) {
  std::vector<std::pair<double, std::pair<int, int>>> all_points;
  find_all_points(lhs, all_points);
  find_all_points(rhs, all_points);

  std::sort(all_points.begin(), all_points.end());

  double prev = all_points[0].first - 1;
  int idx = 0;
  int MAXX = 0, MAXY = 0;

  for (auto &it : all_points) {

    if (it.first > prev) {
      idx++;
    }
    if (it.second.first == 0 || it.second.first == 2) {
      MAXX = idx;
    } else {
      MAXY = idx;
    }
    compress[{it.second}] = idx;
    prev = it.first;
  }
  return {MAXX, MAXY};
}

void mbr_combined_sweep(
    std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
    std::pair<std::vector<Polygon>, std::vector<Polygon>> &final_result) {
  std::set<int> result;
  std::map<std::pair<int, int>, int> compress;

  std::pair<int, int> coords = coordinate_compression(lhs, rhs, compress);
  int MAXX = coords.first, MAXY = coords.second;

  mbr_no_outside(lhs, rhs, compress, MAXX, MAXY, result);
  mbr_no_inside(lhs, rhs, compress, MAXX, MAXY, result);

  printf("%zu\n", result.size());

  for (auto &it : result) {
    if (it > 0) {
      final_result.first.push_back(lhs[it - 1]);
    } else {
      final_result.second.push_back(rhs[-it - 1]);
    }
  }
}
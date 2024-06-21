#include "../../include/utils/compressing.h"

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

std::pair<int, int>
coordinate_compression(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                       std::map<std::pair<int, int>, int> &compress) {
  std::vector<std::pair<double, std::pair<int, int>>> all_points;
  find_all_points(lhs, all_points);
  find_all_points(rhs, all_points);

  // Sort all points to assign compressed coordinates
  std::sort(all_points.begin(), all_points.end());

  double prev = all_points[0].first - 1;
  int idx = 0;
  int MAXX = 0, MAXY = 0;

  // Assign indices based on sorted order to compress coordinates
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

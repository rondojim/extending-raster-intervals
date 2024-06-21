#include "../../include/utils/minimum_bounding_rectangle.h"

mbr::mbr(std::pair<int, int> min_corner, std::pair<int, int> max_corner,
         int polygon_id)
    : minc(min_corner), maxc(max_corner), id(polygon_id) {}

bool mbr::check_intersection(const mbr &other) const {
  if ((maxc.first < other.minc.first) || (minc.first > other.maxc.first) ||
      (maxc.second < other.minc.second) || (minc.second > other.maxc.second)) {
    return false;
  }
  return true;
}

void create_mbr_vectors(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                        std::vector<mbr> &mlhs, std::vector<mbr> &mrhs,
                        std::map<std::pair<int, int>, int> &compress) {

  // Compress coordinates for lhs polygons and store in mlhs
  for (auto &it : lhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];
    mbr cur({min_x, min_y}, {max_x, max_y}, it.polygon_id);
    mlhs.push_back(cur);
  }

  // Compress coordinates for rhs polygons and store in mrhs
  for (auto &it : rhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];
    mbr cur({min_x, min_y}, {max_x, max_y}, it.polygon_id);
    mrhs.push_back(cur);
  }
}
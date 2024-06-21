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
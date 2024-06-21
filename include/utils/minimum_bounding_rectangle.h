#ifndef MINIMUM_BOUNDING_RECTANGLE_H
#define MINIMUM_BOUNDING_RECTANGLE_H

#include <utility>

// Structure representing a Minimum Bounding Rectangle (MBR)
struct mbr {
  std::pair<int, int> minc; // Minimum corner coordinates (x, y)
  std::pair<int, int> maxc; // Maximum corner coordinates (x, y)
  int id;                   // Identifier for the polygon

  mbr(std::pair<int, int> min_corner, std::pair<int, int> max_corner,
      int polygon_id);

  bool check_intersection(const mbr &other) const;
};

#endif // MINIMUM_BOUNDING_RECTANGLE_H
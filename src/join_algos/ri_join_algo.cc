#include "../../include/join_algos/ri_join_algo.h"
#include "../../include/join_algos/serialize_polygon.h"
#include "../../include/utils/geometry_types.h"
#include <algorithm>
#include <iostream>
#include <vector>

static void truncate_bitset(boost::dynamic_bitset<> &bitmask, int start,
                            int end, int ls, int le) {

  printf("start: %d, end: %d, ls: %d, le: %d\n", start, end, ls, le);

  assert(start >= ls);
  assert(end <= le);

  bitmask >>= (start - ls) * 3;
  bitmask.resize(bitmask.size() - (start - ls) * 3);

  bitmask <<= (le - end) * 3;
  bitmask >>= (le - end) * 3;
  bitmask.resize(bitmask.size() - (le - end) * 3);
}

bool ri_join_pair(SerializedPolygon &lhs, SerializedPolygon &rhs, int lhs_idx,
                  int rhs_idx, std::vector<std::pair<int, int>> &indecisive) {
  bool overlap = false;
  int i = 0, j = 0;

  printf("Serialized polygon lhs:\n");
  for (auto &bitmask : lhs.bitmasks) {
    printf("start: %d, end: %d\n", bitmask.start, bitmask.end);
    std::cout << bitmask.bitmask << std::endl;
  }

  printf("Serialized polygon rhs:\n");
  for (auto &bitmask : rhs.bitmasks) {
    printf("start: %d, end: %d\n", bitmask.start, bitmask.end);
    std::cout << bitmask.bitmask << std::endl;
  }
  while (i < lhs.bitmasks.size() && j < rhs.bitmasks.size()) {
    printf("lhs: %d %d, rhs: %d %d\n", lhs.bitmasks[i].start,
           lhs.bitmasks[i].end, rhs.bitmasks[j].start, rhs.bitmasks[j].end);
    if (!(lhs.bitmasks[i].end < rhs.bitmasks[j].start ||
          rhs.bitmasks[j].end < lhs.bitmasks[i].start)) {
      overlap = true;

      // compute the overlap interval endpoints
      int start = std::max(lhs.bitmasks[i].start, rhs.bitmasks[j].start);
      int end = std::min(lhs.bitmasks[i].end, rhs.bitmasks[j].end);

      printf("start: %d, end: %d\n", start, end);

      boost::dynamic_bitset<> lhs_bitmask = lhs.bitmasks[i].bitmask;

      truncate_bitset(lhs_bitmask, start, end, lhs.bitmasks[i].start,
                      lhs.bitmasks[i].end);

      boost::dynamic_bitset<> rhs_bitmask = rhs.bitmasks[j].bitmask;
      truncate_bitset(rhs_bitmask, start, end, rhs.bitmasks[j].start,
                      rhs.bitmasks[j].end);

      // check if and of bitset is non zero
      std::cout << "lhs_bitmask: " << lhs_bitmask << std::endl;
      std::cout << "rhs_bitmask: " << rhs_bitmask << std::endl;

      if ((lhs_bitmask & rhs_bitmask).any()) {
        return true;
      }
    }
    if (lhs.bitmasks[i].end <= rhs.bitmasks[j].end) {
      i++;
    } else {
      j++;
    }
  }
  if (overlap) {
    indecisive.push_back({lhs_idx, rhs_idx});
  }
  return false;
}

void ri_join_algo(std::vector<SerializedPolygon> &lhs_serialized_polygons,
                  std::vector<SerializedPolygon> &rhs_serialized_polygons,
                  std::vector<std::pair<int, int>> &result,
                  std::vector<std::pair<int, int>> &indecisive) {

  for (auto &lhs_polygon : lhs_serialized_polygons) {
    for (auto &rhs_polygon : rhs_serialized_polygons) {
      if (ri_join_pair(lhs_polygon, rhs_polygon, lhs_polygon.polygon_id,
                       rhs_polygon.polygon_id, indecisive)) {
        result.push_back({lhs_polygon.polygon_id, rhs_polygon.polygon_id});
      }
    }
  }
}
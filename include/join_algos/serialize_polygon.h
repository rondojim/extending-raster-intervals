#ifndef SERIALIZE_POLYGON_H
#define SERIALIZE_POLYGON_H

#include "inter_filter.h"
#include <boost/dynamic_bitset.hpp>
#include <map>
#include <vector>

struct IntervalBitmask {
  int start, end;
  boost::dynamic_bitset<> bitmask;

  IntervalBitmask(int start, int end, boost::dynamic_bitset<> bitmask)
      : start(start), end(end), bitmask(bitmask) {}

  // make operator = for copying IntervalBitmask
  IntervalBitmask &operator=(const IntervalBitmask &rhs) {
    start = rhs.start;
    end = rhs.end;
    bitmask = rhs.bitmask;
    return *this;
  }
};

struct SerializedPolygon {
  std::vector<IntervalBitmask> bitmasks;
  int polygon_id;

  SerializedPolygon(std::vector<IntervalBitmask> bitmasks, int polygon_id)
      : bitmasks(bitmasks), polygon_id(polygon_id) {}
};

void serialize_polygons(
    std::vector<RasterPolygonInfo> &lhs_i_j_to_rpoly_info,
    std::vector<RasterPolygonInfo> &rhs_i_j_to_rpoly_info,
    std::map<std::pair<int, int>, int> &hilbert_map,
    std::vector<SerializedPolygon> &lhs_serialized_polygons,
    std::vector<SerializedPolygon> &rhs_serialized_polygons);

#endif // SERIALIZE_POLYGON_H
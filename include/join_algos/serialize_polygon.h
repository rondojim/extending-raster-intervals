#ifndef SERIALIZE_POLYGON_H
#define SERIALIZE_POLYGON_H

#include "rasterization.h"
#include <boost/dynamic_bitset.hpp>
#include <map>
#include <vector>

// Structure representing an interval with a bitmask
struct IntervalBitmask {
  int start, end;                  // Start and end of the interval
  boost::dynamic_bitset<> bitmask; // Bitmask for the interval

  // Constructor to initialize IntervalBitmask
  IntervalBitmask(int start, int end, boost::dynamic_bitset<> bitmask)
      : start(start), end(end), bitmask(bitmask) {}

  // Assignment operator for copying IntervalBitmask
  IntervalBitmask &operator=(const IntervalBitmask &rhs) {
    start = rhs.start;
    end = rhs.end;
    bitmask = rhs.bitmask;
    return *this;
  }
};

// Structure representing a serialized polygon
struct SerializedPolygon {
  std::vector<IntervalBitmask> bitmasks; // Vector of interval bitmasks
  int polygon_id;                        // Identifier for the polygon

  // Constructor to initialize SerializedPolygon
  SerializedPolygon(std::vector<IntervalBitmask> bitmasks, int polygon_id)
      : bitmasks(bitmasks), polygon_id(polygon_id) {}
};

// Function to serialize polygons for Raster Interval join
// @param lhs_i_j_to_rpoly_info: Vector of raster polygon info for the left-hand
// side collection
// @param rhs_i_j_to_rpoly_info: Vector of raster polygon info for the
// right-hand side collection
// @param hilbert_map: Map for Hilbert curve indexing
// @param lhs_serialized_polygons: Vector to store serialized polygons from the
// left-hand side
// @param rhs_serialized_polygons: Vector to store serialized polygons from the
// right-hand side
void serialize_polygons(
    std::vector<RasterPolygonInfo> &lhs_i_j_to_rpoly_info,
    std::vector<RasterPolygonInfo> &rhs_i_j_to_rpoly_info,
    std::map<std::pair<int, int>, int> &hilbert_map,
    std::vector<SerializedPolygon> &lhs_serialized_polygons,
    std::vector<SerializedPolygon> &rhs_serialized_polygons);

#endif // SERIALIZE_POLYGON_H

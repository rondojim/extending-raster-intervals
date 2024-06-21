#ifndef RI_JOIN_ALGO_H
#define RI_JOIN_ALGO_H

#include "inter_filter.h"
#include "serialize_polygon.h"
#include <vector>

// Function to perform a Raster Interval (RI) join on a pair of serialized
// polygons
// @param lhs: Serialized polygon from the left-hand side collection
// @param rhs: Serialized polygon from the right-hand side collection
// @param lhs_idx: Index of the left-hand side polygon
// @param rhs_idx: Index of the right-hand side polygon
// @param indecisive: Vector to store pairs of indices that require further
// processing
// @return: Boolean indicating if the polygons intersect
bool ri_join_pair(SerializedPolygon &lhs, SerializedPolygon &rhs, int lhs_idx,
                  int rhs_idx, std::vector<std::pair<int, int>> &indecisive);

// Function to perform a Raster Interval (RI) join on collections of serialized
// polygons
// @param lhs_serialized_polygons: Vector of serialized polygons from the
// left-hand side collection
// @param rhs_serialized_polygons: Vector of serialized polygons from the
// right-hand side collection
// @param result: Vector to store pairs of indices of intersecting polygons
// @param indecisive: Vector to store pairs of indices that require further
// processing
void ri_join_algo(std::vector<SerializedPolygon> &lhs_serialized_polygons,
                  std::vector<SerializedPolygon> &rhs_serialized_polygons,
                  std::vector<std::pair<int, int>> &result,
                  std::vector<std::pair<int, int>> &indecisive);

#endif // RI_JOIN_ALGO_H

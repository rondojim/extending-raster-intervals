#ifndef RI_JOIN_ALGO_H
#define RI_JOIN_ALGO_H

#include "inter_filter.h"
#include "serialize_polygon.h"
#include <vector>

bool ri_join_pair(SerializedPolygon &lhs, SerializedPolygon &rhs, int lhs_idx,
                  int rhs_idx, std::vector<std::pair<int, int>> &indecisive);

void ri_join_algo(std::vector<SerializedPolygon> &lhs_serialized_polygons,
                  std::vector<SerializedPolygon> &rhs_serialized_polygons,
                  std::vector<std::pair<int, int>> &result,
                  std::vector<std::pair<int, int>> &indecisive);

#endif // RI_JOIN_ALGO_H
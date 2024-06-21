#ifndef MBR_SUBOPTIMAL_2D_SEGT_H
#define MBR_SUBOPTIMAL_2D_SEGT_H

#include "../utils/geometry_types.h"
#include <map>
#include <set>
#include <vector>

// Function to perform a 2D segment tree MBR intersection
// @param lhs: Vector of polygons representing the left-hand side collection
// @param rhs: Vector of polygons representing the right-hand side collection
// @param result: Pair of vectors to store the intersecting polygons from lhs
// and rhs
void mbr_2D_segt(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                 std::pair<std::vector<Polygon>, std::vector<Polygon>> &result);

#endif // MBR_SUBOPTIMAL_2D_SEGT_H

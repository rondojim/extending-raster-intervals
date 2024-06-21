#ifndef COMBINED_SWEEP_H
#define COMBINED_SWEEP_H

#include "../utils/geometry_types.h"
#include <map>
#include <set>
#include <vector>

// Function to perform the combined sweep algorithm for MBR intersection
// @param lhs: Vector of polygons representing the left-hand side collection
// @param rhs: Vector of polygons representing the right-hand side collection
// @param final_result: Pair of vectors to store the final intersecting polygons
// from lhs and rhs
void mbr_combined_sweep(
    std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
    std::pair<std::vector<Polygon>, std::vector<Polygon>> &final_result);

#endif // COMBINED_SWEEP_H

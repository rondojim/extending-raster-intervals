#ifndef MBR_SWEEP_BRUTE_H
#define MBR_SWEEP_BRUTE_H

#include "../utils/geometry_types.h"
#include <map>
#include <set>
#include <vector>

// Function to perform a forward scan for MBR intersection using a brute force
// approach
// @param lhs: Vector of polygons representing the left-hand side collection
// @param rhs: Vector of polygons representing the right-hand side collection
// @param result: Pair of vectors to store the intersecting polygons from lhs
// and rhs
void forward_scan(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                  std::pair<std::vector<Polygon>, std::vector<Polygon>> &result,
                  int partitions = 1);

#endif // MBR_SWEEP_BRUTE_H

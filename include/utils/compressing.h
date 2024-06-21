#ifndef COMPRESSING_H
#define COMPRESSING_H

#include "../utils/geometry_types.h"
#include <algorithm>
#include <map>
#include <vector>

// Function to perform coordinate compression for a set of polygons
// @param lhs: Vector of polygons representing the left-hand side collection
// @param rhs: Vector of polygons representing the right-hand side collection
// @param compress: Map to store the compressed coordinates
// @return: A pair of integers representing the maximum compressed coordinates
// for X and Y
std::pair<int, int>
coordinate_compression(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                       std::map<std::pair<int, int>, int> &compress);

#endif // COMPRESSING_H

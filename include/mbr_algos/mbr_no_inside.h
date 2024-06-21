#ifndef MBR_NO_INSIDE_H
#define MBR_NO_INSIDE_H

#include "../utils/geometry_types.h"
#include <map>
#include <set>
#include <vector>

// Function to perform MBR intersection without detecting rectangles that are
// completely inside another
// @param lhs: Vector of polygons representing the left-hand side collection
// @param rhs: Vector of polygons representing the right-hand side collection
// @param compress: Map to compress and map coordinates for efficient storage
// @param MAXX: Maximum X coordinate after compression
// @param MAXY: Maximum Y coordinate after compression
// @param result: Set to store the IDs of intersecting polygons
void mbr_no_inside(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                   std::map<std::pair<int, int>, int> &compress, int MAXX,
                   int MAXY, std::set<int> &result);

#endif // MBR_NO_INSIDE_H

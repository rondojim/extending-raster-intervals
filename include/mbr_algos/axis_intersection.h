#ifndef AXIS_INTERSECTION_H
#define AXIS_INTERSECTION_H

#include "../utils/geometry_types.h"
#include "../utils/minimum_bounding_rectangle.h"

#include <map>
#include <set>
#include <vector>

long long calculate_axis_intersection(std::vector<mbr> &lhs,
                                      std::vector<mbr> &rhs, bool Xaxis,
                                      int MAXC);

#endif // MBR_NO_INSIDE_H

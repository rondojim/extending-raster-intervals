#ifndef MBR_NO_INSIDE_H
#define MBR_NO_INSIDE_H

#include "../utils/geometry_types.h"
#include <map>
#include <set>
#include <vector>

void mbr_no_inside(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                   std::map<std::pair<int, int>, int> &compress, int MAXX,
                   int MAXY, std::set<int> &result);

#endif // MBR_NO_INSIDE_H

#ifndef MBR_NO_OUTSIDE_H
#define MBR_NO_OUTSIDE_H

#include "../utils/geometry_types.h"
#include <map>
#include <set>
#include <vector>

void mbr_no_outside(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                    std::map<std::pair<int, int>, int> &compress, int MAXX,
                    int MAXY, std::set<int> &result);

#endif // MBR_NO_OUTSIDE_H

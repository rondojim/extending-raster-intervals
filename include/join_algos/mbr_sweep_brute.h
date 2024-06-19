#ifndef MBR_SWEEP_BRUTE_H
#define MBR_SWEEP_BRUTE_H

#include "../utils/geometry_types.h"
#include <map>
#include <set>
#include <vector>

void mbr_brute(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
               std::map<std::pair<int, int>, int> &compress, int MAXX, int MAXY,
               std::pair<std::vector<Polygon>, std::vector<Polygon>> &result);

#endif // MBR_SWEEP_BRUTE_H

#ifndef COMBINED_SWEEP_H
#define COMBINED_SWEEP_H

#include "../utils/geometry_types.h"
#include <map>
#include <set>
#include <vector>

void mbr_combined_sweep(
    std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
    std::pair<std::vector<Polygon>, std::vector<Polygon>> &final_result);

#endif // COMBINED_SWEEP_H
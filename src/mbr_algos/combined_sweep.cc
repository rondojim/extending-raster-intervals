#include "../../include/mbr_algos/combined_sweep.h"
#include "../../include/mbr_algos/mbr_no_inside.h"
#include "../../include/mbr_algos/mbr_no_outside.h"
#include "../../include/utils/compressing.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/interval.h"
#include <algorithm>
#include <map>
#include <set>
#include <vector>

// Function to perform the combined sweep algorithm for MBR intersection
void mbr_combined_sweep(
    std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
    std::pair<std::vector<Polygon>, std::vector<Polygon>> &final_result) {
  std::set<int> result; // Set to store the result of intersecting polygons
  std::map<std::pair<int, int>, int> compress; // Map for coordinate compression

  // Perform coordinate compression and get the maximum compressed coordinates
  std::pair<int, int> coords = coordinate_compression(lhs, rhs, compress);
  int MAXX = coords.first, MAXY = coords.second;

  // Perform MBR intersection algorithms
  mbr_no_outside(lhs, rhs, compress, MAXX, MAXY, result);
  mbr_no_inside(lhs, rhs, compress, MAXX, MAXY, result);

  // Print the number of results
  printf("%zu\n", result.size());

  // Populate the final result with intersecting polygons
  for (auto &it : result) {
    if (it > 0) {
      final_result.first.push_back(lhs[it - 1]);
    } else {
      final_result.second.push_back(rhs[-it - 1]);
    }
  }
}

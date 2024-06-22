#include "../../include/join_algos/ri_join_algo.h"
#include "../../include/join_algos/serialize_polygon.h"
#include "../../include/utils/geometry_types.h"
#include <algorithm>
#include <iostream>
#include <vector>

// Function to truncate a bitset to the specified range
static void truncate_bitset(boost::dynamic_bitset<> &bitmask, int start,
                            int end, int ls, int le) {

  // Ensure the start and end are within bounds
  assert(start >= ls);
  assert(end <= le);

  // Shift right to remove bits before the start of the interval
  bitmask >>= (start - ls) * 3;
  bitmask.resize(bitmask.size() - (start - ls) * 3);

  // Shift left and then right to remove bits after the end of the interval
  bitmask <<= (le - end) * 3;
  bitmask >>= (le - end) * 3;
  bitmask.resize(bitmask.size() - (le - end) * 3);
}

// Function to join two serialized polygons and check for overlaps
bool ri_join_pair(SerializedPolygon &lhs, SerializedPolygon &rhs, int lhs_idx,
                  int rhs_idx, std::set<std::pair<int, int>> &indecisive) {
  bool overlap = false; // Flag to indicate if there is any overlap
  int i = 0, j = 0;     // Indices for iterating over bitmasks

  // Iterate through the bitmasks of both polygons to check for overlaps
  while (i < lhs.bitmasks.size() && j < rhs.bitmasks.size()) {

    // Check if the current intervals overlap
    if (!(lhs.bitmasks[i].end < rhs.bitmasks[j].start ||
          rhs.bitmasks[j].end < lhs.bitmasks[i].start)) {
      overlap = true;

      // Compute the overlap interval endpoints
      int start = std::max(lhs.bitmasks[i].start, rhs.bitmasks[j].start);
      int end = std::min(lhs.bitmasks[i].end, rhs.bitmasks[j].end);

      // Truncate the bitsets to the overlap interval
      boost::dynamic_bitset<> lhs_bitmask = lhs.bitmasks[i].bitmask;
      truncate_bitset(lhs_bitmask, start, end, lhs.bitmasks[i].start,
                      lhs.bitmasks[i].end);

      boost::dynamic_bitset<> rhs_bitmask = rhs.bitmasks[j].bitmask;
      truncate_bitset(rhs_bitmask, start, end, rhs.bitmasks[j].start,
                      rhs.bitmasks[j].end);

      // Check if the AND of the bitsets is non-zero, indicating an overlap
      if ((lhs_bitmask & rhs_bitmask).any()) {
        return true;
      }
    }

    // Move to the next interval in the bitmask that ends first
    if (lhs.bitmasks[i].end <= rhs.bitmasks[j].end) {
      i++;
    } else {
      j++;
    }
  }

  // If there was any overlap, but did not return true, add the pair to the
  // indecisive list
  if (overlap) {
    indecisive.insert({lhs_idx, rhs_idx});
  }

  return false;
}

// Function to join all serialized polygons from two collections and store the
// results
void ri_join_algo(std::vector<SerializedPolygon> &lhs_serialized_polygons,
                  std::vector<SerializedPolygon> &rhs_serialized_polygons,
                  std::set<std::pair<int, int>> &result,
                  std::set<std::pair<int, int>> &indecisive) {

  // Iterate through all pairs of polygons from both collections
  for (auto &lhs_polygon : lhs_serialized_polygons) {
    for (auto &rhs_polygon : rhs_serialized_polygons) {
      // Check if the current pair of polygons overlap
      if (ri_join_pair(lhs_polygon, rhs_polygon, lhs_polygon.polygon_id,
                       rhs_polygon.polygon_id, indecisive)) {
        // If they overlap, add the pair to the result vector
        result.insert({lhs_polygon.polygon_id, rhs_polygon.polygon_id});
      }
    }
  }
}

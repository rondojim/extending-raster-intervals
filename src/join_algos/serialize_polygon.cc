#include "../../include/join_algos/serialize_polygon.h"
#include "../../include/join_algos/rasterization.h"
#include "../../include/utils/geometry_types.h"
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <functional>
#include <map>
#include <vector>

// Function to serialize polygons based on their raster polygon info and a
// Hilbert map
void serialize_polygons(
    std::vector<RasterPolygonInfo> &lhs_i_j_to_rpoly_info,
    std::vector<RasterPolygonInfo> &rhs_i_j_to_rpoly_info,
    std::map<std::pair<int, int>, int> &hilbert_map,
    std::vector<SerializedPolygon> &lhs_serialized_polygons,
    std::vector<SerializedPolygon> &rhs_serialized_polygons) {

  // Create vectors of references to the input and output data structures for
  // easier iteration
  std::vector<std::reference_wrapper<std::vector<RasterPolygonInfo>>> infomaps =
      {lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info};
  std::vector<std::reference_wrapper<std::vector<SerializedPolygon>>>
      serialized_polygons = {lhs_serialized_polygons, rhs_serialized_polygons};

  // Loop through the two sets of polygon info (lhs and rhs)
  for (int i = 0; i < 2; i++) {

    // Iterate through each RasterPolygonInfo in the current set
    for (auto &cur_polygon_info : infomaps[i].get()) {
      // Vector to store Hilbert values and their corresponding BinaryCellCodes
      std::vector<std::pair<int, BinaryCellCode>> hilbert_values;

      // Iterate through each (i, j) pair in the current polygon info
      for (auto &cur_pair : cur_polygon_info.i_j_to_rcell_info) {
        // Find the corresponding Hilbert value for the current (i, j) pair
        int hilbert_val =
            hilbert_map[{cur_pair.first.first, cur_pair.first.second}];

        // Add the Hilbert value and BinaryCellCode to the vector
        hilbert_values.push_back({hilbert_val, cur_pair.second.cell_type});
      }

      // Sort the vector of Hilbert values based on the Hilbert value
      std::sort(hilbert_values.begin(), hilbert_values.end(),
                [](const std::pair<int, BinaryCellCode> &a,
                   const std::pair<int, BinaryCellCode> &b) {
                  return a.first < b.first;
                });

      // Vector to store the intervals of Hilbert values and their corresponding
      // bitmasks
      std::vector<IntervalBitmask> bitmasks;
      int start = -10, end = -10; // Initialize start and end of the interval
      boost::dynamic_bitset<> bitmask; // Initialize the bitmask

      // Iterate through the sorted Hilbert values to create intervals and
      // bitmasks
      for (auto &cur_pair : hilbert_values) {
        // Check if the current Hilbert value extends the current interval
        if (end + 1 == cur_pair.first) {
          end = cur_pair.first;
        } else {
          // If the interval ends, push the current interval to the vector
          if (end != -10) {
            bitmasks.push_back({start, end, bitmask});
          }
          // Start a new interval
          start = cur_pair.first;
          end = cur_pair.first;
          bitmask.clear(); // Clear the bitmask for the new interval
        }

        // If processing the second set of polygons, XOR the bitset with 110
        if (i == 1) {
          std::bitset<3> temp(6);
          cur_pair.second.value ^= temp;
        }

        // Append the bits from the BinaryCellCode to the bitmask
        for (int j = 0; j < cur_pair.second.value.size(); ++j) {
          bitmask.push_back(cur_pair.second.value[j]);
        }
      }

      // Push the last interval to the vector
      bitmasks.push_back({start, end, bitmask});

      // Add the serialized polygon to the output vector
      serialized_polygons[i].get().emplace_back(bitmasks, cur_polygon_info.idx);
    }
  }
}

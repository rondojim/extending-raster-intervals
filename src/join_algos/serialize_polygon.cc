#include "../../include/join_algos/serialize_polygon.h"
#include "../../include/join_algos/inter_filter.h"
#include "../../include/utils/geometry_types.h"
#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <functional>
#include <map>
#include <vector>

void serialize_polygons(
    std::vector<RasterPolygonInfo> &lhs_i_j_to_rpoly_info,
    std::vector<RasterPolygonInfo> &rhs_i_j_to_rpoly_info,
    std::map<std::pair<int, int>, int> &hilbert_map,
    std::vector<SerializedPolygon> &lhs_serialized_polygons,
    std::vector<SerializedPolygon> &rhs_serialized_polygons) {

  std::vector<std::reference_wrapper<std::vector<RasterPolygonInfo>>> infomaps =
      {lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info};
  std::vector<std::reference_wrapper<std::vector<SerializedPolygon>>>
      serialized_polygons = {lhs_serialized_polygons, rhs_serialized_polygons};

  for (int i = 0; i < 2; i++) {
    // loop through each map of infomaps
    printf("in raster polygon info\n");
    for (auto &cur_polygon_info : infomaps[i].get()) {
      // create a vector of hilbert values for each pair in the map and store it
      // along with the BinaryCellCode of the RasterCellInfo
      std::vector<std::pair<int, BinaryCellCode>> hilbert_values;
      // for each pair i, j in map find the corresponding hilbert value
      for (auto &cur_pair : cur_polygon_info.i_j_to_rcell_info) {
        // find the corresponding hilbert value

        printf("i: %d, j: %d\n", cur_pair.first.first, cur_pair.first.second);

        int hilbert_val =
            hilbert_map[{cur_pair.first.first, cur_pair.first.second}];
        // add the hilbert value to the vector
        hilbert_values.push_back({hilbert_val, cur_pair.second.cell_type});
      }

      // sort the vector of hilbert values based on the hilbert value
      std::sort(hilbert_values.begin(), hilbert_values.end(),
                [](const std::pair<int, BinaryCellCode> &a,
                   const std::pair<int, BinaryCellCode> &b) {
                  return a.first < b.first;
                });
      // print hilbert values of the current polygon
      // for (auto &cur_pair : hilbert_values) {
      //   std::cout << cur_pair.first << " " << cur_pair.second.to_type()
      //             << std::endl;
      // }

      // for each continuous interval of hilbert values, create a bitset which
      // is the union of all the BinaryCellCodes in that interval
      std::vector<IntervalBitmask> bitmasks;
      int start = -10, end = -10;
      boost::dynamic_bitset<> bitmask;
      for (auto &cur_pair : hilbert_values) {

        if (end + 1 == cur_pair.first) {
          end = cur_pair.first;
        } else {
          // interval ends here so push already created interval to the vector
          // if end is not -10
          if (end != -10) {
            bitmasks.push_back({start, end, bitmask});
          }
          // start a new interval
          start = cur_pair.first;
          end = cur_pair.first;
          bitmask.clear();
        }

        if (i == 1) {
          // xor bitset cur_pair.second.value with bitset 110
          std::bitset<3> temp(6);
          cur_pair.second.value ^= temp;
        }

        for (int j = 0; j < cur_pair.second.value.size(); ++j) {
          bitmask.push_back(cur_pair.second.value[j]);
        }
      }
      bitmasks.push_back({start, end, bitmask});
      serialized_polygons[i].get().emplace_back(bitmasks, cur_polygon_info.idx);
    }
  }
}
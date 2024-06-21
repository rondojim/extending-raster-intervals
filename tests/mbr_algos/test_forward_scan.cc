#include "../../include/mbr_algos/mbr_forward_scan.h"
#include "../../include/utils/data_reader.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/progressbar.hpp"

#include <algorithm>
#include <chrono>
#include <map>
#include <random>
#include <set>
#include <time.h>
#include <vector>

static void id_polygons(std::vector<Polygon> &polygons, bool lhs) {
  int id = 1;
  for (auto &it : polygons) {
    it.polygon_id = (lhs ? id : -id);
    id++;
  }
}

int main() {
  std::vector<Polygon> lhs;
  std::vector<Polygon> rhs;

  std::string filename_lhs("../../dataset_files/OSM_by_continent/O5AF");
  std::string filename_rhs("../../dataset_files/OSM_by_continent/O6AF");

  Point lhsMinCorner, lhsMaxCorner, rhsMinCorner, rhsMaxCorner;

  lhs = read_data_find_MBR(filename_lhs, lhsMinCorner, lhsMaxCorner);
  printf("Read first set\n");

  id_polygons(lhs, true);

  rhs = read_data_find_MBR(filename_rhs, rhsMinCorner, rhsMaxCorner);

  id_polygons(rhs, false);

  std::vector<Polygon> lhs_copy = lhs;
  std::vector<Polygon> rhs_copy = rhs;

  std::pair<std::vector<Polygon>, std::vector<Polygon>> result;
  forward_scan(lhs, rhs, result, 1);

  printf("With 1 partition: %zu\n", result.first.size() + result.second.size());

  std::pair<std::vector<Polygon>, std::vector<Polygon>> result2;
  forward_scan(lhs_copy, rhs_copy, result2, 3);

  printf("With 3 partitions: %zu\n",
         result2.first.size() + result2.second.size());

  return 0;
}
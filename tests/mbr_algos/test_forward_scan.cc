#include "../../include/mbr_algos/mbr_forward_scan.h"
#include "../../include/utils/compressing.h"
#include "../../include/utils/data_reader.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/minimum_bounding_rectangle.h"
#include "../../include/utils/progressbar.hpp"
#include <algorithm>
#include <cassert>
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

  std::string filename_lhs("../../dataset_files/OSM_filtered_datasets/O5OC");
  std::string filename_rhs("../../dataset_files/OSM_filtered_datasets/O6OC");

  Point lhsMinCorner, lhsMaxCorner, rhsMinCorner, rhsMaxCorner;

  lhs = read_data_find_MBR(filename_lhs, lhsMinCorner, lhsMaxCorner);

  id_polygons(lhs, true);

  rhs = read_data_find_MBR(filename_rhs, rhsMinCorner, rhsMaxCorner);

  id_polygons(rhs, false);

  std::vector<mbr> mlhs;
  std::vector<mbr> mrhs;

  std::map<std::pair<int, int>, int> compress;
  coordinate_compression(lhs, rhs, compress);

  std::vector<std::pair<int, int>> result;
  forward_scan(lhs, rhs, result, 3);

  return 0;
}
#include "../../include/mbr_algos/mbr_2d_segt.h"
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

  std::string filename_lhs("../../dataset_files/OSM_by_continent/O5OC");
  std::string filename_rhs("../../dataset_files/OSM_by_continent/O6OC");

  Point lhsMinCorner, lhsMaxCorner, rhsMinCorner, rhsMaxCorner;

  lhs = read_data_find_MBR(filename_lhs, lhsMinCorner, lhsMaxCorner);

  id_polygons(lhs, true);

  rhs = read_data_find_MBR(filename_rhs, rhsMinCorner, rhsMaxCorner);

  id_polygons(rhs, false);

  std::vector<mbr> mlhs;
  std::vector<mbr> mrhs;

  std::map<std::pair<int, int>, int> compress;
  coordinate_compression(lhs, rhs, compress);

  std::pair<std::vector<Polygon>, std::vector<Polygon>> result;
  mbr_2D_segt(lhs, rhs, result);

  create_mbr_vectors(result.first, result.second, mlhs, mrhs, compress);

  // test for every mbr in mlhs if it intersects with at least one mbr in mrhs
  for (auto &lhs_mbr : mlhs) {
    bool intersects = false;
    for (auto &rhs_mbr : mrhs) {
      if (lhs_mbr.check_intersection(rhs_mbr)) {
        intersects = true;
        break;
      }
    }
    assert(intersects && "Intersection not found");
  }

  // test for every mbr in mrhs if it intersects with at least one mbr in mlhs
  for (auto &rhs_mbr : mrhs) {
    bool intersects = false;
    for (auto &lhs_mbr : mlhs) {
      if (rhs_mbr.check_intersection(lhs_mbr)) {
        intersects = true;
        break;
      }
    }
    assert(intersects && "Intersection not found");
  }

  assert((result.first.size() + result.second.size()) == 48891 &&
         "Incorrect number of polygons");
  return 0;
}
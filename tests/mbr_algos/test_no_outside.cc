#include "../../include/mbr_algos/mbr_no_outside.h"
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
  std::pair<int, int> coords = coordinate_compression(lhs, rhs, compress);
  int MAXX = coords.first;
  int MAXY = coords.second;

  std::pair<std::vector<Polygon>, std::vector<Polygon>> result;
  std::set<int> resultset;
  mbr_no_outside(lhs, rhs, compress, MAXX, MAXY, resultset);

  for (auto &p : resultset) {
    if (p > 0)
      result.first.push_back(lhs[p - 1]);
    else
      result.second.push_back(rhs[-p - 1]);
  }

  create_mbr_vectors(lhs, rhs, mlhs, mrhs, compress);

  // for every polygon in result.first, check if there exists an intersection in
  // mrhs
  for (auto &lhs_mbr : result.first) {
    bool intersects = false;
    for (auto &rhs_mbr : mrhs) {
      if (mlhs[lhs_mbr.polygon_id - 1].check_intersection(rhs_mbr)) {
        intersects = true;
        break;
      }
    }
    assert(intersects && "Intersection not found");
  }
  // for every polygon in result.second, check if there exists an intersection
  // in mlhs
  for (auto &rhs_mbr : result.second) {
    bool intersects = false;
    for (auto &lhs_mbr : mlhs) {
      if (mrhs[-rhs_mbr.polygon_id - 1].check_intersection(lhs_mbr)) {
        intersects = true;
        break;
      }
    }
    assert(intersects && "Intersection not found");
  }

  return 0;
}
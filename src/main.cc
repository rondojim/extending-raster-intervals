#include "../include/join_algos/combined_sweep.h"
#include "../include/join_algos/inter_filter.h"
#include "../include/join_algos/ri_join_algo.h"
#include "../include/join_algos/serialize_polygon.h"
#include "../include/utils/data_reader.h"
#include "../include/utils/geometry_types.h"
#include "../include/utils/hilbert.h"
#include <algorithm>
#include <chrono>
#include <set>
#include <vector>

static void id_polygons(std::vector<Polygon> &polygons, bool lhs) {
  int id = 1;
  for (auto &it : polygons) {
    it.polygon_id = (lhs ? id : -id);
    id++;
  }
}

static void get_preprocessed_polygons(std::vector<Polygon> &lhs,
                                      std::vector<Polygon> &rhs,
                                      const std::string &filename_lhs,
                                      const std::string &filename_rhs,
                                      Point &gridMinCorner,
                                      Point &gridMaxCorner, int n = -1) {
  Point lhsMinCorner, lhsMaxCorner, rhsMinCorner, rhsMaxCorner;
  lhs = read_data_find_MBR(filename_lhs, lhsMinCorner, lhsMaxCorner, n);
  printf("Read first set\n");

  id_polygons(lhs, true);

  rhs = read_data_find_MBR(filename_rhs, rhsMinCorner, rhsMaxCorner, n);

  id_polygons(rhs, false);

  printf("Read second set\n");

  gridMinCorner = {std::min(lhsMinCorner.x, rhsMinCorner.x),
                   std::min(lhsMinCorner.y, rhsMinCorner.y)};

  gridMaxCorner = {std::max(lhsMaxCorner.x, rhsMaxCorner.x),
                   std::max(lhsMaxCorner.y, rhsMaxCorner.y)};
}

int find_interesctions(
    unsigned int N, std::string lhs_f_name, std::string rhs_f_name,
    std::pair<std::vector<Polygon>, std::vector<Polygon>> intersections) {
  std::vector<Polygon> lhs_polygons;
  std::vector<Polygon> rhs_polygons;
  Point gridMinCorner, gridMaxCorner;

  get_preprocessed_polygons(lhs_polygons, rhs_polygons, lhs_f_name, rhs_f_name,
                            gridMinCorner, gridMaxCorner, -1);

  std::vector<Polygon> lhs_polygons_copy = lhs_polygons;
  std::vector<Polygon> rhs_polygons_copy = rhs_polygons;

  // printf all polygon ids in lhs_polygons
  for (auto &p : lhs_polygons) {
    printf("Polygon id: %d\n", p.polygon_id);
  }
  // printf all polygon ids in rhs_polygons
  for (auto &p : rhs_polygons) {
    printf("Polygon id: %d\n", p.polygon_id);
  }

  RasterGrid grid = RasterGrid(N, gridMinCorner, gridMaxCorner);
  std::cout << grid.min_corner.to_str() << std::endl;
  std::cout << grid.max_corner.to_str() << std::endl;
  // run mbr combined sweep ------------------------------------------------
  std::pair<std::vector<Polygon>, std::vector<Polygon>> final_result;
  printf("%s %s %d\n", lhs_polygons[0].minCorner.to_str().c_str(),
         lhs_polygons[0].maxCorner.to_str().c_str(),
         lhs_polygons[0].polygon_id);
  printf("%s %s %d\n", rhs_polygons[0].minCorner.to_str().c_str(),
         rhs_polygons[0].maxCorner.to_str().c_str(),
         rhs_polygons[0].polygon_id);
  mbr_combined_sweep(lhs_polygons, rhs_polygons, final_result);
  printf("%d %d\n", final_result.first.size(), final_result.second.size());

  printf("Combined sweep done\n");
  // // run interval filter ----------------------------------------------------
  // RasterGrid grid = RasterGrid(N, gridMinCorner, gridMaxCorner);

  std::vector<RasterPolygonInfo> lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info;

  if (!rasterize_polygons(grid, final_result.first, final_result.second,
                          lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info,
                          "bin_error.txt", false)) {
    printf("Error in rasterization\n");
    return -1;
  }

  printf("Rasterization done\n");

  std::vector<SerializedPolygon> lhs_serialized_polygons;
  std::vector<SerializedPolygon> rhs_serialized_polygons;

  // load hilbert map
  std::map<std::pair<int, int>, int> hilbert_map;

  std::string base_path =
      "/Users/dimronto/Desktop/DSIT/Spring/Database "
      "Systems/project/phase1-paper-presentation-roadmap-2024-raster-intervals/"
      "hilbert_maps/hilbert_";
  std::string hilbert_map_f_name = base_path + std::to_string(N) + ".txt";

  loadMapFromFile(hilbert_map, hilbert_map_f_name);
  printf("Hilbert map loaded\n");
  serialize_polygons(lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info, hilbert_map,
                     lhs_serialized_polygons, rhs_serialized_polygons);
  printf("Serialization done\n");

  std::vector<std::pair<int, int>> result;
  std::vector<std::pair<int, int>> indecisive;

  printf("Starting RI join\n");
  ri_join_algo(lhs_serialized_polygons, rhs_serialized_polygons, result,
               indecisive);

  // print result
  printf("Result\n");
  for (auto &r : result) {
    printf("lhs polygon id: %d, rhs polygon id: %d\n", r.first, r.second);
  }
  // print indecisive
  printf("Indecisive\n");
  for (auto &r : indecisive) {
    printf("lhs polygon id: %d, rhs polygon id: %d\n", r.first, r.second);
  }

  // check indecisive
  for (auto &r : indecisive) {
    Polygon lhs_p = lhs_polygons_copy[r.first - 1];
    Polygon rhs_p = rhs_polygons_copy[-r.second - 1];

    if (lhs_p.intersects(rhs_p)) {
      printf("Pair (%d, %d) intersects \n", r.first, r.second);
      result.push_back(r);
    }
  }

  // print lhs serialized polygons
  // for (auto &sp : lhs_serialized_polygons) {
  //   printf("Polygon id: %d\n", sp.polygon_id);
  //   for (auto &ib : sp.bitmasks) {
  //     printf("Interval: %d %d\n", ib.start, ib.end);
  //     printf("Bitmask: ");
  //     std::cout << "Bitmask: " << ib.bitmask << std::endl;
  //   }
  // }
  // print rhs serialized polygons
  // for (auto &sp : rhs_serialized_polygons) {
  //   printf("Polygon id: %d\n", sp.polygon_id);
  //   for (auto &ib : sp.bitmasks) {
  //     printf("Interval: %d %d\n", ib.start, ib.end);
  //     printf("Bitmask: ");
  //     std::cout << "Bitmask: " << ib.bitmask << std::endl;
  //   }
  // }

  return 0;
}

int main() {
  unsigned int N = 2;
  // std::string rhs_f_name("../../dataset_files/OSM_by_continent/O5OC");
  // std::string lhs_f_name("../../dataset_files/OSM_by_continent/O6OC");
  std::string lhs_f_name("../datasets/O5AF.txt");
  std::string rhs_f_name("../datasets/O6AF.txt");
  std::pair<std::vector<Polygon>, std::vector<Polygon>> intersections;

  int result = find_interesctions(N, lhs_f_name, rhs_f_name, intersections);
  // get polygons -----------------------------------------------------------

  return 0;
}
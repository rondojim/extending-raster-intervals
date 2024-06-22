#include "../../include/join_algos/inter_filter.h"
#include "../../include/join_algos/ri_join_algo.h"
#include "../../include/join_algos/serialize_polygon.h"
#include "../../include/mbr_algos/combined_sweep.h"
#include "../../include/utils/compressing.h"
#include "../../include/utils/data_reader.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/hilbert.h"
#include "../../include/utils/minimum_bounding_rectangle.h"
#include <algorithm>
#include <chrono>
#include <set>
#include <vector>
#include "../../include/utils/pipeline.h"

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
                                      Point &gridMaxCorner, int n) {
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

void find_interesctions(
    unsigned int N, std::string lhs_f_name, std::string rhs_f_name,
    std::pair<std::vector<Polygon>, std::vector<Polygon>> &intersections) {
  std::vector<Polygon> lhs_polygons;
  std::vector<Polygon> rhs_polygons;
  Point gridMinCorner, gridMaxCorner;

  get_preprocessed_polygons(lhs_polygons, rhs_polygons, lhs_f_name, rhs_f_name,
                            gridMinCorner, gridMaxCorner);

  std::vector<Polygon> lhs_polygons_copy = lhs_polygons;
  std::vector<Polygon> rhs_polygons_copy = rhs_polygons;

  RasterGrid grid = RasterGrid(N, gridMinCorner, gridMaxCorner);

  // run mbr combined sweep ------------------------------------------------
  std::pair<std::vector<Polygon>, std::vector<Polygon>> final_result;
  mbr_combined_sweep(lhs_polygons, rhs_polygons, final_result);
  printf("%d %d\n", final_result.first.size(), final_result.second.size());

  std::vector<RasterPolygonInfo> lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info;

  std::set<int> null_cell_code_poly_idxs, error_poly_idxs;
  if (rasterize_polygons(grid, final_result.first, final_result.second,
                         lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info,
                         null_cell_code_poly_idxs, error_poly_idxs,
                         "bin_error.txt", false)) {
    printf("Error in rasterization\n");
  }

  printf("Ratio of wrong polygons: %f\n",
         (double)(null_cell_code_poly_idxs.size() + error_poly_idxs.size()) /
             (final_result.first.size() + final_result.second.size()));

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

  std::set<std::pair<int, int>> result;
  std::set<std::pair<int, int>> indecisive;

  printf("Starting RI join\n");
  ri_join_algo(lhs_serialized_polygons, rhs_serialized_polygons, result,
               indecisive);

  printf("Indecisive ratio: %f\n",
         (double)indecisive.size() /
             (final_result.first.size() * final_result.second.size()));

  printf("Ended RI join\n");
  // check indecisive
  for (auto &r : indecisive) {
    Polygon lhs_p = lhs_polygons_copy[r.first - 1];
    Polygon rhs_p = rhs_polygons_copy[-r.second - 1];

    if (lhs_p.intersects(rhs_p)) {
      result.insert(r);
    }
  }
  printf("Found indecisive\n");

  // create a set with the pairs that intersect
  std::set<std::pair<int, int>> result_set;
  for (auto &r : result) {
    result_set.insert(r);
  }

  // for each pair of polygon from final_result check if they intersect
  for (int i = 0; i < final_result.first.size(); i++) {
    // check if id in null_cell_code_poly_idxs
    if (null_cell_code_poly_idxs.find(final_result.first[i].polygon_id) !=
        null_cell_code_poly_idxs.end()) {
      continue;
    }
    // check if id in error_poly_idxs
    if (error_poly_idxs.find(final_result.first[i].polygon_id) !=
        error_poly_idxs.end()) {
      continue;
    }
    for (int j = 0; j < final_result.second.size(); j++) {
      if (null_cell_code_poly_idxs.find(final_result.first[j].polygon_id) !=
          null_cell_code_poly_idxs.end()) {
        continue;
      }
      if (error_poly_idxs.find(final_result.first[j].polygon_id) !=
          error_poly_idxs.end()) {
        continue;
      }
      Polygon lhs_p = lhs_polygons_copy[final_result.first[i].polygon_id - 1];
      Polygon rhs_p = rhs_polygons_copy[-final_result.second[j].polygon_id - 1];

      if (lhs_p.intersects(rhs_p)) {
        std::pair<int, int> p = {final_result.first[i].polygon_id,
                                 final_result.second[j].polygon_id};
        if (result_set.find(p) == result_set.end()) {
          printf("ERROR ids: %d %d\n", final_result.first[i].polygon_id,
                 final_result.second[j].polygon_id);
          std::vector<Polygon> lhs_p_vec = {final_result.first[i]};
          std::vector<Polygon> rhs_p_vec = {final_result.second[j]};
          save_polygons_to_csv(lhs_p_vec, "../lhs_p.csv");
          save_polygons_to_csv(rhs_p_vec, "../rhs_p.csv");
          exit(0);
        }
      }
    }
  }
}
#include "../../include/utils/pipeline.h"
#include "../../include/join_algos/inter_filter.h"
#include "../../include/join_algos/ri_join_algo.h"
#include "../../include/join_algos/serialize_polygon.h"
#include "../../include/mbr_algos/combined_sweep.h"
#include "../../include/mbr_algos/mbr_forward_scan.h"
#include "../../include/utils/compressing.h"
#include "../../include/utils/data_reader.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/hilbert.h"
#include "../../include/utils/minimum_bounding_rectangle.h"
#include <algorithm>
#include <chrono>
#include <set>
#include <vector>

void id_polygons(std::vector<Polygon> &polygons, bool lhs) {
  int id = 1;
  for (auto &it : polygons) {
    it.polygon_id = (lhs ? id : -id);
    id++;
  }
}

void get_preprocessed_polygons(std::vector<Polygon> &lhs,
                               std::vector<Polygon> &rhs,
                               const std::string &filename_lhs,
                               const std::string &filename_rhs,
                               Point &gridMinCorner, Point &gridMaxCorner,
                               int n) {
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

  // std::vector<Polygon> lhs_polygons_vec = {lhs_polygons[82950 - 1]};
  // std::vector<Polygon> rhs_polygons_vec = {rhs_polygons[75624 - 1]};

  // save_polygons_to_csv(lhs_polygons_vec, "../lhs_p.csv");
  // save_polygons_to_csv(rhs_polygons_vec, "../rhs_p.csv");

  // printf("%d\n", lhs_polygons[82950 - 1].intersects(rhs_polygons[75624 -
  // 1]));

  RasterGrid grid = RasterGrid(N, gridMinCorner, gridMaxCorner, 1e-7);

  // // print grid corners
  // printf("Grid corners: %lf %lf %lf %lf\n", gridMinCorner.x, gridMinCorner.y,
  //        gridMaxCorner.x, gridMaxCorner.y);

  // exit(0);
  std::vector<Polygon> lhs_polygons_copy = lhs_polygons;
  std::vector<Polygon> rhs_polygons_copy = rhs_polygons;

  // run mbr combined sweep ------------------------------------------------
  std::vector<std::pair<int, int>> final_result;
  forward_scan(lhs_polygons, rhs_polygons, final_result);

  std::set<int> lhs_ids;
  std::set<int> rhs_ids;

  for (auto &pair : final_result) {
    lhs_ids.insert(pair.first);
    rhs_ids.insert(pair.second);
  }

  std::vector<Polygon> lhs_filtered;
  std::vector<Polygon> rhs_filtered;

  for (auto &it : lhs_polygons) {
    if (lhs_ids.find(it.polygon_id) != lhs_ids.end()) {
      lhs_filtered.push_back(it);
    }
  }

  for (auto &it : rhs_polygons) {
    if (rhs_ids.find(it.polygon_id) != rhs_ids.end()) {
      rhs_filtered.push_back(it);
    }
  }

  std::vector<RasterPolygonInfo> lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info;

  std::set<int> null_cell_code_poly_idxs, error_poly_idxs;
  if (rasterize_polygons(grid, lhs_filtered, rhs_filtered,
                         lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info,
                         null_cell_code_poly_idxs, error_poly_idxs,
                         "bin_error.txt", false)) {
    printf("Error in rasterization\n");
  }

  printf("Ratio of wrong polygons: %f\n",
         (double)(null_cell_code_poly_idxs.size() + error_poly_idxs.size()) /
             (lhs_filtered.size() + rhs_filtered.size()));

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

  std::vector<std::pair<int, int>> final_filtered;
  // for each pair in final_result check if they are in null_cell_code_poly_idxs
  // or error_poly_idxs
  for (auto &r : final_result) {
    if (null_cell_code_poly_idxs.find(r.first) !=
            null_cell_code_poly_idxs.end() ||
        error_poly_idxs.find(r.first) != error_poly_idxs.end()) {
      continue;
    }
    if (null_cell_code_poly_idxs.find(r.second) !=
            null_cell_code_poly_idxs.end() ||
        error_poly_idxs.find(r.second) != error_poly_idxs.end()) {
      continue;
    }
    final_filtered.push_back(r);
  }

  // for each serialized polygon in lhs make a map from id to index
  std::map<int, int> lhs_id_to_idx;
  for (int i = 0; i < lhs_serialized_polygons.size(); i++) {
    lhs_id_to_idx[lhs_serialized_polygons[i].polygon_id] = i;
  }

  // for each serialized polygon in rhs make a map from id to index
  std::map<int, int> rhs_id_to_idx;
  for (int i = 0; i < rhs_serialized_polygons.size(); i++) {
    rhs_id_to_idx[rhs_serialized_polygons[i].polygon_id] = i;
  }

  std::set<std::pair<int, int>> result;
  std::set<std::pair<int, int>> indecisive;

  printf("Starting RI join\n");
  ri_join_algo(lhs_serialized_polygons, rhs_serialized_polygons, lhs_id_to_idx,
               rhs_id_to_idx, final_filtered, result, indecisive);

  // for each pair in result check if they really intersect
  int total_errors = 0;
  for (auto &r : result) {
    Polygon lhs_p = lhs_polygons_copy[r.first - 1];
    Polygon rhs_p = rhs_polygons_copy[-r.second - 1];

    if (!lhs_p.intersects(rhs_p)) {
      printf("ERROR ids: %d %d\n", r.first, r.second);
      std::vector<Polygon> lhs_p_vec = {lhs_p};
      std::vector<Polygon> rhs_p_vec = {rhs_p};

      save_polygons_to_csv(lhs_p_vec, "../lhs_p.csv");
      save_polygons_to_csv(rhs_p_vec, "../rhs_p.csv");
      break;
      total_errors++;
    }
  }

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

  printf("Indecisive ratio: %f\n",
         (double)indecisive.size() / final_result.size());

  printf("False positive: %f\n", (double)total_errors / final_filtered.size());
  // for each pair in final_filtered, if it does not belong in result, check if
  // it intersects
  int total_errors2 = 0;
  for (auto &r : final_filtered) {
    if (result.find(r) == result.end()) {
      Polygon lhs_p = lhs_polygons_copy[r.first - 1];
      Polygon rhs_p = rhs_polygons_copy[-r.second - 1];

      if (lhs_p.intersects(rhs_p)) {
        printf("ERROR ids: %d %d\n", r.first, r.second);
        std::vector<Polygon> lhs_p_vec = {lhs_p};
        std::vector<Polygon> rhs_p_vec = {rhs_p};

        save_polygons_to_csv(lhs_p_vec, "../lhs_p2.csv");
        save_polygons_to_csv(rhs_p_vec, "../rhs_p2.csv");
        total_errors2++;
      }
    }
  }
  printf("False negative: %f\n", (double)total_errors2 / final_filtered.size());
  printf("Total intersections: %d\n", result.size());
}

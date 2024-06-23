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

  id_polygons(lhs, true);

  rhs = read_data_find_MBR(filename_rhs, rhsMinCorner, rhsMaxCorner, n);

  id_polygons(rhs, false);

  gridMinCorner = {std::min(lhsMinCorner.x, rhsMinCorner.x),
                   std::min(lhsMinCorner.y, rhsMinCorner.y)};

  gridMaxCorner = {std::max(lhsMaxCorner.x, rhsMaxCorner.x),
                   std::max(lhsMaxCorner.y, rhsMaxCorner.y)};
}

void find_interesctions(
    unsigned int N, std::string lhs_f_name, std::string rhs_f_name,
    std::pair<std::vector<Polygon>, std::vector<Polygon>> &intersections,
    int n) {

  std::vector<Polygon> lhs_polygons;
  std::vector<Polygon> rhs_polygons;
  Point gridMinCorner, gridMaxCorner;

  get_preprocessed_polygons(lhs_polygons, rhs_polygons, lhs_f_name, rhs_f_name,
                            gridMinCorner, gridMaxCorner, n);

  int initial_total_polygons = lhs_polygons.size() + rhs_polygons.size();
  int initial_total_vertices = 0;
  for (auto &it : lhs_polygons) {
    initial_total_vertices += it.vertices.size();
  }
  for (auto &it : rhs_polygons) {
    initial_total_vertices += it.vertices.size();
  }

  RasterGrid grid = RasterGrid(N, gridMinCorner, gridMaxCorner, 1e-10);

  std::vector<Polygon> lhs_polygons_copy = lhs_polygons;
  std::vector<Polygon> rhs_polygons_copy = rhs_polygons;

  std::vector<std::pair<int, int>> final_result;

  // compute time take for forward scan to run

  auto start = std::chrono::high_resolution_clock::now();
  forward_scan(lhs_polygons, rhs_polygons, final_result);
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed_seconds = end - start;
  double mbr_filter_duration = elapsed_seconds.count();

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

  // time rasterization
  start = std::chrono::high_resolution_clock::now();
  if (rasterize_polygons(grid, lhs_filtered, rhs_filtered,
                         lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info,
                         null_cell_code_poly_idxs, error_poly_idxs,
                         "bin_error.txt", false)) {
    printf("Error in rasterization\n");
  }
  end = std::chrono::high_resolution_clock::now();

  elapsed_seconds = end - start;
  double rasterization_duration = elapsed_seconds.count();

  printf("Ratio of wrong polygons: %f\n",
         (double)(null_cell_code_poly_idxs.size() + error_poly_idxs.size()) /
             (lhs_filtered.size() + rhs_filtered.size()));

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

  // time serialization
  start = std::chrono::high_resolution_clock::now();
  serialize_polygons(lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info, hilbert_map,
                     lhs_serialized_polygons, rhs_serialized_polygons);
  end = std::chrono::high_resolution_clock::now();

  elapsed_seconds = end - start;
  double serialization_duration = elapsed_seconds.count();

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

  // time ri join
  start = std::chrono::high_resolution_clock::now();
  ri_join_algo(lhs_serialized_polygons, rhs_serialized_polygons, lhs_id_to_idx,
               rhs_id_to_idx, final_filtered, result, indecisive);
  end = std::chrono::high_resolution_clock::now();

  elapsed_seconds = end - start;
  double ri_join_duration = elapsed_seconds.count();

  int true_hits = result.size();
  double true_hits_ratio = (double)true_hits / final_filtered.size();
  int false_hits = final_filtered.size() - true_hits - indecisive.size();
  double false_hits_ratio = (double)false_hits / final_filtered.size();

  double indecisive_ratio = (double)indecisive.size() / final_filtered.size();

  // measure time for indecisive
  start = std::chrono::high_resolution_clock::now();
  for (auto &r : indecisive) {
    Polygon lhs_p = lhs_polygons_copy[r.first - 1];
    Polygon rhs_p = rhs_polygons_copy[-r.second - 1];

    if (lhs_p.intersects(rhs_p)) {
      result.insert(r);
    }
  }
  end = std::chrono::high_resolution_clock::now();

  elapsed_seconds = end - start;
  double indecisive_duration = elapsed_seconds.count();

  // print all metrics
  printf("Total polygons: %d\n", initial_total_polygons);
  printf("Average number of vertices: %d\n",
         initial_total_vertices / initial_total_polygons);
  printf("MBR filter duration: %.2f\n", mbr_filter_duration);
  printf("Rasterization duration: %.2f\n", rasterization_duration);
  printf("Serialization duration: %.2f\n", serialization_duration);
  printf("RI join duration: %.2f\n", ri_join_duration);
  printf("Refiment step duration: %.2f\n", indecisive_duration);
  printf("True hits percentage: %.2f\n", true_hits_ratio * 100);
  printf("False hits percentage: %.2f\n", false_hits_ratio * 100);
  printf("Indecisive percentage: %.2f\n", indecisive_ratio * 100);

  // // check if true hits are correct
  // int total_errors_true_hits = 0;
  // for (auto &r : result) {
  //   Polygon lhs_p = lhs_polygons_copy[r.first - 1];
  //   Polygon rhs_p = rhs_polygons_copy[-r.second - 1];

  //   if (!lhs_p.intersects(rhs_p)) {
  //     printf("Id %d and %d\n", r.first, r.second);
  //     printf("Error in true hits\n");
  //     total_errors_true_hits++;
  //   }
  // }

  // // check if false hits are correct
  // int total_errors_false_hits = 0;
  // for (auto &r : final_filtered) {
  //   if (result.find(r) != result.end()) {
  //     continue;
  //   }

  //   Polygon lhs_p = lhs_polygons_copy[r.first - 1];
  //   Polygon rhs_p = rhs_polygons_copy[-r.second - 1];

  //   if (lhs_p.intersects(rhs_p)) {
  //     printf("Id %d and %d\n", r.first, r.second);
  //     printf("Error in false hits\n");
  //     total_errors_false_hits++;
  //   }
  // }

  // // print total errors and with ratios
  // printf("Total errors in true hits: %d\n", total_errors_true_hits);
  // printf("Error in true hits ratio: %.2f\n",
  //        (double)total_errors_true_hits / true_hits * 100);
  // printf("Total errors in false hits: %d\n", total_errors_false_hits);
  // printf("Error in false hits ratio: %.2f\n",
  //        (double)total_errors_false_hits / false_hits * 100);
}
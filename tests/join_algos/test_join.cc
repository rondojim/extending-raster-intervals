#include "../../include/join_algos/rasterization.h"
#include "../../include/join_algos/ri_join_algo.h"
#include "../../include/join_algos/serialize_polygon.h"
#include "../../include/mbr_algos/combined_sweep.h"
#include "../../include/mbr_algos/mbr_forward_scan.h"
#include "../../include/utils/compressing.h"
#include "../../include/utils/data_reader.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/hilbert.h"
#include "../../include/utils/minimum_bounding_rectangle.h"
#include "../../include/utils/pipeline.h"
#include <algorithm>
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <thread>

void display_loading_bar(double progress) {
  int barWidth = 70;
  std::cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos)
      std::cout << "=";
    else if (i == pos)
      std::cout << ">";
    else
      std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
}

bool test_polygon_area(
    double orig_poly_area,
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        i_j_to_rcell_info) {

  double clipped_poly_area = 0.0;
  int full_cnt = 0;
  for (const auto &it : i_j_to_rcell_info) {
    RasterCellInfo rcell_info = it.second;
    BinaryCellCode c_code = rcell_info.cell_type;
    std::vector<std::vector<const Point *>> vertices_vectors =
        rcell_info.cell_polygons;
    double cur_area;

    if (strcmp(c_code.to_type(), "FULL") == 0) {
      cur_area = get_polygons_area(vertices_vectors);
    } else {
      cur_area = get_polygons_area(vertices_vectors);
    }

    clipped_poly_area += cur_area;
  }

  if (are_equal(orig_poly_area, clipped_poly_area, 1.0e-6)) {
    return true;
  }

  i_j_to_rcell_info.clear();
  std::cout << std::fixed << std::setprecision(10);
  std::cout << "orig_poly_are != clipped_poly_area: " << orig_poly_area
            << "!=" << clipped_poly_area << std::endl;

  return false;
}

int test_join(unsigned int N, std::string rhs_f_name, std::string lhs_f_name,
              std::string info_file_name = "",
              std::string err_poly_f_name = "") {

  std::vector<Polygon> lhs_polygons;
  std::vector<Polygon> rhs_polygons;
  Point gridMinCorner, gridMaxCorner;
  int n = -1;
  get_preprocessed_polygons(lhs_polygons, rhs_polygons, lhs_f_name, rhs_f_name,
                            gridMinCorner, gridMaxCorner, n);

  //   lhs_polygons[0].save_poly("lhs.txt");
  //   rhs_polygons[0].save_poly("rhs.txt");

  std::cout << "First dataset sz: : " << lhs_polygons.size() << std::endl;
  std::cout << "Second dataset sz: : " << rhs_polygons.size() << std::endl;

  // calculate area of each polygon to test after rasetrization
  // // if the sum of the clipped polygon area is equal
  // std::vector<double> lhs_polygons_area;
  // for (int i = 0; i < lhs_polygons.size(); i++) {
  //     lhs_polygons_area.push_back(polygon_area(lhs_polygons[i].vertices));
  // }
  // std::vector<double> rhs_polygons_area;
  // for (int i = 0; i < rhs_polygons.size(); i++) {
  //     rhs_polygons_area.push_back(polygon_area(rhs_polygons[i].vertices));
  // }

  // Construct Grid
  RasterGrid grid = RasterGrid(N, gridMinCorner, gridMaxCorner);
  std::cout << "Grid mbr: " << grid.max_corner.to_str()
            << grid.min_corner.to_str() << std::endl;

  std::set<int> null_cell_code_poly_idxs, error_poly_idxs;

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

  if (rasterize_polygons(grid, lhs_filtered, rhs_filtered,
                         lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info,
                         null_cell_code_poly_idxs, error_poly_idxs)) {
    printf("Error in rasterization\n");
  }

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

  printf("Ratio of wrong polygons: %f\n",
         (double)(null_cell_code_poly_idxs.size() + error_poly_idxs.size()) /
             (lhs_filtered.size() + rhs_filtered.size()));

  printf("Rasterization done\n");

  // load hilbert map

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

  printf("Indecisive ratio: %f\n",
         (double)indecisive.size() / final_result.size());

  printf("Ended RI join\n");

  std::set<std::pair<int, int>> join_ctype_result;
  std::set<std::pair<int, int>> join_ctype_indecisive;
  join_poly_cell_types(lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info,
                       join_ctype_result, join_ctype_indecisive);

  // for each join in join_ctype_result check if it exists in result
  for (const auto &it : join_ctype_result) {
    if (result.find(it) == result.end()) {
      std::cout << "Join in join_ctype_result not found in result\n";
      return 0;
    }
  }

  // for each join in join_ctype_indesicive check if it exists in indecisive
  for (const auto &it : join_ctype_indecisive) {
    if (indecisive.find(it) == indecisive.end()) {
      std::cout << "Join in join_ctype_indecisive not found in indecisive\n";
      return 0;
    }
  }

  // check if sizes are equal
  if (join_ctype_result.size() != result.size()) {
    std::cout << "join_ctype_result size != result size\n";
    return 0;
  }
  // check if sizes are equal
  if (join_ctype_indecisive.size() != indecisive.size()) {
    std::cout << "join_ctype_indecisive size != indecisive size\n";
    return 0;
  }
}

int main() {
  unsigned int N = 2;
  std::string rhs_f_name("../../dataset_files/archive/T1.csv");
  std::string lhs_f_name("../../dataset_files/archive/T2.csv");

  std::vector<RasterPolygonInfo> i_j_to_rpoly_info;
  int success = test_join(N, lhs_f_name, rhs_f_name, "", "error_polygons.csv");

  if (!success) {
    return success;
  }

  return 0;
}

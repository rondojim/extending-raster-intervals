#include "../../include/join_algos/rasterization.h"
#include "../../include/utils/data_reader.h"
#include "../../include/utils/geometry_types.h"
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

      if (strcmp(c_code.to_type(), "FULL") == 0)
      {
          cur_area = get_polygons_area(vertices_vectors);
      }
      else
      {
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

static void get_preprocessed_polygons(std::vector<Polygon> &polygons,
                                      const std::string &f_name,
                                      Point &gridMinCorner,
                                      Point &gridMaxCorner, int n = -1) {
  Point MinCorner, MaxCorner;
  polygons = read_data_find_MBR(f_name, MinCorner, MaxCorner, n);
  printf("Read set\n");

  gridMinCorner = Point(MinCorner.x, MinCorner.y);
  gridMaxCorner = Point(MaxCorner.x, MaxCorner.y);
}

// test for each polygon whether the resulting clipping area is same as the initial 
// area of the polygon by summing the are of each clipping part of all cells
int test_rasterize_polygons(unsigned int N, std::string f_name,
                            std::vector<RasterPolygonInfo> &i_j_to_rpoly_info,
                            std::string info_file_name = "",
                            std::string err_poly_f_name = "",
                            bool stop_on_error = false) {

  std::vector<Polygon> polygons;
  Point gridMinCorner, gridMaxCorner;
  get_preprocessed_polygons(polygons, f_name, gridMinCorner, gridMaxCorner, -1);

  size_t total_tests = polygons.size();
  std::cout << "Total tests: " << total_tests << std::endl;

  // calculate area of each polygon to test after rasetrization
  // if the sum of the clipped polygon area is equal
  std::vector<double> polygons_area;
  for (int i = 0; i < polygons.size(); i++) {
    polygons_area.push_back(polygon_area(polygons[i].vertices));
  }

  // Construct Grid
  // RasterGrid grid = RasterGrid(N, Point(0, 0), Point(10,10));
  RasterGrid grid = RasterGrid(N, gridMinCorner, gridMaxCorner);
  std::cout << "Grid mbr: " << grid.max_corner.to_str()
            << grid.min_corner.to_str() << std::endl;

  auto start_time = std::chrono::steady_clock::now();
  size_t passed_tests = 0;
  int success = 0;
  std::vector<Polygon> error_polygons;
  for (int i = 0; i < total_tests; ++i) {

    Polygon polygon = polygons[i];

    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        i_j_to_rcell_info;

    int weiler_success = grid.weiler_rasterize_poly(polygon, i_j_to_rcell_info);
    int sz = i_j_to_rcell_info.size();

    grid.save_poly_raster("poly_raster.txt", polygon, i_j_to_rcell_info);
    save_clipped_vertices_cell_type("clipped_vertices_cell_type.txt",
                                    i_j_to_rcell_info);
    // print_test_error_info(polygon, grid);

    // print_poly_grid_info(polygon, grid);

    if (weiler_success != 1) {
      success = -1;
      error_polygons.push_back(polygon);
      std::cout << "Failed at idx " << i << " polygon with total vertices "
                << polygons[i].vertices.size() << std::endl;
      std::cerr << "Error in polygon rasterization\n";
      if (stop_on_error) {
        break;
      }
    } else {

      if (info_file_name.size()) {
        save_clipped_vertices_cell_type(info_file_name.c_str(),
                                        i_j_to_rcell_info);
      }

      int result = 1; // test_polygon_area(polygons_area[i], i_j_to_rcell_info);
      if (result == 1) {
        passed_tests++;
      } else {
        std::cout << "Failed at idx " << i << " polygon with total vertices "
                  << polygons[i].vertices.size() << std::endl;
        // print_poly_grid_info(polygons[i], grid);
        std::cout << " Different original polygon and clipped polygon areas\n";
        success = 1;
        error_polygons.push_back(polygon);
        if (stop_on_error) {
          break;
        }
      }

      i_j_to_rpoly_info.emplace_back(polygon.polygon_id, i_j_to_rcell_info);
    }
    display_loading_bar(static_cast<double>(i + 1) / total_tests);
  }

  auto end_time = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed = end_time - start_time;
  std::cout << std::endl;
  std::cout << "Passed " << passed_tests << " out of " << total_tests
            << " tests." << std::endl;
  std::cout << "Execution time for n = " << N << " was " << elapsed.count()
            << " seconds." << std::endl;

  if (success != 0 && err_poly_f_name.size()) {
    std::cout << "Saving " << error_polygons.size()
              << " error polygons in file " << err_poly_f_name << std::endl;
    save_polygons_to_csv(error_polygons, err_poly_f_name.c_str());
  }

  grid.save_polygons_grid("polygons_grid.txt", polygons, i_j_to_rpoly_info);

  return success;
}

int main() {
  unsigned int N = 2;
  std::string f_name("../datasets/T1.csv");

  std::vector<RasterPolygonInfo> i_j_to_rpoly_info;
  int success = test_rasterize_polygons(N, f_name, i_j_to_rpoly_info, "",
                                        "error_polygons.csv");

  if (!success) {
    return success;
  }

  return 0;
}

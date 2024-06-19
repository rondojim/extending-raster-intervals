#include "../../include/join_algos/inter_filter.h"
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

int test_polygon_area(RasterGrid &grid, Polygon polygon, double orig_poly_area,
                      double &avg_full, bool fprint_unequal = true,
                      std::string grid_output_f_name = "", bool debug = false) {

    double grid_cell_area = grid.cell_area;
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        i_j_to_rcell_info;

    bool success = grid.weiler_rasterize_poly(polygon, i_j_to_rcell_info, debug);

    // if (debug)
    // {
        // save_vertices_vectors("inter_result.txt", j_to_vertices_vectors);
        if (grid_output_f_name.size())
        {
            // save_vertices_vectors("result.txt", i_j_to_vertices_vectors);
            save_rasterization(grid_output_f_name.c_str(), i_j_to_rcell_info);
            polygon.save_poly("polygon.txt");
        }
    //}

  if (!success) {
    return -1;
  }

  double clipped_poly_area = 0.0;
  int full_cnt = 0;
  for (const auto &it : i_j_to_rcell_info) {
    RasterCellInfo rcell_info = it.second;
    BinaryCellCode c_code = rcell_info.cell_type;
    std::vector<std::vector<const Point *>> vertices_vectors =
        rcell_info.cell_polygons;
    double cur_area;

        if (strcmp(c_code.to_type(true), "FULL") == 0)
        {
            cur_area = get_polygons_area(vertices_vectors);
            // if ( !are_equal(grid_cell_area, cur_area) ){
            //     std::cout << "["<< it.first.first << ", " << it.first.second << "]: " << "Should not be full size: grid_cell_area != cur_area -> " << grid_cell_area << " != " << cur_area << std::endl;
            //     return 0;
            // }
            // full_cnt++;
            // cur_area = grid_cell_area;
        }
        else
        {
            cur_area = get_polygons_area(vertices_vectors);
        }

    clipped_poly_area += cur_area;
  }

  avg_full = (double)full_cnt / (double)i_j_to_rcell_info.size();

  if (are_equal(orig_poly_area, clipped_poly_area, 1.0e-6)) {
    return 1;
  }

  if (fprint_unequal) {
    save_vertices(polygon.vertices, "error.txt");
    i_j_to_rcell_info.clear();
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "orig_poly_are != clipped_poly_area: " << orig_poly_area
              << "!=" << clipped_poly_area << std::endl;
  }

  return 0;
}

void print_test_error_info(Polygon polygon, RasterGrid grid)
{
    print_vec(polygon.vertices);
    std::cout << std::endl;

    Point min_corner = polygon.minCorner;
    Point max_corner = polygon.maxCorner;
    std::cout << "with mbr: min_corner: (" << min_corner.x << ", " << min_corner.y << "), (" << max_corner.x << ", " << max_corner.y << ")\n\n";
    std::cout << "X partition:\n";
    for (double x = grid.min_corner.x; is_less_or_equal(x, grid.max_corner.x); x += grid.step)
    {
        std::cout << x << "\t";
    }

  std::cout << std::endl;
  std::cout << "Y partition:\n";
  for (double y = grid.min_corner.y; is_less_or_equal(y, grid.max_corner.y);
       y += grid.step) {
    std::cout << y << "\t";
  }
  std::cout << std::endl;
  std::cout << "grid mbr:(g_xmin, g_ymin), (g_xmax, g_ymax), step: ("
            << grid.min_corner.x << ", " << grid.min_corner.y << "), ("
            << grid.max_corner.x << ", " << grid.max_corner.y << "), "
            << grid.step << "\n";
}

int main() {
  std::string grid_output_f_name("grid.txt");
  std::vector<const Point *> complex_poly_1 = {
      new Point(2.0, 9.0), new Point(7.0, 8.0), new Point(4.0, 7.5),
      new Point(7.0, 7.0), new Point(4.0, 4.0), new Point(7.0, 6.0),
      new Point(8.0, 4.0), new Point(8.5, 6.0), new Point(9.0, 4.0),
      new Point(2.0, 2.0), new Point(2.0, 9.0)};
  Point complex_poly_1_p1 = Point(5.0, 2.0);
  Point complex_poly_1_p2 = Point(5.0, 9.0);

  std::vector<const Point *> complex_poly_2 = {
      new Point(2.0, 9.0), new Point(7.0, 8.0), new Point(4.0, 7.5),
      new Point(7.0, 7.0), new Point(2.0, 7.0), new Point(2.0, 9.0)};
  Point complex_poly_2_p1 = Point(5.0, 7.0);
  Point complex_poly_2_p2 = Point(5.0, 9.0);

  std::vector<const Point *> snail_poly = {
      new Point(7.0, 10.0), new Point(0.0, 10.0), new Point(0.0, 0.0),
      new Point(8.0, 0.0),  new Point(8.0, 7.0),  new Point(2.0, 7.0),
      new Point(2.0, 4.0),  new Point(5.0, 4.0),  new Point(5.0, 5.0),
      new Point(3.0, 5.0),  new Point(3.0, 6.0),  new Point(7.0, 6.0),
      new Point(7.0, 1.0),  new Point(1.0, 1.0),  new Point(1.0, 9.0),
      new Point(7.0, 9.0),  new Point(7.0, 10.0)};

  // points for snail_poly
  /*
  "POLYGON((7 10, 0 10, 0 0, 8 0, 8 7, 2 7, 2 4, 5 4, 5 5, 3 5, 3 6, 7 6, 7 1, 1
  1, 1 9, 7 9, 7 10))",05OC
  */

  // reversed polygon of snail_poly
  /*
  "POLYGON((7 10, 7 9, 1 9, 1 1, 7 1, 7 6, 3 6, 3 5, 5 5, 5 4, 2 4, 2 7, 8 7, 8
  0, 0 0, 0 10, 7 10))",05OC
  */

  // points for complex_poly_1
  /*
  "POLYGON((2 9, 7 8, 4 7.5, 7 7, 4 4, 7 6, 8 4, 8.5 6, 9 4, 2 2, 2 9))",06OC
  */

  std::reverse(snail_poly.begin(), snail_poly.end());
  Point snail_poly_p1 = Point(4.0, 0.0);
  Point snail_poly_p2 = Point(4.0, 10.0);
  Point snail_poly_p3 = Point(4.0, 0.0);
  Point snail_poly_p4 = Point(4.0, 10.0);

  std::vector<const Point *> triangle = {
      new Point(5.0, 8.0), new Point(5.0, 9.0), new Point(7.0, 10.0),
      new Point(5.0, 8.0)};
  Point tr_p1 = Point(6, 8);
  Point tr_p2 = Point(6, 10);

  std::vector<const Point *> bug_poly = {
      new Point(6.375, 7.895833), new Point(5.9375, 7.822917),
      new Point(5.9375, 8.2125), new Point(6.375, 8.125),
      new Point(6.375, 7.895833)};
  Point bug_p1 = Point(5.9375, 8.125);
  Point bug_p2 = Point(6.375, 8.125);

    std::vector<std::vector<const Point *>> my_test_vectors = {
        // complex_poly_1,
        // complex_poly_2,
        snail_poly,
        // triangle,
        // bug_poly
        };

    std::vector<Polygon> polygons;
    Point gridMinCorner, gridMaxCorner;
    int cnt = 0;

    for (std::vector<const Point*> test_vertices : my_test_vectors) {
        Polygon polygon;
        polygon.vertices = test_vertices;
        polygon.findMBR();
        
        //std::cout << "mbr: " << polygon.maxCorner.to_str() << polygon.minCorner.to_str() << std::endl;
        polygons.push_back(polygon);
        if (cnt == 0)
        {
            gridMinCorner.x = polygon.minCorner.x, gridMinCorner.y = polygon.minCorner.y,
            gridMaxCorner.x = polygon.maxCorner.x, gridMaxCorner.y = polygon.maxCorner.y;
        }
        else
        {
            if (polygon.minCorner.x < gridMinCorner.x)
            {
                gridMinCorner.x = polygon.minCorner.x;
            }
            if (polygon.minCorner.y < gridMinCorner.y)
            {
                gridMinCorner.y = polygon.minCorner.y;
            }
            if (polygon.maxCorner.x > gridMaxCorner.x)
            {
                gridMaxCorner.x = polygon.maxCorner.x;
            }
            if (polygon.maxCorner.y > gridMaxCorner.y)
            {
                gridMaxCorner.y = polygon.maxCorner.y;
            }
        }
        cnt++;
    }

    // std::string filename = "../../datasets/T2.csv";
    std::string filename = "error_poly.txt";
    // std::string filename = "../../datasets/sample_simple_T1.csv";

        // std::string filename = "../../datasets/simple_T1.csv";
        // std::string filename = "../../datasets/sample_multi_polys_T1.csv";
        // std::string filename = "../../datasets/test.csv";

        
    // polygons = read_data_find_MBR(filename, gridMinCorner, gridMaxCorner, -1);

    size_t total_tests = polygons.size();
    std::cout << "Total tests: " << total_tests << std::endl;

    std::vector<double> polygons_area;
    for (int i = 0; i < polygons.size(); i++) {
        // std::cout << "sz: " << polygons[i].vertices.size() << std::endl;
        // print_vec(polygons[i].vertices);
        polygons_area.push_back(polygon_area(polygons[i].vertices));
    }

    bool error = false;
    unsigned int min_n = 2, max_n = 2;
    for (unsigned int n = min_n; n < max_n + 1; n += 5)
    {
        auto start_time = std::chrono::steady_clock::now();
        RasterGrid grid = RasterGrid(n, gridMinCorner, gridMaxCorner);
        std::cout << "Grid size: " << grid.n << " X " << grid.n << std::endl;
        std::cout << "Grid mbr: " << grid.max_corner.to_str() << grid.min_corner.to_str() << std::endl;
        size_t passed_tests = 0;

        double total_avgs = 0;
        for (int i = 0; i < total_tests; ++i)
        {
            grid.save_mbr_grid("mbr_grid.txt", polygons[i].minCorner, polygons[i].maxCorner);
            
            double cur_avg = 0;
            int result = test_polygon_area(grid, polygons[i], polygons_area[i], cur_avg, true, grid_output_f_name, false); // + std::to_string(i) + ".txt");
            total_avgs += cur_avg;

            if (result == 1)
            {
                passed_tests++;
            }
            else
            {
                std::cout << "Failed at idx " << i << " polygon with total vertices " << polygons[i].vertices.size() << ": ";
                polygons[i].save_poly_to_wkt("error_poly.txt");

                if (result == -1)
                {
                    std::cerr << "An unexpected error occured\n";
                }
                else
                {
                    //print_test_error_info(polygons[i], grid);
                    std::cout << " Different original polygon and clipped polygon areas\n";
                }
                // error = true;
                // break;
            }
            display_loading_bar(static_cast<double>(i + 1) / total_tests);
        }

    auto end_time = std::chrono::steady_clock::now();
    // Calculate duration
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << std::endl;
    std::cout << "Passed " << passed_tests << " out of " << total_tests
              << " tests." << std::endl;
    std::cout << "On average polygon include " << total_avgs / total_tests
              << " full cells\n";
    std::cout << "Execution time for n = " << n << " was " << elapsed.count()
              << " seconds." << std::endl;

    // if (error) {
    //   break;
    // }
  }

  return 0;
}

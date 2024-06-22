#include "../../include/join_algos/inter_filter.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/progressbar.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
const double EPSILON = 1e-24;

RasterGrid::RasterGrid(double n_, Point min_corner_, Point max_corner_)
    : n(std::pow(2.0, n_)), min_corner(min_corner_), max_corner(max_corner_) {

  double vert_size = max_corner.y - min_corner.y;
  double hor_size = max_corner.x - min_corner.x;
  if (vert_size >= hor_size) {
    step = vert_size / n;
    max_corner.x = min_corner.x + vert_size;
  } else {
    step = hor_size / n;
    max_corner.y = min_corner.y + hor_size;
  }
  std::cout << "Grid of size " << n << " X " << n << ", step: " << step
            << " initialized\n";
  cell_area = step * step;
}

double RasterGrid::min_above_or_equal_k(double xmin, double xmax, double k) {
  if (step <= 0) {
    std::cerr << "Step must be positive.\n";
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (k < xmin || k > xmax) {
    std::cerr << "Invalid input, given k is not in range [xmin,xmax]\n";
    return std::numeric_limits<double>::quiet_NaN();
  }

  double adjustment = fmod(k - xmin, step);
  if (adjustment != 0) {
    adjustment = step - adjustment;
  }

  double result = k + adjustment;

  // should be unneccessary
  double epsilon = 1e-8;
  if (result < xmin - epsilon || result > xmax + epsilon) {
    std::cout << std::fixed << std::setprecision(16);
    std::cerr << "Unexpected result not in [xmin, xmax] = [" << xmin << ", "
              << xmax << "], k = " << k << ", result = " << result << std::endl;
    return std::numeric_limits<double>::quiet_NaN();
  }

  return result;
}

double RasterGrid::max_below_or_equal_k(double xmin, double xmax, double k) {
  if (step <= 0) {
    std::cerr << "Step must be positive.\n";
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (k < xmin || k > xmax) {
    std::cerr << "Invalid input, given k is not in range [xmin,xmax]\n";
    std::cerr << "k: " << k << ", xmin: " << xmin << ", xmax: " << xmax
              << std::endl;
    return std::numeric_limits<double>::quiet_NaN();
  }

  double adjustment = fmod(k - xmin, step);
  double result = k - adjustment;

  // should be unneccessary
  double epsilon = 1e-8;
  if (result < xmin - epsilon || result > xmax + epsilon) {
    std::cerr << std::fixed << std::setprecision(20);
    std::cerr << "Unexpected result not in [xmin, xmax] = [" << xmin << ", "
              << xmax << "], k = " << k << ", result = " << result << std::endl;
    return std::numeric_limits<double>::quiet_NaN();
  }
  return result;
}

unsigned int RasterGrid::sequence_idx(double x, double xmin) {
  return std::round((x - xmin) / step);
}

int RasterGrid::k_belongs_in_sequence(double xmin, double xmax, double k,
                                      double epsilon) {
  if (step <= 0) {
    std::cerr << "Step must be positive.\n";
    return -1;
  }

  if (k < xmin || k > xmax) {
    std::cerr << "Invalid input, given k is not in range [xmin,xmax]\n";
    return -1;
  }

  // Calculate the modulus and check if it is close to 0 or step
  double remainder = std::fmod(k - xmin, step);
  if (std::fabs(remainder) < epsilon) {
    return 1;
  }

  return 0;
}

bool RasterGrid::set_segment_borders_types(
    const Point &p1, const Point &p2,
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        &i_j_to_rcell_info) {
  unsigned int parallelism = checkParallelism(p1, p2);
  int min_col_idx = -1, max_col_idx = -1, min_row_idx = -1, max_row_idx = -1;
  double xmin = min_corner.x, xmax = max_corner.y, ymin = min_corner.x,
         ymax = max_corner.y;

  if (!parallelism) {
    return true;
  }

  if (parallelism == 1) {
    // parallel to x axis
    int y_is_on_grid_row = k_belongs_in_sequence(xmin, xmax, p1.y);
    if (y_is_on_grid_row == -1) {
      return false;
    }
    if (y_is_on_grid_row) {
      int row_idx = sequence_idx(p1.y, ymin);
      min_row_idx, max_row_idx;
      if (row_idx == 0) {
        min_row_idx = 0;
        max_row_idx = 0;
      } else if (p1.y == ymax) {
        min_row_idx = row_idx - 1;
        max_row_idx = row_idx - 1;
      } else {
        min_row_idx = row_idx - 1;
        max_row_idx = row_idx;
      }

      double seg_min_x = p1.x;
      double seg_max_x = p2.x;
      if (p1.x > p2.x) {
        seg_min_x = p2.x;
        seg_max_x = p1.x;
      }

      double on_grid_min_x, on_grid_max_x;
      on_grid_min_x = max_below_or_equal_k(xmin, xmax, seg_min_x);
      on_grid_max_x = min_above_or_equal_k(xmin, xmax, seg_max_x);

      if (std::isnan(on_grid_min_x) || std::isnan(on_grid_max_x)) {
        std::cerr << "Error in set_segment_borders_types for parallel in x "
                     "axis: min_above_or_equal_k or "
                     "max_below_or_equal_k returned NaN\n";
        return false;
      }

      min_col_idx = sequence_idx(on_grid_min_x, xmin);
      max_col_idx = sequence_idx(on_grid_max_x, xmin);
      // if we have columns idx i, j, i < j then the cells'
      // position in map are i, i+1, ..., j-1
      max_col_idx--;

      if (are_equal(seg_min_x, on_grid_min_x) && min_col_idx > 0) {
        min_col_idx--;
      }

      if (are_equal(seg_max_x, on_grid_max_x) && on_grid_max_x < xmax) {
        max_col_idx++;
      }

    } else {
      return true;
    }
  } else if (parallelism == 2) {
    // parallel to y axis
    int y_is_on_grid_col = k_belongs_in_sequence(ymin, ymax, p1.x);
    if (y_is_on_grid_col == -1) {
      return false;
    }
    // std::cout << "y_is_on_grid_col: " << y_is_on_grid_col << std::endl;
    if (y_is_on_grid_col) {
      int col_idx = sequence_idx(p1.x, xmin);
      min_col_idx, max_col_idx;
      if (col_idx == 0) {
        min_col_idx = 0;
        max_col_idx = 0;
      } else if (p1.x == xmax) {
        min_col_idx = col_idx - 1;
        max_col_idx = col_idx - 1;
      } else {
        min_col_idx = col_idx - 1;
        max_col_idx = col_idx;
      }

      double seg_min_y = p1.y;
      double seg_max_y = p2.y;
      if (p1.y > p2.y) {
        seg_min_y = p2.y;
        seg_max_y = p1.y;
      }

      double on_grid_min_y, on_grid_max_y;
      on_grid_min_y = max_below_or_equal_k(ymin, ymax, seg_min_y);
      on_grid_max_y = min_above_or_equal_k(ymin, ymax, seg_max_y);

      if (std::isnan(on_grid_min_y) || std::isnan(on_grid_max_y)) {
        std::cerr << "Error in set_segment_borders_types for parallel in y "
                     "axis: min_above_or_equal_k or "
                     "max_below_or_equal_k returned NaN\n";
        return false;
      }

      min_row_idx = sequence_idx(on_grid_min_y, ymin);
      max_row_idx = sequence_idx(on_grid_max_y, ymin);

      // if we have rows min, max idxs i, j, i < j then the cells'
      // position in map are i, i+1, ..., j-1
      max_row_idx--;

      if (are_equal(seg_min_y, on_grid_min_y) && min_row_idx > 0) {
        min_row_idx--;
      }

      if (are_equal(seg_max_y, on_grid_max_y) && on_grid_max_y < ymax) {
        max_row_idx++;
      }
      // std::cout << "min_col_idx: " << min_col_idx << ", max_col_idx: " <<
      // max_col_idx << ", min_row_idx: " << min_row_idx << ", max_row_idx: " <<
      // max_row_idx << std::endl;

    } else {
      return true;
    }
  }

  // in case the segment is laying on a cell edge
  // having found the neighboor cells we set them the weak type

  BinaryCellCode weak_cell_code = BinaryCellCode(BinaryCellCode::WEAK_R);
  std::vector<std::vector<const Point *>> empty_vec;
  RasterCellInfo weak_rcell_info = RasterCellInfo(empty_vec, weak_cell_code);

  for (int j = min_row_idx; j <= max_row_idx; j++) {
    // double y1 = min_corner.y + j * step, y2 = min_corner.y + (j+1) * step;
    for (int i = min_col_idx; i <= max_col_idx; i++) {
      // double x1 = min_corner.x + i * step, x2 = min_corner.x + (i+1) * step;
      std::pair<int, int> j_i = {j, i};
      i_j_to_rcell_info[j_i] = weak_rcell_info;
    }
  }
  return true;
}

bool RasterGrid::set_polygon_borders_cell_types(
    std::vector<const Point *> &vertices,
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        &i_j_to_rcell_info) {
  for (int i = 0; i < vertices.size() - 1; i++) {

    int nxt_i = (i + 1) % vertices.size();
    const Point *p_i = vertices[i];
    const Point *p_nxt_i = vertices[nxt_i];
    if (!set_segment_borders_types(*p_i, *p_nxt_i, i_j_to_rcell_info)) {
      return false;
    }
  }

  return true;
}

bool RasterGrid::is_certain_full_cell(
    const std::vector<const Point *> &vertices, const Point &p1,
    const Point &p2, bool debug) {
  const Point p3(p2.x, p2.y + step);
  const Point p4(p1.x, p1.y + step);

  if (vertices.size() != 5) {
    return false;
  }

  std::array<const Point, 4> cell_corners = {p1, p2, p3, p4};
  // Check if all cell corners are in the vertices vector
  for (const Point &corner : cell_corners) {
    bool found = false;
    for (const Point *p : vertices) {
      if (*p == corner) {
        found = true;
        break; // Exit the inner loop early if the corner is found
      }
    }
    if (!found) {
      if (debug) {
        std::cout << "Not found corner: (" << corner.x << ", " << corner.y
                  << ")\n";
      }
      return false;
    }
  }

  return true;
}

BinaryCellCode
RasterGrid::encode(std::vector<std::vector<const Point *>> &vertices_vectors,
                   const Point &p1, const Point &p2, bool debug) {
  if (debug)
    std::cout << "in encode...\n";

  if (vertices_vectors.size() == 0) {
    std::cout << "aaaaaa\n";

    return BinaryCellCode(BinaryCellCode::NULL_CODE);
  }

  if (vertices_vectors.size() == 1 &&
      is_certain_full_cell(vertices_vectors[0], p1, p2)) {
    double cell_poly_area = get_polygons_area(vertices_vectors);
    if (!are_equal(cell_poly_area, cell_area, 1e-6)) {
      std::cout << "cell_poly_area: " << cell_poly_area
                << ", cell_area: " << cell_area << std::endl;
      is_certain_full_cell(vertices_vectors[0], p1, p2, true);
      std::cout << "cell_poly_area != cell_area \n";
      return BinaryCellCode(BinaryCellCode::NULL_CODE);
    }

    return BinaryCellCode(BinaryCellCode::FULL_R);
  }

  double cell_poly_area = get_polygons_area(vertices_vectors);

  if (are_equal(cell_poly_area, cell_area)) {

    return BinaryCellCode(BinaryCellCode::FULL_R);

  } else {
    if (cell_poly_area - cell_area > 1e-5) {

      std::cout << "cell_poly_area: " << cell_poly_area
                << ", cell_area: " << cell_area << std::endl;
      is_certain_full_cell(vertices_vectors[0], p1, p2, true);
      return BinaryCellCode(BinaryCellCode::NULL_CODE);
    }

    if (2.0 * cell_poly_area >= cell_area) {
      return BinaryCellCode(BinaryCellCode::STRONG_R);

    } else {
      if (debug) {
        std::cout << "cell_poly_area: " << cell_poly_area
                  << ", cell_area: " << cell_area << std::endl;
        std::cout << "WEAK\n";
      }
      return BinaryCellCode(BinaryCellCode::WEAK_R);
    }
  }
  return BinaryCellCode(BinaryCellCode::NULL_CODE);
}

bool are_equal(double a, double b, double epsilon) {
  return std::fabs(a - b) < epsilon;
}

bool is_greater_or_equal(double a, double b, double epsilon) {
  return (a > b) || (std::fabs(a - b) < epsilon);
}

bool is_less_or_equal(double a, double b, double epsilon) {
  return (a < b) || (std::fabs(a - b) < epsilon);
}

BinaryCellCode::BinaryCellCode(unsigned char val, bool type_R_)
    : value(val & 0b111), type_R(type_R_) {}

BinaryCellCode::BinaryCellCode(Code code, bool type_R_)
    : value(code & 0b111), type_R(type_R_) {}

// Get value
unsigned long BinaryCellCode::to_ulong() const { return value.to_ulong(); }

const char *BinaryCellCode::to_type(bool debug) const {
  unsigned long val = value.to_ulong();
  if (type_R) {
    if (val == FULL_R) {
      return "FULL";
    } else if (val == STRONG_R) {
      return "STRONG";
    } else if (val == WEAK_R) {
      return "WEAK";
    } else {
      return "NULL";
    }
  }

  if (val == FULL_S) {
    return "FULL";
  } else if (val == STRONG_S) {
    return "STRONG";
  } else if (val == WEAK_S) {
    return "WEAK";
  }

  return "NULL";
}

// Equality check
bool BinaryCellCode::equals(const BinaryCellCode &other) const {
  return (this->value == other.value) && (this->type_R == other.type_R);
}

// Overload equality operator
bool operator==(const BinaryCellCode &lhs, const BinaryCellCode &rhs) {
  return lhs.equals(rhs);
}

// Overload inequality operator
bool operator!=(const BinaryCellCode &lhs, const BinaryCellCode &rhs) {
  return !(lhs == rhs);
}

bool rasterize_polygons(RasterGrid &grid, std::vector<Polygon> &lhs_polygons,
                        std::vector<Polygon> &rhs_polygons,
                        std::vector<RasterPolygonInfo> &lhs_i_j_to_rpoly_info,
                        std::vector<RasterPolygonInfo> &rhs_i_j_to_rpoly_info,
                        std::set<int> &null_cell_code_poly_idxs,
                        std::set<int> &error_poly_idxs,
                        std::string err_poly_f_name, bool debug) {

  progressbar bar(lhs_polygons.size() + rhs_polygons.size());

  bool error = false;
  for (Polygon &polygon : lhs_polygons) {
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        i_j_to_rcell_info;
    int success = grid.weiler_rasterize_poly(polygon, i_j_to_rcell_info, debug);
    if (!success) {
      if (err_poly_f_name.size()) {
        polygon.save_poly(err_poly_f_name.c_str());
      }
      error_poly_idxs.insert(polygon.polygon_id);
      error = true;
    } else if (success == 2) {
      null_cell_code_poly_idxs.insert(polygon.polygon_id);
      error = true;
    } else {
      lhs_i_j_to_rpoly_info.emplace_back(polygon.polygon_id, i_j_to_rcell_info);
    }
    bar.update();
  }

  printf("done with lhs\n");
  for (Polygon &polygon : rhs_polygons) {
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        i_j_to_rcell_info;

    int success = grid.weiler_rasterize_poly(polygon, i_j_to_rcell_info, debug);
    if (!success) {
      if (err_poly_f_name.size()) {
        polygon.save_poly(err_poly_f_name.c_str());
      }
      error_poly_idxs.insert(polygon.polygon_id);
      error = true;
    } else if (success == 2) {
      null_cell_code_poly_idxs.insert(polygon.polygon_id);
      error = true;
    } else {
      rhs_i_j_to_rpoly_info.emplace_back(polygon.polygon_id, i_j_to_rcell_info);
    }
    bar.update();
  }

  return error;
}

int get_double_sigh(double x, double epsilon) {
  if (std::fabs(x - 0.0) < epsilon) {
    return 0;
  } else if (x < 0.0) {
    return -1;
  }
  return 1;
}

// 1 true hit
// 2 indecisive
int join_cell_types(BinaryCellCode l_cell_type, BinaryCellCode r_cell_type) {
  if (r_cell_type == BinaryCellCode::FULL_S ||
      l_cell_type == BinaryCellCode::FULL_R ||
      (r_cell_type == BinaryCellCode::STRONG_S &&
       l_cell_type == BinaryCellCode::STRONG_R)) {
    return 1;
  }

  return 2;
}

void join_poly_cell_types(std::vector<RasterPolygonInfo> &lhs_i_j_to_rpoly_info,
                          std::vector<RasterPolygonInfo> &rhs_i_j_to_rpoly_info,
                          std::set<std::pair<int, int>> &result,
                          std::set<std::pair<int, int>> &indecisive) {

  for (RasterPolygonInfo &l_rpoly_info : lhs_i_j_to_rpoly_info) {
    int l_idx = l_rpoly_info.idx;
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        &l_rcells_info = l_rpoly_info.i_j_to_rcell_info;

    for (RasterPolygonInfo &r_rpoly_info : rhs_i_j_to_rpoly_info) {
      int r_idx = r_rpoly_info.idx;
      std::pair<int, int> l_r_idx = {l_idx, r_idx};
      int inter = 0;

      std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
          &r_rcells_info = r_rpoly_info.i_j_to_rcell_info;

      for (const auto &l_cell : l_rcells_info) {
        std::pair<unsigned int, unsigned int> l_cell_j_i = l_cell.first;
        RasterCellInfo l_cell_info = l_cell.second;

        auto r_cell_it = r_rcells_info.find(l_cell_j_i);
        if (r_cell_it != r_rcells_info.end()) {

          RasterCellInfo r_cell_info = r_cell_it->second;
          inter = join_cell_types(l_cell_info.cell_type, r_cell_info.cell_type);
          if (inter == 1) {
            std::cout << "interesct in i,j: " << l_cell_j_i.first << ", "
                      << l_cell_j_i.second << std::endl;
            break;
          }
        }
      }

      if (inter == 1) {
        result.insert(l_r_idx);
      } else if (inter == 2) {
        indecisive.insert(l_r_idx);
      }
    }
  }
}
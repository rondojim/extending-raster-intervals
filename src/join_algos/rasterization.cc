#include "../../include/join_algos/rasterization.h"
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

RasterGrid::RasterGrid(double n_, Point min_corner_, Point max_corner_,
                       double prec_epsilon_)
    : n(std::pow(2.0, n_)), min_corner(min_corner_), max_corner(max_corner_),
      prec_epsilon(prec_epsilon_) {

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
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (k < xmin || k > xmax) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double adjustment = fmod(k - xmin, step);
  if (adjustment != 0) {
    adjustment = step - adjustment;
  }

  double result = k + adjustment;

  // should be unneccessary
  if (result < xmin - prec_epsilon || result > xmax + prec_epsilon) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return result;
}

double RasterGrid::max_below_or_equal_k(double xmin, double xmax, double k) {
  if (step <= 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (k < xmin || k > xmax) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double adjustment = fmod(k - xmin, step);
  double result = k - adjustment;

  // should be unneccessary
  if (result < xmin - prec_epsilon || result > xmax + prec_epsilon) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return result;
}

unsigned int RasterGrid::sequence_idx(double x, double xmin) {
  return std::round((x - xmin) / step);
}

int RasterGrid::k_belongs_in_sequence(double xmin, double xmax, double k) {
  if (step <= 0) {
    return -1;
  }

  if (k < xmin || k > xmax) {
    return -1;
  }

  // Calculate the modulus and check if it is close to 0 or step
  double remainder = std::fmod(k - xmin, step);
  if (std::fabs(remainder) < prec_epsilon) {
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
  double xmin = min_corner.x, xmax = max_corner.x, ymin = min_corner.y,
         ymax = max_corner.y;
  if (!parallelism) {
    return true;
  }
  if (parallelism == 1) {

    // parallel to x axis
    int y_is_on_grid_row = k_belongs_in_sequence(ymin, ymax, p1.y);
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
        return false;
      }

      min_col_idx = sequence_idx(on_grid_min_x, xmin);
      max_col_idx = sequence_idx(on_grid_max_x, xmin);
      // if we have columns idx i, j, i < j then the cells'
      // position in map are i, i+1, ..., j-1
      max_col_idx--;

      if (are_equal(seg_min_x, on_grid_min_x, prec_epsilon) &&
          min_col_idx > 0) {
        // the max x of the line which overalaps a grid column is on a corner
        // include the lower column too
        min_col_idx--;
      }

      if (are_equal(seg_max_x, on_grid_max_x, prec_epsilon) &&
          on_grid_max_x < xmax) {
        // the min x of the line which overalaps a grid column is on a corner
        // include the higher column too
        max_col_idx++;
      }

    } else {
      return true;
    }
  } else if (parallelism == 2) {
    // parallel to y axis
    int y_is_on_grid_col = k_belongs_in_sequence(xmin, xmax, p1.x);
    if (y_is_on_grid_col == -1) {
      return false;
    }
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
        return false;
      }

      min_row_idx = sequence_idx(on_grid_min_y, ymin);
      max_row_idx = sequence_idx(on_grid_max_y, ymin);

      // if we have rows min, max idxs i, j, i < j then the cells'
      // position in map are i, i+1, ..., j-1
      max_row_idx--;

      if (are_equal(seg_min_y, on_grid_min_y, prec_epsilon) &&
          min_row_idx > 0) {
        // the min y of the line which overalaps a grid row is on a corner
        // include the lower row too
        min_row_idx--;
      }

      if (are_equal(seg_max_y, on_grid_max_y, prec_epsilon) &&
          on_grid_max_y < ymax) {
        // the max y of the line which overalaps a grid row is on a corner
        // include the higher row too
        max_row_idx++;
      }
    } else {
      return true;
    }
  }

  // in case the segment is laying on a cell edge
  // having found the neighboor cells we set them the weak type

  BinaryCellCode weak_cell_code = BinaryCellCode(BinaryCellCode::WEAK_R);
  std::vector<std::vector<const Point *>> empty_vec;
  RasterCellInfo weak_rcell_info = RasterCellInfo(empty_vec, weak_cell_code);

  for (int i = min_row_idx; i <= max_row_idx; i++) {
    for (int j = min_col_idx; j <= max_col_idx; j++) {
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
    const Point &p2) {
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
      if (corner.equal_points(*p, prec_epsilon)) {
        found = true;
        break; // Exit the inner loop early if the corner is found
      }
    }
    if (!found) {
      return false;
    }
  }

  return true;
}

BinaryCellCode
RasterGrid::encode(std::vector<std::vector<const Point *>> &vertices_vectors,
                   const Point &p1, const Point &p2) {

  if (vertices_vectors.size() == 0) {
    return BinaryCellCode(BinaryCellCode::NULL_CODE);
  }

  if (vertices_vectors.size() == 1 &&
      is_certain_full_cell(vertices_vectors[0], p1, p2)) {
    double cell_poly_area = get_polygons_area(vertices_vectors);
    if (!are_equal(cell_poly_area, cell_area, prec_epsilon)) {
      return BinaryCellCode(BinaryCellCode::NULL_CODE);
    }

    return BinaryCellCode(BinaryCellCode::FULL_R);
  }

  double cell_poly_area = get_polygons_area(vertices_vectors);

  if (are_equal(cell_poly_area, cell_area, prec_epsilon)) {

    return BinaryCellCode(BinaryCellCode::FULL_R);

  } else {
    if (cell_poly_area - cell_area > clipped_cell_area_tol) {
      return BinaryCellCode(BinaryCellCode::NULL_CODE);
    }
    if (2.0 * cell_poly_area >= cell_area) {
      return BinaryCellCode(BinaryCellCode::STRONG_R);
    } else {
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

const char *BinaryCellCode::to_type() const {
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
                        std::set<int> &error_poly_idxs, int algo) {

  progressbar bar(lhs_polygons.size() + rhs_polygons.size());

  bool error = false;
  for (Polygon &polygon : lhs_polygons) {
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        i_j_to_rcell_info;
    

    int success; 
    if (algo == 0) {
      success = grid.weiler_rasterize_poly(polygon, i_j_to_rcell_info);
    }
    else {
      success = grid.hodgman_rasterize_poly(polygon, i_j_to_rcell_info);
    }
    
    if (!success) {
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

    int success; 
    if (algo == 0) {
      success = grid.weiler_rasterize_poly(polygon, i_j_to_rcell_info);
    }
    else {
      success = grid.hodgman_rasterize_poly(polygon, i_j_to_rcell_info);
    }
    
    if (!success) {
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

int RasterGrid::get_double_sigh(double x) {
  if (std::fabs(x - 0.0) < prec_epsilon) {
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

void print_poly_grid_info(Polygon polygon, RasterGrid grid) {
  // print_vec(polygon.vertices);
  std::cout << std::endl;

  Point min_corner = polygon.minCorner;
  Point max_corner = polygon.maxCorner;
  std::cout << "with mbr: min_corner: (" << min_corner.x << ", " << min_corner.y
            << "), (" << max_corner.x << ", " << max_corner.y << ")\n\n";
  std::cout << "X partition:\n";
  for (double x = grid.min_corner.x;
       is_less_or_equal(x, grid.max_corner.x, 1e-6); x += grid.step) {
    std::cout << x << "\t";
  }

  std::cout << std::endl;
  std::cout << "Y partition:\n";
  for (double y = grid.min_corner.y;
       is_less_or_equal(y, grid.max_corner.y, 1e-6); y += grid.step) {
    std::cout << y << "\t";
  }
  std::cout << std::endl;
  std::cout << "grid mbr:(g_xmin, g_ymin), (g_xmax, g_ymax), step: ("
            << grid.min_corner.x << ", " << grid.min_corner.y << "), ("
            << grid.max_corner.x << ", " << grid.max_corner.y << "), "
            << grid.step << "\n";
}
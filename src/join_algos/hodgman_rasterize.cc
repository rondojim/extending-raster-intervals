#include "../../include/join_algos/inter_filter.h"
#include "../../include/utils/geometry_types.h"
#include <cmath>
#include <iostream>
#include <map>

void RasterGrid::hodgman_clip_segment(
    int i_pos_sigh, int nxt_i_pos_sigh, const Point *p_inter, const Point *p_nxt_i,
    std::vector<const Point *> &left_vert_clipped_points,
    std::vector<const Point *> &right_vert_clipped_points, bool debug) {
  // Case 1: For the right side: When both points are inside
  if (i_pos_sigh < 0 && nxt_i_pos_sigh < 0) {
    if (debug) {
      std::cout << "case 1.1\n";
    }
    // Only second point is added
    right_vert_clipped_points.push_back(p_nxt_i);
  }
  // Case 2: For the right side:
  // First point is outside or on line
  // Second is inside
  else if (i_pos_sigh >= 0 && nxt_i_pos_sigh < 0) {
    if (debug) {
      std::cout << "case 1.2\n";
    }
    // Point of intersection with edge and the second point is added
    right_vert_clipped_points.push_back(p_inter);
    right_vert_clipped_points.push_back(p_nxt_i);
  }
  // Case 3: For the right side:
  // First point is inside
  // Second point is outside or on line
  else if (i_pos_sigh < 0 && nxt_i_pos_sigh >= 0) {
    if (debug) {
      std::cout << "case 1.3\n";
    }
    // Only point of intersection with edge is added
    right_vert_clipped_points.push_back(p_inter);
  }
  // Case 4: For the right side: When both points are outside, no points are
  // added

  // Case 1: For the left side: Both points are inside
  if (i_pos_sigh > 0 && nxt_i_pos_sigh > 0) {
    if (debug) {
      std::cout << "case 2.1\n";
    }
    left_vert_clipped_points.push_back(p_nxt_i);
  }
  // Case 2: First vertex is outside while second one is inside
  else if (i_pos_sigh <= 0  && nxt_i_pos_sigh > 0) {
    if (debug) {
      std::cout << "case 2.2\n";
    }
    // For the left clipped point of intersection with edge and the second point
    // is added
    left_vert_clipped_points.push_back(p_inter);
    left_vert_clipped_points.push_back(p_nxt_i);
  }
  // Case 3: First vertex is inside while second one is outside
  else if (i_pos_sigh > 0 && nxt_i_pos_sigh <= 0) {
    if (debug) {
      std::cout << "case 2.3\n";
    }
    // For the fully clipped only point of intersection with edge is added
    left_vert_clipped_points.push_back(p_inter);
  }
}  

bool RasterGrid::hodgman_clip(const std::vector<const Point *> &vertices,
      std::vector<const Point *> &left_vert_clipped_points,
      std::vector<const Point *> &right_vert_clipped_points,
      const Point &p1, const Point &p2, bool debug){

  int poly_size = vertices.size();
  for (int i = 0; i < poly_size - 1; i++) {
    // i and nxt_i form a line in polygon
    int nxt_i = (i + 1) % poly_size;
    const Point *p_i = vertices[i];
    const Point *p_nxt_i = vertices[nxt_i];

    // Calculating position of first, second point w.r.t. clipper line
    double i_pos = cross_product(p1, p2, *p_i);
    double nxt_i_pos = cross_product(p1, p2, *p_nxt_i);

    int i_pos_sigh = get_double_sigh(i_pos);
    int nxt_i_pos_sigh = get_double_sigh(nxt_i_pos);


    Point *p_inter = nullptr;
    // if the points are on different sides wrt to the p1-p2 line or at least
    // one of them is on the line
    if (i_pos_sigh * nxt_i_pos_sigh <= 0) {
      p_inter = p_intersect(p1, p2, *p_i, *p_nxt_i, prec_epsilon);
      if (p_inter == nullptr) {
        std::cerr << "Nan intersection\n";
        return false;
      }
    }

    hodgman_clip_segment(i_pos_sigh, nxt_i_pos_sigh, p_inter, p_nxt_i, left_vert_clipped_points,
                 right_vert_clipped_points, debug);
  }

  // Append first point as last to include the edge between them
  if (left_vert_clipped_points.size()) {
    left_vert_clipped_points.push_back(left_vert_clipped_points[0]);
  }

  if (right_vert_clipped_points.size()) {
    right_vert_clipped_points.push_back(right_vert_clipped_points[0]);
  }

  return true;
}

bool RasterGrid::hodgman_rasterize_poly(Polygon &polygon,
                             std::map<std::pair<unsigned int, unsigned int>,
                                      RasterCellInfo> &i_j_to_rcell_info,
                             bool debug){

    if (!set_polygon_borders_cell_types(polygon.vertices, i_j_to_rcell_info)) {
    std::cerr << "Error in set_polygon_borders_cell_types\n";
    return false;
  }

  // Find the corners of the grid 
  // that enclose the polygon mbr 
  double xmin =
      max_below_or_equal_k(min_corner.x, max_corner.x, polygon.minCorner.x);
  double ymin =
      max_below_or_equal_k(min_corner.y, max_corner.y, polygon.minCorner.y);
  double xmax =
      min_above_or_equal_k(min_corner.x, max_corner.x, polygon.maxCorner.x);
  double ymax =
      min_above_or_equal_k(min_corner.y, max_corner.y, polygon.maxCorner.y);

  if (std::isnan(xmin) || std::isnan(ymin) || std::isnan(xmax) ||
      std::isnan(ymax)) {
    std::cerr << "Error in hodgman_rasterize_poly: min_above_or_equal_k or "
                 "max_below_or_equal_k returned NaN\n";
    return false;
  }

  // Find the rows, columns indexes of previous coordinates on the mbr
  int xmin_idx = sequence_idx(xmin, min_corner.x);
  int xmax_idx = sequence_idx(xmax, min_corner.x);

  int ymin_idx = sequence_idx(ymin, min_corner.y);
  int ymax_idx = sequence_idx(ymax, min_corner.y);

  std::vector<const Point*> fully_vert_clipped_vertices, semi_vert_clipped_vertices, 
      cur_vert_clipped_vertices, fully_hori_clipped_vertices, 
      semi_hori_clipped_vertices, cur_hori_clipped_vertices; 

  std::map<std::pair<int, int>, std::vector<std::vector<const Point *>>>
      i_j_to_vertical_clipped, i_j_to_hori_clipped;

  cur_vert_clipped_vertices = polygon.vertices;
  semi_vert_clipped_vertices = polygon.vertices; 

  double x_j, y_i;
  for (int j = xmin_idx; j < xmax_idx; j++) {
      x_j = min_corner.x + (j + 1) * step;
      fully_vert_clipped_vertices.clear();

      if (j == xmax_idx - 1) {
          // in last column the semi_vert_clipped_vertices, ie
          // the remaining clipped vertices are fully clipped 
          cur_hori_clipped_vertices = semi_vert_clipped_vertices;
      }
      else {
          semi_vert_clipped_vertices.clear();
          Point p1_vert(x_j, ymin);
          Point p2_vert(x_j, ymax);
          if (!hodgman_clip(cur_vert_clipped_vertices,
              fully_vert_clipped_vertices, semi_vert_clipped_vertices, p1_vert,
              p2_vert, debug)) 
          {
              return false;
          }

          cur_vert_clipped_vertices = semi_vert_clipped_vertices;
          cur_hori_clipped_vertices = fully_vert_clipped_vertices;
      }

      semi_hori_clipped_vertices = cur_hori_clipped_vertices;
      
      // FOR TESTING
      if (!semi_hori_clipped_vertices.empty()) {
          std::pair<int, int> i_j = {j, j};
          i_j_to_vertical_clipped[i_j].push_back(cur_hori_clipped_vertices);
      }


      for (int i = ymax_idx - 1; i >= ymin_idx; i--) {
          fully_hori_clipped_vertices.clear();
          std::pair<int, int> i_j = {j, i};
           y_i = min_corner.y + i * step;
          Point p1_hori(x_j - step, y_i);
          Point p2_hori(x_j, y_i);

          if (i == ymin_idx) {
            // in the lowest row the semi_hori_clipped_vertices, ie
            // the remaining clipped vertices are fully clipped (being up to the last row)

              if (semi_hori_clipped_vertices.size()) {
                  std::vector<std::vector<const Point*>> semi_hori_clipped_vertices_vec = {semi_hori_clipped_vertices};
                  BinaryCellCode cell_code =
                      encode(semi_hori_clipped_vertices_vec, p1_hori, p2_hori);
                  if (cell_code.value == BinaryCellCode::NULL_CODE) {
                      save_vertices_vectors_seg(semi_hori_clipped_vertices_vec, p1_hori, p2_hori,
                                      "wrong_hodgman_clip.txt");
                      std::cerr << "NUll Cell Code\n";
                      std::cout << "Failed in [" << i << ", " << j << "]\n";
                      std::cout << "cell polygons: ";
                      for (std::vector<const Point *> &vec : semi_hori_clipped_vertices_vec) {
                        print_vec(vec);
                      }
                      return false;
                  }
                  i_j_to_rcell_info[i_j] = RasterCellInfo(semi_hori_clipped_vertices_vec, cell_code);
                  i_j_to_hori_clipped[i_j].push_back(semi_hori_clipped_vertices);
              }
          }
          else {
              semi_hori_clipped_vertices.clear();
              if (!hodgman_clip(cur_hori_clipped_vertices,
                  fully_hori_clipped_vertices, semi_hori_clipped_vertices, p1_hori,
                  p2_hori, debug))
              {
                return false;
              }

              cur_hori_clipped_vertices = semi_hori_clipped_vertices;

              if (fully_hori_clipped_vertices.size()) {
                  std::vector<std::vector<const Point*>> fully_hori_clipped_vertices_vec = {fully_hori_clipped_vertices};
                  BinaryCellCode cell_code =
                      encode(fully_hori_clipped_vertices_vec, p1_hori, p2_hori);
                  if (cell_code.value == BinaryCellCode::NULL_CODE) {
                    save_vertices_vectors_seg(fully_hori_clipped_vertices_vec, p1_hori, p2_hori,
                                    "wrong_hodgman_clip.txt");
                    std::cerr << "NUll Cell Code\n";
                    std::cout << "Failed in [" << i << ", " << j << "]\n";
                    std::cout << "cell polygons: ";
                    for (std::vector<const Point *> &vec :
                        fully_hori_clipped_vertices_vec) {
                      print_vec(vec);
                    }
                    return false;
                  }
                  i_j_to_rcell_info[i_j] =
                      RasterCellInfo(fully_hori_clipped_vertices_vec, cell_code);
                      i_j_to_hori_clipped[i_j].push_back(fully_hori_clipped_vertices);
              }

          }

      }
    }
    save_clipped_vertices_vectors("hodgman_column_rasterization.txt", i_j_to_vertical_clipped, polygon.minCorner,polygon.maxCorner);
    save_clipped_vertices_vectors("hodgman_row_rasterization.txt", i_j_to_hori_clipped, polygon.minCorner,polygon.maxCorner);


    return true;
}

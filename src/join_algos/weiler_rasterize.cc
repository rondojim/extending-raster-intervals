#include "../../include/join_algos/inter_filter.h"
#include "../../include/utils/geometry_types.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>

bool RasterGrid::weiler_clip_vertices_vectors(
    const std::vector<std::vector<const Point *>> &vertices_vectors,
    std::vector<std::vector<const Point *>> &left_vertices_vectors,
    std::vector<std::vector<const Point *>> &right_vertices_vectors,
    const Point &p1, const Point &p2) {

  for (const std::vector<const Point *> &vertices : vertices_vectors) {

    if (!weiler_clip(vertices, left_vertices_vectors, right_vertices_vectors,
                     p1, p2, 100000)) {
      return false;
    }
  }
  return true;
}

void RasterGrid::sort_inter_points(std::vector<size_t> &indexes,
                                   std::vector<InterPointInfo> &inter_pts_info,
                                   bool sort_by_y) {
  // Sort the indexes vector based on the boolean flag
  std::sort(indexes.begin(), indexes.end(),
            [&inter_pts_info, sort_by_y](size_t i1, size_t i2) {
              if (sort_by_y) {
                return inter_pts_info[i1].inter_point->y <
                       inter_pts_info[i2].inter_point->y;
              } else {
                return inter_pts_info[i1].inter_point->x <
                       inter_pts_info[i2].inter_point->x;
              }
            });

  // Sort the actual inter_pts_info vector based on the boolean flag
  std::sort(
      inter_pts_info.begin(), inter_pts_info.end(),
      [sort_by_y](const InterPointInfo &info1, const InterPointInfo &info2) {
        if (sort_by_y) {
          return info1.inter_point->y < info2.inter_point->y;
        } else {
          return info1.inter_point->x < info2.inter_point->x;
        }
      });
}

int RasterGrid::weiler_scan_polygon(const std::vector<const Point *> &vertices,
                                    const Point &p1, const Point &p2,
                                    std::vector<PointInfo> &vertices_with_info,
                                    std::vector<InterPointInfo> &inter_pts_info) {
  int poly_size = vertices.size();
  int total_inter_pts = 0;
  bool point_on_left = false, point_on_right = false;

  for (int i = 0; i < poly_size - 1; i++) {

    int nxt_i = (i + 1) % poly_size;
    const Point *p_i = vertices[i];
    const Point *p_nxt_i = vertices[nxt_i];

    // Calculating position of first, second point w.r.t. clipper line
    double i_pos = cross_product(p1, p2, *p_i);
    double nxt_i_pos = cross_product(p1, p2, *p_nxt_i);
    int i_pos_sigh = get_double_sigh(i_pos);
    int nxt_i_pos_sigh = get_double_sigh(nxt_i_pos);

    if (i_pos_sigh < 0) {
      // we have at least one point on the right of the segment p1,p2
      point_on_right = true;
    } else if (i_pos_sigh > 0) {
      // we have at least one point on left of the segment p1,p2
      point_on_left = true;
    }

    const Point *p_inter = nullptr;
    Direction dir = Direction::NONE;

    // indicates the index on the intersection point
    // in case we have interesection and we are entering/exiting
    // the polygon vertices
    int inter_p_idx = inter_pts_info.size();

    bool are_collinear = (i_pos_sigh == 0 && nxt_i_pos_sigh == 0);

    // indicates the index on the polygon point
    // in case we are entering/exiting the polygon
    int nxt_point_idx = vertices_with_info.size();
    if (are_collinear) {
      p_inter = p_i;
      if (i == poly_size - 2) {
        nxt_point_idx = 0;
      }
      InterPointInfo inter_p_info = InterPointInfo(p_inter, dir, nxt_point_idx);
      inter_pts_info.push_back(inter_p_info);
    } else if (i_pos_sigh * nxt_i_pos_sigh <= 0) {
      // we are sure the current i, nxt_i segment intersects the p1,p2

      // left or right direction wrt to sef p1, p2 based on orientation
      if ((i_pos_sigh == 0 && nxt_i_pos_sigh > 0) ||
          (i_pos_sigh < 0 && nxt_i_pos_sigh == 0) ||
          (i_pos_sigh < 0 && nxt_i_pos_sigh > 0)) {
        dir = Direction::LEFT;
      } else {
        dir = Direction::RIGHT;
      }

      // we are sure there is interesction and
      // we find the interesction point
      p_inter = p_intersect(p1, p2, *p_i, *p_nxt_i, prec_epsilon);

      if (nxt_i_pos_sigh != 0) {
        if (p_inter == nullptr) {
          return -1;
        }

        if (i != poly_size - 2) {
          if (i_pos_sigh != 0) {
            nxt_point_idx++;
          }
        } else {
          nxt_point_idx = 0; // no other vertice will be added, the interestion
                             // tail is the first point inserted
        }

        InterPointInfo inter_p_info =
            InterPointInfo(p_inter, dir, nxt_point_idx);
        inter_pts_info.push_back(inter_p_info);
      } else if (i == poly_size - 2) {
        inter_p_idx = 0;
      }
    }

    // Do not insert if current point is collinear to the line
    // Because it will be inserted multiple times in the final result
    if (i_pos_sigh != 0) {
      PointInfo p_info = PointInfo(p_i, dir, inter_p_idx);
      vertices_with_info.push_back(p_info);
    }
  }

  if (!(point_on_left && point_on_right)) {
    if (!point_on_left) {
      return 1; // all points on the right of the line
    }
    return 2; // all points on the left of the line
  }
  return 0;
}

bool RasterGrid::weiler_clip(
    const std::vector<const Point *> &vertices,
    std::vector<std::vector<const Point *>> &left_vertices_vectors,
    std::vector<std::vector<const Point *>> &right_vertices_vectors,
    const Point &p1, const Point &p2, unsigned int max_iters) {

  std::vector<PointInfo> vertices_with_info;
  std::vector<InterPointInfo> inter_pts_info;

  int points_side = weiler_scan_polygon(vertices, p1, p2, vertices_with_info,
                                        inter_pts_info);
  if (points_side == -1) {
    return false; // for error
  }

  if (points_side == 1) {
    // if all points on the right
    right_vertices_vectors.push_back(vertices);
  } else if (points_side == 2) {
    // if all points on the left
    left_vertices_vectors.push_back(vertices);
  } else {
    // indexes corresponding to inter_pts_info initi position
    std::vector<size_t> indexes(inter_pts_info.size());
    std::iota(indexes.begin(), indexes.end(), 0);

    // if we do vertical clipping sort by x
    // else by y
    bool vertical_clip = (p1.x == p2.x);
    // we sort indexes to, in order to have a map
    // from the initial pos to the pos after sorting
    sort_inter_points(indexes, inter_pts_info, vertical_clip);

    // old_to_new role is the maping of the old inter_pts_info
    // before sorting to the new after sorting
    std::vector<size_t> old_to_new(inter_pts_info.size());
    for (size_t new_index = 0; new_index < indexes.size(); ++new_index) {
      old_to_new[indexes[new_index]] = new_index;
    }

    // a concave polygon after clipping might result
    // in multiple polygons wrt to a line
    // eg the shape of letter M closed in the bottom
    // clipped wrt the a segment parallel on the x axis
    // interesecting it somewhere higher than the middle y
    std::vector<const Point *> cur_polygon;
    int total_iters = 0;
    for (InterPointInfo &inter_p_info : inter_pts_info) {

      // visit each entering point only once and a NONE dir indicated
      // collinear point the segment which we do not claim as
      if (!inter_p_info.visited && inter_p_info.dir != Direction::NONE) {
        // entering point belong the clipped polygon
        cur_polygon.push_back(inter_p_info.inter_point);
        inter_p_info.visited = true;
        const Point *cur_p_inter = nullptr;
        int exit_inter_pts_idx = -1;
        int enter_vertices_idx = inter_p_info.nxt_point_idx;

        Direction enter_dir = inter_p_info.dir;
        Direction exit_dir;
        std::vector<std::vector<const Point *>> *polygon_vector_to_insert =
            nullptr;

        // if the enter dir is right the inter points
        // should be traversed from bottom to top
        // and the opposite if enter dir is left
        int inter_pts_step;

        if (inter_p_info.dir == Direction::RIGHT) {
          inter_pts_step = 1;
          exit_dir = Direction::LEFT;
          polygon_vector_to_insert = &right_vertices_vectors;
        } else if (inter_p_info.dir == Direction::LEFT) {
          inter_pts_step = -1;
          exit_dir = Direction::RIGHT;
          polygon_vector_to_insert = &left_vertices_vectors;
        } 

        while (cur_p_inter != inter_p_info.inter_point &&
               total_iters++ < max_iters) {

          exit_inter_pts_idx = walk_on_polygon_vertices(
              cur_polygon, vertices_with_info, exit_dir, enter_vertices_idx,
              max_iters);
          if (exit_inter_pts_idx == -1) {
            std::cerr << "Error in walk_on_polygon_vertices\n";
            return false;
          }

          int exit_inter_pts_maped_idx = old_to_new[exit_inter_pts_idx];

          enter_vertices_idx = walk_on_inter_vertices(
              cur_polygon, inter_pts_info, enter_dir, exit_inter_pts_maped_idx,
              inter_pts_step, &cur_p_inter, max_iters);

          if (enter_vertices_idx == -1) {
            std::cerr << "Error in walk_on_inter_vertices\n";
            return false;
          }
        }

        polygon_vector_to_insert->push_back(cur_polygon);
        cur_polygon.clear();

        if (total_iters > max_iters) {
          std::cout << "Error in weiler: max_iters == total_iters\n";
          return false;
        }
      }
    }
  }

  return true;
}

int RasterGrid::walk_on_polygon_vertices(
    std::vector<const Point *> &polygon_to_fill,
    std::vector<PointInfo> &vertices_with_info, Direction exit_dir,
    int start_idx, unsigned int max_iters) {

  Direction cur_dir = vertices_with_info[start_idx].dir;
  int total_iters = 0;
  int cur_idx = start_idx;

  while (cur_dir != exit_dir && total_iters++ != max_iters) {
    polygon_to_fill.push_back(vertices_with_info[cur_idx].point);
    cur_idx = (cur_idx + 1) % vertices_with_info.size();
    cur_dir = vertices_with_info[cur_idx].dir;
  }

  if (total_iters > max_iters) {
    // for case of infinite loop
    return -1;
  }

  polygon_to_fill.push_back(vertices_with_info[cur_idx].point);
  return vertices_with_info[cur_idx].inter_point_idx;
}

int RasterGrid::walk_on_inter_vertices(
    std::vector<const Point *> &polygon_to_fill,
    std::vector<InterPointInfo> &inter_pts_info, Direction exit_dir,
    int start_idx, int step, const Point **exit_point, unsigned int max_iters) {

  Direction cur_dir = inter_pts_info[start_idx].dir;
  int total_iters = 0;
  int cur_idx = start_idx;

  while (cur_dir != exit_dir && total_iters++ != max_iters) {
    polygon_to_fill.push_back(inter_pts_info[cur_idx].inter_point);
    cur_idx = (cur_idx + step) < 0 ? inter_pts_info.size() - 1
                                   : (cur_idx + step) % inter_pts_info.size();
    cur_dir = inter_pts_info[cur_idx].dir;
  }

  if (total_iters > max_iters) {
    // for case of infinite loop
    return -1;
  }

  inter_pts_info[cur_idx].visited = true;
  *exit_point = inter_pts_info[cur_idx].inter_point;
  polygon_to_fill.push_back(*exit_point);

  return inter_pts_info[cur_idx].nxt_point_idx;
}

int RasterGrid::weiler_rasterize_poly(
    Polygon &polygon,
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        &i_j_to_rcell_info) {

  if (!set_polygon_borders_cell_types(polygon.vertices, i_j_to_rcell_info)) {
    return 0;
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
    return 0;
  }

  // Find the rows, columns indexes of previous coordinates on the mbr
  int xmin_idx = sequence_idx(xmin, min_corner.x);
  int xmax_idx = sequence_idx(xmax, min_corner.x);

  int ymin_idx = sequence_idx(ymin, min_corner.y);
  int ymax_idx = sequence_idx(ymax, min_corner.y);

  std::vector<std::vector<const Point *>> vertices_vectors;
  vertices_vectors.push_back(polygon.vertices);

  std::vector<std::vector<const Point *>> fully_vert_clipped_vertices,
      semi_vert_clipped_vertices, cur_vert_clipped_vertices,
      fully_hori_clipped_vertices, semi_hori_clipped_vertices,
      cur_hori_clipped_vertices;

  std::map<std::pair<int, int>, std::vector<std::vector<const Point *>>>
      i_j_to_vertical_clipped, i_j_to_hori_clipped;

  cur_vert_clipped_vertices = vertices_vectors;
  semi_vert_clipped_vertices = vertices_vectors;
  double x_j, y_i;

  for (int j = xmin_idx; j < xmax_idx; j++) {
    x_j = min_corner.x + (j + 1) * step;
    fully_vert_clipped_vertices.clear();

    if (j == xmax_idx - 1) {
      // in the last column the semi_vert_clipped_vertices, ie
      // the remaining clipped vertices are fully clipped
      cur_hori_clipped_vertices = semi_vert_clipped_vertices;
    } else {
      semi_vert_clipped_vertices.clear();
      Point p1_vert(x_j, ymin);
      Point p2_vert(x_j, ymax);

      bool success = weiler_clip_vertices_vectors(
          cur_vert_clipped_vertices, fully_vert_clipped_vertices,
          semi_vert_clipped_vertices, p1_vert, p2_vert);
      if (!success) {
        return 0;
      }

      cur_vert_clipped_vertices = semi_vert_clipped_vertices;
      cur_hori_clipped_vertices = fully_vert_clipped_vertices;
    }

    semi_hori_clipped_vertices = cur_hori_clipped_vertices;

    // FOR TESTING
    if (!semi_hori_clipped_vertices.empty()) {
      std::vector<std::vector<const Point *>> cur_vec;

      for (const auto &polygon : semi_hori_clipped_vertices) {
        std::vector<const Point *> new_polygon;
        for (const auto *point : polygon) {
          if (point) {
            // Create a new Point instance for deep copy
            Point *new_point = new Point(point->x, point->y);
            new_polygon.push_back(new_point);
          }
        }
        cur_vec.push_back(new_polygon);
      }

      std::pair<int, int> i_j = {j, j};
      i_j_to_vertical_clipped[i_j] = cur_vec;
    }

    for (int i = ymax_idx - 1; i >= ymin_idx; i--) {
      fully_hori_clipped_vertices.clear();
      std::pair<int, int> i_j = {j, i};

      y_i = min_corner.y + i * step;
      Point p1_hori(x_j - step, y_i);
      Point p2_hori(x_j, y_i);

      if (i == ymin_idx) {
        // in the lowest row the semi_hori_clipped_vertices, ie
        // the remaining clipped vertices are fully clipped (being up to the
        // last row)

        if (semi_hori_clipped_vertices.size()) {
          BinaryCellCode cell_code =
              encode(semi_hori_clipped_vertices, p1_hori, p2_hori);
          if (cell_code.value == BinaryCellCode::NULL_CODE) {
            return 2;
          }
          i_j_to_rcell_info[i_j] =
              RasterCellInfo(semi_hori_clipped_vertices, cell_code);
          i_j_to_hori_clipped[i_j] = semi_hori_clipped_vertices;
        }
      } else {
        semi_hori_clipped_vertices.clear();

        int success = weiler_clip_vertices_vectors(
            cur_hori_clipped_vertices, fully_hori_clipped_vertices,
            semi_hori_clipped_vertices, p1_hori, p2_hori);

        cur_hori_clipped_vertices = semi_hori_clipped_vertices;
        if (fully_hori_clipped_vertices.size()) {
          BinaryCellCode cell_code =
              encode(fully_hori_clipped_vertices, p1_hori, p2_hori);
          if (cell_code.value == BinaryCellCode::NULL_CODE) {
            return 2;
          }
          i_j_to_rcell_info[i_j] =
              RasterCellInfo(fully_hori_clipped_vertices, cell_code);
          i_j_to_hori_clipped[i_j] = fully_hori_clipped_vertices;
        }
      }
    }
  }


  return 1;
}
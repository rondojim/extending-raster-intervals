#include "../../include/join_algos/inter_filter.h"
#include "../../include/utils/geometry_types.h"
#include <cmath>
#include <iostream>
#include <map>
#include <algorithm>
#include <numeric>


bool RasterGrid::weiler_clip_vertices_vectors(
    const std::vector<std::vector<const Point *>> &vertices_vectors,
    std::vector<std::vector<const Point *>> &left_vertices_vectors,
    std::vector<std::vector<const Point *>> &right_vertices_vectors,
    const Point &p1, const Point &p2, bool debug) {

  for (const std::vector<const Point *> &vertices : vertices_vectors) {

    if (!weiler_clip(vertices, left_vertices_vectors, right_vertices_vectors,
                     p1, p2, 100000, debug)) {
      return false;
    }
  }
  return true;
}

void sort_inter_points(std::vector<size_t> &indexes,
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
                                    std::vector<InterPointInfo> &inter_pts_info,
                                    bool debug) {
  int poly_size = vertices.size();
  int total_inter_pts = 0;
  bool point_on_left = false, point_on_right = false;

  for (int i = 0; i < poly_size - 1; i++) {

    int nxt_i = (i + 1) % poly_size;
    const Point *p_i = vertices[i];
    const Point *p_nxt_i = vertices[nxt_i];

    if (debug) {
      std::cout << p_i->to_str() << " - " << p_nxt_i->to_str() << std::endl;
    }

    // Calculating position of first, second point w.r.t. clipper line
    double i_pos = cross_product(p1, p2, *p_i);
    double nxt_i_pos = cross_product(p1, p2, *p_nxt_i);
    int i_pos_sigh = get_double_sigh(i_pos);
    int nxt_i_pos_sigh = get_double_sigh(nxt_i_pos);

    if (i_pos_sigh < 0) {
      point_on_right = true;
    } else if (i_pos_sigh > 0) {
      point_on_left = true;
    }

    if (debug) {
      std::cout << "i_pos: " << i_pos << ", nxt_i_pos: " << nxt_i_pos
                << std::endl;
    }

    const Point *p_inter = nullptr;
    Direction dir = Direction::NONE;
    int inter_p_idx = inter_pts_info.size();

    bool are_collinear = (i_pos_sigh == 0 && nxt_i_pos_sigh == 0);

    int nxt_point_idx = vertices_with_info.size();
    if (are_collinear) {
      p_inter = p_i;
      if (i == poly_size - 2) {
        nxt_point_idx = 0;
      }
      InterPointInfo inter_p_info = InterPointInfo(p_inter, dir, nxt_point_idx);
      inter_pts_info.push_back(inter_p_info);
    } else if (i_pos_sigh * nxt_i_pos_sigh <= 0) {
      if ((i_pos_sigh == 0 && nxt_i_pos_sigh > 0) ||
          (i_pos_sigh < 0 && nxt_i_pos_sigh == 0) ||
          (i_pos_sigh < 0 && nxt_i_pos_sigh > 0)) {
        dir = Direction::LEFT;
        if (debug) {
          std::cout << "left\n";
        }
      } else {
        if (debug) {
          std::cout << "right\n";
        }
        dir = Direction::RIGHT;
      }

      p_inter = p_intersect(p1, p2, *p_i, *p_nxt_i);

      if (nxt_i_pos_sigh != 0) {
        if (p_inter == nullptr) {
          std::cerr << "Nan intersection\n";
          return -1;
        }

        if (i != poly_size - 2) {
          if (i_pos_sigh != 0) {
            nxt_point_idx++;
          }
        } else {
          nxt_point_idx = 0; // no other vertice will be added, the interestion
                             // tail is the first point inseted
        }

        if (debug) {
          std::cout << "nxt_point_idx: " << nxt_point_idx << std::endl;
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
    const Point &p1, const Point &p2, unsigned int max_iters, bool debug) {
  std::vector<PointInfo> vertices_with_info;
  std::vector<InterPointInfo> inter_pts_info;

  int points_side = weiler_scan_polygon(vertices, p1, p2, vertices_with_info,
                                        inter_pts_info, debug);
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
    std::vector<size_t> indexes(inter_pts_info.size());
    std::iota(indexes.begin(), indexes.end(), 0);

    // if we do vertical clipping sort by x
    // else by y
    bool vertical_clip = (p1.x == p2.x);
    sort_inter_points(indexes, inter_pts_info, vertical_clip);

    // old_to_new role is the maping of the old inter_pts_info
    // before sorting to the new after sorting
    std::vector<size_t> old_to_new(inter_pts_info.size());
    for (size_t new_index = 0; new_index < indexes.size(); ++new_index) {
      old_to_new[indexes[new_index]] = new_index;
    }

    if (debug) {
      int inter_vert_size = inter_pts_info.size();

      std::cout << "\nvertices_with_info:\n";
      for (int i = 0; i < vertices_with_info.size(); i++) {
        std::cout << i << ") ";
        vertices_with_info[i].print_(false);
        int inter_point_idx = vertices_with_info[i].inter_point_idx;
        if (vertices_with_info[i].dir == 1) {
          std::cout << " -> ";
        } else if (vertices_with_info[i].dir == 2) {
          std::cout << " <- ";
        } else {
          std::cout << " - ";
        }
        if (vertices_with_info[i].dir) {
          if (vertices_with_info[i].inter_point_idx > inter_vert_size - 1) {
            std::cout << "\nvert_info inter point index > inter_vert_size\n";
            return false;
          }
          std::cout << "maped_idx: " << old_to_new[inter_point_idx] << " "
                    << inter_pts_info[old_to_new[inter_point_idx]]
                           .inter_point->to_str();
        }
        std::cout << std::endl;
      }
    }
    if (debug) {
      int vert_info_size = vertices_with_info.size();

      std::cout << "\ninter_pts_info:\n";
      for (int i = 0; i < inter_pts_info.size(); i++) {
        std::cout << i << ") ";
        inter_pts_info[i].print_(false);
        int nxt_point_idx = inter_pts_info[i].nxt_point_idx;
        if (inter_pts_info[i].dir == 1) {
          std::cout << " -> ";
        } else if (inter_pts_info[i].dir == 2) {
          std::cout << " <- ";
        } else {
          std::cout << " - ";
        }
        if (inter_pts_info[i].dir != Direction::NONE) {
          if (nxt_point_idx > vert_info_size - 1) {
            std::cout << "\ninter_pts_info nxt point index > nxt_point_idx\n";
            return false;
          }
          std::cout << "nxt_point_idx: " << nxt_point_idx << ": "
                    << vertices_with_info[nxt_point_idx].point->to_str();
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }

    std::vector<const Point *> cur_polygon;
    int total_iters = 0;
    for (InterPointInfo &inter_p_info : inter_pts_info) {

      if (!inter_p_info.visited && inter_p_info.dir != Direction::NONE) {
        if (debug) {
          std::cout << "START: ( " << inter_p_info.inter_point->x << ", "
                    << inter_p_info.inter_point->y << ")\n";
        }

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
        } else {
          std::cerr << "Error: There is inter_p_info.dir with NULL direction\n";
        }

        // std::cout << "step: " << inter_pts_step << std::endl;
        while (cur_p_inter != inter_p_info.inter_point &&
               total_iters++ < max_iters) {
          if (debug) {
            std::cout << "enter_vertice: "
                      << vertices_with_info[enter_vertices_idx].point->to_str()
                      << std::endl;
          }

          exit_inter_pts_idx = walk_on_polygon_vertices(
              cur_polygon, vertices_with_info, exit_dir, enter_vertices_idx,
              max_iters, debug);
          if (exit_inter_pts_idx == -1) {
            std::cerr << "Error in walk_on_polygon_vertices\n";
            return false;
          }

          int exit_inter_pts_maped_idx = old_to_new[exit_inter_pts_idx];
          if (debug) {
            std::cout << "before map idx: " << exit_inter_pts_idx << std::endl;
            std::cout << "exit_vertice: "
                      << inter_pts_info[exit_inter_pts_maped_idx]
                             .inter_point->to_str()
                      << std::endl;
          }

          enter_vertices_idx = walk_on_inter_vertices(
              cur_polygon, inter_pts_info, enter_dir, exit_inter_pts_maped_idx,
              inter_pts_step, &cur_p_inter, max_iters, debug);

          if (debug) {
            std::cout << "cur_p_inter: " << cur_p_inter->to_str() << std::endl;
          }

          if (enter_vertices_idx == -1) {
            std::cerr << "Error in walk_on_inter_vertices\n";
            return false;
          }
        }

        if (debug) {
          std::cout << "POLYGON push back: \n";
          print_vec(cur_polygon);
          std::cout << std::endl;
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

  if (debug) {
    std::cout << "left:\n";
    for (std::vector<const Point *> vec : left_vertices_vectors) {
      print_vec(vec);
    }

    std::cout << "right:\n";
    for (std::vector<const Point *> vec : right_vertices_vectors) {
      print_vec(vec);
    }
  }

  return true;
}

int RasterGrid::walk_on_polygon_vertices(
    std::vector<const Point *> &polygon_to_fill,
    std::vector<PointInfo> &vertices_with_info, Direction exit_dir,
    int start_idx, unsigned int max_iters, bool debug) {
  Direction cur_dir = vertices_with_info[start_idx].dir;
  int total_iters = 0;
  int cur_idx = start_idx;
  if (debug) {
    std::cout << "walk_on_polygon: exit dir " << exit_dir << "\n";
  }

  while (cur_dir != exit_dir && total_iters++ != max_iters) {
    polygon_to_fill.push_back(vertices_with_info[cur_idx].point);
    if (debug) {
      std::cout << "(" << vertices_with_info[cur_idx].point->x << ", "
                << vertices_with_info[cur_idx].point->y
                << "), dir: " << vertices_with_info[cur_idx].dir << "\n";
    }
    cur_idx = (cur_idx + 1) % vertices_with_info.size();
    cur_dir = vertices_with_info[cur_idx].dir;
  }

  if (total_iters > max_iters) {
    std::cerr << "Error in walk_on_polygon: max_iters == total_iters\n";
    return -1;
  }

  polygon_to_fill.push_back(vertices_with_info[cur_idx].point);
  if (debug) {
    std::cout << "(" << vertices_with_info[cur_idx].point->x << ", "
              << vertices_with_info[cur_idx].point->y
              << "), dir: " << vertices_with_info[cur_idx].dir << "\n";
  }

  return vertices_with_info[cur_idx].inter_point_idx;
}

int RasterGrid::walk_on_inter_vertices(
    std::vector<const Point *> &polygon_to_fill,
    std::vector<InterPointInfo> &inter_pts_info, Direction exit_dir,
    int start_idx, int step, const Point **exit_point, unsigned int max_iters,
    bool debug) {
  Direction cur_dir = inter_pts_info[start_idx].dir;
  int total_iters = 0;
  int cur_idx = start_idx;
  if (debug) {
    std::cout << "walk_on_intersections: exit dir " << exit_dir << "\n";
  }

  while (cur_dir != exit_dir && total_iters++ != max_iters) {
    polygon_to_fill.push_back(inter_pts_info[cur_idx].inter_point);
    if (debug) {
      std::cout << "(" << inter_pts_info[cur_idx].inter_point->x << ", "
                << inter_pts_info[cur_idx].inter_point->y
                << "), dir: " << inter_pts_info[cur_idx].dir << "\n";
    }
    cur_idx =  (cur_idx + step) < 0 ? inter_pts_info.size() - 1 : (cur_idx + step) % inter_pts_info.size(); 
    cur_dir = inter_pts_info[cur_idx].dir;
  }

  if (total_iters > max_iters) {
    std::cerr << "Error in walk_on_intersections: max_iters == total_iters\n";
    return -1;
  }

  inter_pts_info[cur_idx].visited = true;
  *exit_point = inter_pts_info[cur_idx].inter_point;
  polygon_to_fill.push_back(*exit_point);
  if (debug) {
    std::cout << "exit_point: (" << inter_pts_info[cur_idx].inter_point->x
              << ", " << inter_pts_info[cur_idx].inter_point->y
              << "), dir: " << inter_pts_info[cur_idx].dir << " visited\n";
  }

  return inter_pts_info[cur_idx].nxt_point_idx;
}


int RasterGrid::weiler_rasterize_poly(
    Polygon &polygon,
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        &i_j_to_rcell_info,
    bool debug) 
{

  if (!set_polygon_borders_cell_types(polygon.vertices, i_j_to_rcell_info)) {
    std::cerr << "Error in set_polygon_borders_cell_types\n";
    return 0;
  }

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
    std::cerr << "Error in weiler_rasterize_poly: min_above_or_equal_k or "
                 "max_below_or_equal_k returned NaN\n";
    return 0;
  }

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

  std::map<unsigned int, std::vector<std::vector<const Point *>>>
      j_to_vertices_vectors;

  cur_vert_clipped_vertices = vertices_vectors;
  semi_vert_clipped_vertices =
      vertices_vectors; 
  double x_j, y_i;

  for (int j = xmin_idx; j < xmax_idx; j++) {
    x_j = min_corner.x + (j + 1) * step;
    fully_vert_clipped_vertices.clear();

    if (j == xmax_idx - 1) {
      cur_hori_clipped_vertices = semi_vert_clipped_vertices;
    } else {
      semi_vert_clipped_vertices.clear();
      Point p1_vert(x_j, ymin);
      Point p2_vert(x_j, ymax);

      bool success = weiler_clip_vertices_vectors(
          cur_vert_clipped_vertices, fully_vert_clipped_vertices,
          semi_vert_clipped_vertices, p1_vert, p2_vert, debug);
      if (!success) {
        std::cout << "Failed in j " << j << "\n";
        return 0;
      }

      if (!success) {
        std::cout << "Failed in j: " << j << "\n";
        return 0;
      }

      cur_vert_clipped_vertices = semi_vert_clipped_vertices;
      cur_hori_clipped_vertices = fully_vert_clipped_vertices;
    }

    semi_hori_clipped_vertices = cur_hori_clipped_vertices;

    // FOR TESTING
    j_to_vertices_vectors[j] = semi_hori_clipped_vertices;

    for (int i = ymax_idx - 1; i >= ymin_idx; i--) {
      fully_hori_clipped_vertices.clear();
      std::pair<int, int> i_j = {j, i};

      y_i = min_corner.y + i * step;
      Point p1_hori(x_j - step, y_i);
      Point p2_hori(x_j, y_i);

      if (i == ymin_idx) {
        if (semi_hori_clipped_vertices.size()) {
          BinaryCellCode cell_code =
              encode(semi_hori_clipped_vertices, p1_hori, p2_hori);
          if (cell_code.value == BinaryCellCode::NULL_CODE) {
            save_vertices_vectors_seg(semi_hori_clipped_vertices, p1_hori, p2_hori,
                            "wrong_weiler_clip.txt");
            std::cerr << "NUll Cell Code\n";
            std::cout << "Failed in [" << i << ", " << j << "]\n";
            std::cout << "cell polygons: ";
            for (std::vector<const Point *> &vec : semi_hori_clipped_vertices) {
              print_vec(vec);
            }
            return 2;
          }
          i_j_to_rcell_info[i_j] =
              RasterCellInfo(semi_hori_clipped_vertices, cell_code);
        }
      } else {
        semi_hori_clipped_vertices.clear();

        int success = weiler_clip_vertices_vectors(
            cur_hori_clipped_vertices, fully_hori_clipped_vertices,
            semi_hori_clipped_vertices, p1_hori, p2_hori, debug);
        if (!success) {
          save_vertices_vectors("col_error.txt", j_to_vertices_vectors[j]);
          std::cout << "Failed in [" << i << ", " << j << "]\n";
          return 0;
        }

        cur_hori_clipped_vertices = semi_hori_clipped_vertices;
        if (fully_hori_clipped_vertices.size()) {
          BinaryCellCode cell_code =
              encode(fully_hori_clipped_vertices, p1_hori, p2_hori);
          if (cell_code.value == BinaryCellCode::NULL_CODE) {
            save_vertices_vectors_seg(fully_hori_clipped_vertices, p1_hori, p2_hori,
                            "wrong_weiler_clip.txt");
            std::cerr << "NUll Cell Code\n";
            std::cout << "Failed in [" << i << ", " << j << "]\n";
            std::cout << "cell polygons: ";
            for (std::vector<const Point *> &vec :
                 fully_hori_clipped_vertices) {
              print_vec(vec);
            }
            return 2;
          }
          i_j_to_rcell_info[i_j] =
              RasterCellInfo(fully_hori_clipped_vertices, cell_code);
        }
      }
    }
  }

  return 1;
}
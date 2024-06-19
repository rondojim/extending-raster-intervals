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

// Calculate the smallest number greater than or equal to k in a
// sequence that starts at xmin ends to xmax and increments by step
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

// Calculate the largest number less than or equal to k in a
// sequence that starts at xmin ends to xmax and increments by step
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

// Compute the index of a given value x in a sequence
// that starts at xmin and increments by step
unsigned int RasterGrid::sequence_idx(double x, double xmin) {
  return std::round((x - xmin) /
                    step); // (x-xmin) / step; // std::round(index); Maybe
                           // usefull for floating-point precision error
}

double RasterGrid::idx_to_coord(int idx, double xmin) {
  return xmin + idx * step;
}

void save_vertices(std::vector<const Point *> vertices, const char *output_file,
                   const char *mode) {
  FILE *fp = fopen(output_file, mode);
  for (int i = 0; i < vertices.size(); ++i) {
    fprintf(fp, "%lf %lf\n", vertices[i]->x, vertices[i]->y);
  }
  fclose(fp);
}

int get_nxt_idx(int idx, int size, int step) {
  return (idx + step) < 0 ? size - 1 : (idx + step) % size;
}

bool RasterGrid::weiler_clip_vertices_vectors(
    const std::vector<std::vector<const Point *>> &vertices_vectors,
    std::vector<std::vector<const Point *>> &left_vertices_vectors,
    std::vector<std::vector<const Point *>> &right_vertices_vectors,
    const Point &p1, const Point &p2, bool debug) {

  for (const std::vector<const Point *> &vertices : vertices_vectors) {

    if (!weiler_clip(vertices, left_vertices_vectors, right_vertices_vectors,
                     p1, p2, 100000, debug)) {
      // if (debug) {
      //  FILE *fp = fopen("weiler_error.txt", "w");
      //  fprintf(fp, "segment %lf %lf %lf %lf\n", p1.x, p1.y, p2.x, p2.y);
      //  fclose(fp);

      // save_vertices(vertices, "weiler_error.txt", "a");
      // std::cout << std::fixed << std::setprecision(16);
      // std::cout << "p1: (" << p1.x << ", " << p1.y << "), p2: (" << p2.x
      //           << ", " << p2.y << ")\n";
      // // std::cout << "p1: " << p1.to_str() << ", p2: " << p2.to_str() <<
      // // std::endl;
      // weiler_clip(vertices, left_vertices_vectors, right_vertices_vectors,
      // p1,
      //             p2, 10, true);
      //}

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

int get_double_sigh(double x, double epsilon = 1e-24) {
  if (std::fabs(x - 0.0) < epsilon) {
    return 0;
  } else if (x < 0.0) {
    return -1;
  }
  return 1;
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
    double i_pos = cross_product_(p1, p2, *p_i);
    double nxt_i_pos = cross_product_(p1, p2, *p_nxt_i);
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
    int start_idx, int step, const Point **cur_p_inter, unsigned int max_iters,
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
    cur_idx = get_nxt_idx(cur_idx, inter_pts_info.size(), step);
    cur_dir = inter_pts_info[cur_idx].dir;
  }

  if (total_iters > max_iters) {
    std::cerr << "Error in walk_on_intersections: max_iters == total_iters\n";
    return -1;
  }

  inter_pts_info[cur_idx].visited = true;
  *cur_p_inter = inter_pts_info[cur_idx].inter_point;
  polygon_to_fill.push_back(*cur_p_inter);
  if (debug) {
    std::cout << "cur_p_inter: (" << inter_pts_info[cur_idx].inter_point->x
              << ", " << inter_pts_info[cur_idx].inter_point->y
              << "), dir: " << inter_pts_info[cur_idx].dir << " visited\n";
  }

  return inter_pts_info[cur_idx].nxt_point_idx;
}

void save_error_info(std::vector<std::vector<const Point *>> vertices_vectors,
                     const Point &p1, const Point &p2,
                     const char *output_file) {

  std::cout << "Total polygons in cell: " << vertices_vectors.size() << "\n";
  std::map<unsigned int, std::vector<std::vector<const Point *>>>
      i_j_to_vertices_vector;
  std::pair<unsigned int, unsigned int> i_j = {0, 0};
  i_j_to_vertices_vector[0] = vertices_vectors;
  save_vertices_vectors("wrong_clip.txt", vertices_vectors, "w");

  for (const std::vector<const Point *> &vertices : vertices_vectors) {

    FILE *fp = fopen("weiler_error.txt", "w");
    fprintf(fp, "segment %lf %lf %lf %lf\n", p1.x, p1.y, p2.x, p2.y);
    fclose(fp);

    save_vertices(vertices, "weiler_error.txt", "a");
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "p1: (" << p1.x << ", " << p1.y << "), p2: (" << p2.x << ", "
              << p2.y << ")\n";
  }
}

bool RasterGrid::weiler_rasterize_poly(
    Polygon &polygon,
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        &i_j_to_rcell_info,
    bool debug) {
  double double_null = std::numeric_limits<double>::quiet_NaN();
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
    return false;
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
      vertices_vectors; // in case the whole polygon is contained in the first
                        // columns ( xmin + step > xmax)
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
        // save_vertices_vectors("col_error.txt", j_to_vertices_vectors[j]);
        std::cout << "Failed in j " << j << "\n";
        return false;
      }

      if (!success) {
        std::cout << "Failed in j: " << j << "\n";
        return false;
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
        // i_j_to_vertices_vectors[i_j] = semi_hori_clipped_vertices;
        if (semi_hori_clipped_vertices.size()) {
          BinaryCellCode cell_code =
              encode(semi_hori_clipped_vertices, p1_hori, p2_hori);
          if (cell_code.value == BinaryCellCode::NULL_CODE) {
            save_error_info(semi_hori_clipped_vertices, p1_hori, p2_hori,
                            "wrong_clip.txt");
            std::cerr << "NUll Cell Code\n";
            std::cout << "Failed in [" << i << ", " << j << "]\n";
            std::cout << "cell polygons: ";
            for (std::vector<const Point *> &vec : semi_hori_clipped_vertices) {
              print_vec(vec);
            }
            return false;
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
          return false;
        }

        cur_hori_clipped_vertices = semi_hori_clipped_vertices;
        if (fully_hori_clipped_vertices.size()) {
          BinaryCellCode cell_code =
              encode(fully_hori_clipped_vertices, p1_hori, p2_hori);
          if (cell_code.value == BinaryCellCode::NULL_CODE) {
            save_error_info(fully_hori_clipped_vertices, p1_hori, p2_hori,
                            "wrong_clip.txt");
            std::cerr << "NUll Cell Code\n";
            std::cout << "Failed in [" << i << ", " << j << "]\n";
            std::cout << "cell polygons: ";
            for (std::vector<const Point *> &vec :
                 fully_hori_clipped_vertices) {
              print_vec(vec);
            }
            return false;
          }
          i_j_to_rcell_info[i_j] =
              RasterCellInfo(fully_hori_clipped_vertices, cell_code);
          // i_j_to_vertices_vectors[i_j] = fully_hori_clipped_vertices;
        }
      }
    }
  }

  return true;
}

bool RasterGrid::save_mbr_grid(const char *output_file, Point &mbr_minCorner,
                               Point &mbr_maxCorner, const char *mode) {
  double xmin =
      max_below_or_equal_k(min_corner.x, max_corner.x, mbr_minCorner.x);
  double ymin =
      max_below_or_equal_k(min_corner.y, max_corner.y, mbr_minCorner.y);
  double xmax =
      min_above_or_equal_k(min_corner.x, max_corner.x, mbr_maxCorner.x);
  double ymax =
      min_above_or_equal_k(min_corner.y, max_corner.y, mbr_maxCorner.y);

  if (std::isnan(xmin) || std::isnan(ymin) || std::isnan(xmax) ||
      std::isnan(ymax)) {
    std::cerr << "Error in weiler_rasterize_poly: min_above_or_equal_k or "
                 "max_below_or_equal_k returned NaN\n";
    return false;
  }

  int xmin_idx = sequence_idx(xmin, min_corner.x);
  int xmax_idx = sequence_idx(xmax, min_corner.x);

  int ymin_idx = sequence_idx(ymin, min_corner.y);
  int ymax_idx = sequence_idx(ymax, min_corner.y);

  FILE *fp = fopen(output_file, mode);

  double x_j, y_i;
  for (int j = xmin_idx; j <= xmax_idx; j++) {
    x_j = min_corner.x + j * step;
    for (int i = ymax_idx; i >= ymin_idx; i--) {
      y_i = min_corner.y + i * step;

      fprintf(fp, "%lf %lf\n", x_j, y_i);
    }
  }
  fclose(fp);

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
  if (debug) {
    std::cout << std::fixed << std::setprecision(15);
    std::cout << "p1: (" << p1.x << ", " << p1.y << "), p2: (" << p2.x << ", "
              << p2.y << "), p3: (" << p3.x << ", " << p3.y << "), p4: ("
              << p4.x << ", " << p4.y << ")" << std::endl;
    print_vec(vertices);
  }
  // Check if all cell corners are in the vertices vector
  for (const Point &corner : cell_corners) {
    bool found = false;
    for (const Point *p : vertices) {
      if (*p == corner) {
        if (debug) {
          std::cout << "p: (" << p->x << ", " << p->y << "), corner: ("
                    << corner.x << ", " << corner.y << ")\n";
        }
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

double
get_polygons_area(std::vector<std::vector<const Point *>> &vertices_vectors) {
  double polys_area = 0.0;
  for (std::vector<const Point *> &vertices : vertices_vectors) {
    polys_area += polygon_area(vertices);
    ;
  }
  return polys_area;
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

void RasterGrid::split_points(
    double i_pos, double nxt_i_pos, const Point *p_inter, const Point *p_nxt_i,
    std::vector<const Point *> &left_vert_clipped_points,
    std::vector<const Point *> &right_vert_clipped_points, bool debug) {
  // Case 1: For the right side: When both points are inside
  if (i_pos < 0 && nxt_i_pos < 0) {
    if (debug) {
      std::cout << "case 1.1\n";
    }
    // Only second point is added
    right_vert_clipped_points.push_back(p_nxt_i);
  }
  // Case 2: For the right side:
  // First point is outside or on line
  // Second is inside
  else if (is_greater_or_equal(i_pos, 0.0) && nxt_i_pos < 0) {
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
  else if (i_pos < 0 && is_greater_or_equal(nxt_i_pos, 0.0)) {
    if (debug) {
      std::cout << "case 1.3\n";
    }
    // Only point of intersection with edge is added
    right_vert_clipped_points.push_back(p_inter);
  }
  // Case 4: For the right side: When both points are outside, no points are
  // added

  // Case 1: For the left side: Both points are inside
  if (i_pos > 0 && nxt_i_pos > 0) {
    if (debug) {
      std::cout << "case 2.1\n";
    }
    left_vert_clipped_points.push_back(p_nxt_i);
  }
  // Case 2: First vertex is outside while second one is inside
  else if (is_less_or_equal(i_pos, 0.0) && nxt_i_pos > 0) {
    if (debug) {
      std::cout << "case 2.2\n";
    }
    // For the left clipped point of intersection with edge and the second point
    // is added
    left_vert_clipped_points.push_back(p_inter);
    left_vert_clipped_points.push_back(p_nxt_i);
  }
  // Case 3: First vertex is inside while second one is outside
  else if (i_pos > 0 && is_less_or_equal(nxt_i_pos, 0.0)) {
    if (debug) {
      std::cout << "case 2.3\n";
    }
    // For the fully clipped only point of intersection with edge is added
    left_vert_clipped_points.push_back(p_inter);
  }
}

bool RasterGrid::clip(const std::vector<const Point *> &vertices,
                      std::vector<const Point *> &left_vert_clipped_points,
                      std::vector<const Point *> &right_vert_clipped_points,
                      const Point &p1, const Point &p2, bool debug) {

  int poly_size = vertices.size();
  for (int i = 0; i < poly_size - 1; i++) {
    // i and nxt_i form a line in polygon
    int nxt_i = (i + 1) % poly_size;
    const Point *p_i = vertices[i];
    const Point *p_nxt_i = vertices[nxt_i];

    // Calculating position of first, second point w.r.t. clipper line
    double i_pos = cross_product_(p1, p2, *p_i);
    double nxt_i_pos = cross_product_(p1, p2, *p_nxt_i);

    Point *p_inter = nullptr;
    // if the points are on different sides wrt to the p1-p2 line or at least
    // one of them is on the line
    if (i_pos * nxt_i_pos <= 0) {
      p_inter = p_intersect(p1, p2, *p_i, *p_nxt_i);
      if (p_inter == nullptr) {
        std::cerr << "Nan intersection\n";
        return false;
      }
      allocated_points.push_back(p_inter);
    }

    split_points(i_pos, nxt_i_pos, p_inter, p_nxt_i, left_vert_clipped_points,
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

int RasterGrid::clip_check_full(
    const std::vector<const Point *> &vertices,
    std::vector<const Point *> &left_vert_clipped_points,
    std::vector<const Point *> &right_vert_clipped_points, const Point &p1,
    const Point &p2, bool debug, bool debug2) {

  int poly_size = vertices.size();
  double step = p2.x - p1.x;
  Point p3(p2.x, p2.y + step);
  Point p4(p1.x, p1.y + step);
  std::array<bool, 4> point_in_corner = {false, false, false, false};
  bool collinear_with_edge = true;

  for (int i = 0; i < poly_size - 1; i++) {
    // i and nxt_i form a line in polygon
    int nxt_i = (i + 1) % poly_size;
    const Point *p_i = vertices[i];
    const Point *p_nxt_i = vertices[nxt_i];

    if (collinear_with_edge && p_i->y >= p1.y) {
      unsigned int equal_corner_idx = equal_with_corner(*p_i, p1, p2, p3, p4);
      if (equal_corner_idx) {
        point_in_corner[equal_corner_idx - 1] = true;
      } else if (!collinear_with_any_edge(*p_i, p1, p2, p3, p4)) {
        collinear_with_edge = false;
      }
    }

    // Calculating position of first, second point w.r.t. clipper line
    double i_pos = cross_product_(p1, p2, *p_i);
    double nxt_i_pos = cross_product_(p1, p2, *p_nxt_i);

    Point *p_inter = nullptr;
    // if the points are on different sides wrt to the p1-p2 line or at least
    // one of them is on the line
    if (i_pos * nxt_i_pos <= 0) {
      p_inter = p_intersect(p1, p2, *p_i, *p_nxt_i);
      if (p_inter == nullptr) {
        std::cerr << "Nan intersection\n";
        return -1;
      }
      allocated_points.push_back(p_inter);

      if (collinear_with_edge) {
        unsigned int equal_corner_idx =
            equal_with_corner(*p_inter, p1, p2, p3, p4);
        if (equal_corner_idx) {
          point_in_corner[equal_corner_idx - 1] = true;
        }
      }
    }

    split_points(i_pos, nxt_i_pos, p_inter, p_nxt_i, left_vert_clipped_points,
                 right_vert_clipped_points, debug);
  }

  // Append first point as last to include the edge between them
  if (left_vert_clipped_points.size()) {
    left_vert_clipped_points.push_back(left_vert_clipped_points[0]);
  }

  if (right_vert_clipped_points.size()) {
    right_vert_clipped_points.push_back(right_vert_clipped_points[0]);
  }

  if (collinear_with_edge) {
    for (bool &pt_in_corner_i : point_in_corner) {
      if (!pt_in_corner_i) {
        return 0;
      }
    }
    return 1;
  }
  return 0;
}

void print_vec(const std::vector<const Point *> &vec) {
  for (const Point *p : vec) {
    std::cout << std::fixed << std::setprecision(15);
    std::cout << "(" << p->x << ", " << p->y << "), ";
  }
  std::cout << std::endl;
}

// //the whole polygon is contained in a cell and specifically the down left
// bool RasterGrid::rasterize_poly(const std::vector<const Point*>& vertices,
//                Point& mbr_min_corner, Point& mbr_max_corner,
//                std::map<std::pair<unsigned int, unsigned int>,
//                RasterCellInfo>& i_j_to_rcell_info, bool debug){

//     double xmin = max_below_or_equal_k(g_xmin, g_xmax, mbr_min_corner.x);
//     double ymin = max_below_or_equal_k(g_ymin, g_ymax, mbr_min_corner.y);
//     double xmax = min_above_or_equal_k(g_xmin, g_xmax, mbr_max_corner.x);
//     double ymax = min_above_or_equal_k(g_ymin, g_ymax, mbr_max_corner.y);

//     int xmin_idx = sequence_idx(xmin, g_xmin);
//     int xmax_idx = sequence_idx(xmax, g_xmin);

//     int ymin_idx = sequence_idx(ymin, g_ymin);
//     int ymax_idx = sequence_idx(ymax, g_ymin);

//     std::vector<const Point*> fully_vert_clipped_vertices,
//     semi_vert_clipped_vertices, cur_vert_clipped_vertices,
//     fully_hori_clipped_vertices, semi_hori_clipped_vertices,
//     cur_hori_clipped_vertices; cur_vert_clipped_vertices = vertices;
//     semi_vert_clipped_vertices = vertices; // in case the whole polygon is
//     contained in the first columns ( xmin + step > xmax) double x_j, y_i;
//     std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
//     test_grid, test_grid2; bool debug1=false, debug2 = false;

//     for (int j = xmin_idx; j < xmax_idx; j++) {
//         // std::cout << "j: " << j << std::endl;
//         x_j = g_xmin + (j+1)*step;
//         fully_vert_clipped_vertices.clear();

//         if (j == xmax_idx - 1) {
//             cur_hori_clipped_vertices = semi_vert_clipped_vertices;
//         }
//         else {
//             semi_vert_clipped_vertices.clear();
//             Point p1_vert(x_j, ymin);
//             Point p2_vert(x_j, ymax);
//             bool success = clip(cur_vert_clipped_vertices,
//             fully_vert_clipped_vertices, semi_vert_clipped_vertices, p1_vert,
//             p2_vert, debug);

//             if (!success) {
//                 return false;
//             }

//             cur_vert_clipped_vertices = semi_vert_clipped_vertices;
//             cur_hori_clipped_vertices = fully_vert_clipped_vertices;
//         }

//         semi_hori_clipped_vertices = cur_hori_clipped_vertices;

//         // if (debug) {
//             std::pair<int, int> i_j = {j, j};
//             std::vector<const Point*> test = semi_hori_clipped_vertices;
//             test_grid[i_j] = RasterCellInfo(test,
//             BinaryCellCode(BinaryCellCode::NULL_CODE));
//         // }

//         for (int i = ymax_idx - 1; i >= ymin_idx; i--) {
//             bool full_type_cell = true;
//             // if (i==3 && j==2){ // && j==2 || i == 2 && j == 2){
//             //     std::cout << "i: " << i << ", j: " << j << std::endl;
//             //     // print_vec(cur_hori_clipped_vertices);
//             //     debug2 = true;
//             //     debug1 = true;
//             // }
//             // else {
//             //     debug1 = false;
//             //     debug2 = false;
//             // }
//             y_i = g_ymin + i * step;
//             fully_hori_clipped_vertices.clear();

//             std::pair<int, int> i_j = {i, j};

//             if (i == ymin_idx) {
//                 // i_j_to_clipped_vertices[i_j] = semi_hori_clipped_vertices;
//                 // i_j_to_cell_type[i_j] =
//                 encode(semi_hori_clipped_vertices); i_j_to_rcell_info[i_j] =
//                 RasterCellInfo(semi_hori_clipped_vertices,
//                 encode(semi_hori_clipped_vertices));
//                 // if (debug) {
//                 //     test_grid2[i_j] = semi_hori_clipped_vertices;
//                 // }
//             }
//             else {
//                 semi_hori_clipped_vertices.clear();
//                 Point p1_hori(x_j - step, y_i);
//                 Point p2_hori(x_j, y_i);

//                 int full_cell_type =
//                 clip_check_full(cur_hori_clipped_vertices,
//                 fully_hori_clipped_vertices, semi_hori_clipped_vertices,
//                 p1_hori, p2_hori, debug1, debug2); if (full_cell_type == -1)
//                 {
//                     return false;
//                 }

//                 cur_hori_clipped_vertices = semi_hori_clipped_vertices;
//                 BinaryCellCode real_cell_code =
//                 encode(fully_hori_clipped_vertices); BinaryCellCode
//                 cell_code;
//                 // i_j_to_cell_type[i_j] = cell_code;
//                 if (full_cell_type == 1 && real_cell_code.to_type() !=
//                 "FULL") {
//                     std::cout << "i: " << i << ", j: " << j << " -> should
//                     NOT be FULL\n";

//                 }

//                 if (real_cell_code.to_type() == "FULL" && !full_cell_type){
//                     std::cout << "i: " << i << ", j: " << j << " -> should be
//                     FULL\n";
//                 }

//                 if (full_cell_type) {
//                     if (true) {
//                         cell_code = BinaryCellCode(BinaryCellCode::FULL_R);
//                     }
//                     else {
//                         cell_code = BinaryCellCode(BinaryCellCode::FULL_S);
//                     }
//                 }
//                 else{
//                     cell_code = encode(fully_hori_clipped_vertices);
//                 }

//                 i_j_to_rcell_info[i_j] =
//                 RasterCellInfo(fully_hori_clipped_vertices, cell_code);

//                 //
//                 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//                 // for (auto p : allocated_points) {
//                 //     delete p;
//                 // }
//                 // allocated_points.clear();
//                 if (debug) {
//                     test_grid2[i_j] =
//                     RasterCellInfo(fully_hori_clipped_vertices,
//                     BinaryCellCode(BinaryCellCode::NULL_CODE));
//                 }
//             }

//         }

//     }
//     save_rasterization("test.txt", test_grid);

//     return true;

//     // save_rasterization("test2.txt", test_grid2);
// }

void save_rasterization(
    const char *output_file,
    const std::map<std::pair<unsigned int, unsigned int>,
                   const std::vector<const Point *>> &i_j_to_clipped_vertices,
    const char *mode) {

  FILE *fp = fopen(output_file, mode);
  for (const auto &entry : i_j_to_clipped_vertices) {
    for (const Point *point : entry.second) {
      fprintf(fp, "%lf %lf\n", point->x, point->y);
    }
    fprintf(fp, "poly %d %d\n", entry.first.first, entry.first.second);
  }
  fclose(fp);
}

void save_rasterization(const char *output_file,
                        const std::map<std::pair<unsigned int, unsigned int>,
                                       RasterCellInfo> &i_j_to_rcell_info,
                        const char *mode) {

  FILE *fp = fopen(output_file, mode);
  for (const auto &entry : i_j_to_rcell_info) {
    unsigned int i = entry.first.first;
    unsigned int j = entry.first.second;
    const RasterCellInfo &rcell_info = entry.second;

    for (const std::vector<const Point *> &vec : rcell_info.cell_polygons) {
      for (const Point *point : vec) {
        fprintf(fp, "%lf %lf\n", point->x, point->y);
      }

      fprintf(fp, "poly\n");
    }

    fprintf(fp, "cell %d %d\n", i, j);
    fprintf(fp, "%s\n", rcell_info.cell_type.to_type());
  }
  fclose(fp);
}

void save_vertices_vectors(
    const char *output_file,
    std::map<std::pair<unsigned int, unsigned int>,
             std::vector<std::vector<const Point *>>> &i_j_to_vertices_vectors,
    const char *mode) {

  FILE *fp = fopen(output_file, mode);

  for (const auto &entry : i_j_to_vertices_vectors) {
    unsigned int i = entry.first.first;
    unsigned int j = entry.first.second;
    std::vector<std::vector<const Point *>> vertices_vectors = entry.second;

    for (std::vector<const Point *> &vec : vertices_vectors) {
      for (const Point *point : vec) {
        fprintf(fp, "%lf %lf\n", point->x, point->y);
      }

      fprintf(fp, "poly\n");
    }
    fprintf(fp, "cell %d %d\n", i, j);
  }
  fclose(fp);
}

void save_vertices_vectors(
    const char *output_file,
    std::map<unsigned int, std::vector<std::vector<const Point *>>>
        &i_j_to_vertices_vectors,
    const char *mode) {

  FILE *fp = fopen(output_file, mode);

  for (const auto &entry : i_j_to_vertices_vectors) {
    unsigned int j = entry.first;
    std::vector<std::vector<const Point *>> vertices_vectors = entry.second;

    for (std::vector<const Point *> &vec : vertices_vectors) {
      for (const Point *point : vec) {
        fprintf(fp, "%lf %lf\n", point->x, point->y);
      }

      fprintf(fp, "poly\n");
    }
    fprintf(fp, "cell %d %d\n", j, j);
  }
  fclose(fp);
}

void save_vertices_vectors(
    const char *output_file,
    std::vector<std::vector<const Point *>> &vertices_vectors,
    const char *mode) {

  FILE *fp = fopen(output_file, mode);

  for (int i = 0; i < vertices_vectors.size(); i++) {
    std::vector<const Point *> vec = vertices_vectors[i];
    for (const Point *point : vec) {
      fprintf(fp, "%lf %lf\n", point->x, point->y);
    }

    fprintf(fp, "poly\n");
  }
  fprintf(fp, "cell %d %d\n", 0, 0);
  fclose(fp);
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

/////////////////////////////////////////////// BinaryCellCode Functions
/////////////////////////////////////////

BinaryCellCode::BinaryCellCode(unsigned char val, bool type_R_)
    : value(val & 0b111), type_R(type_R_) {}

BinaryCellCode::BinaryCellCode(Code code, bool type_R_)
    : value(code & 0b111), type_R(type_R_) {}

// Bitwise AND
BinaryCellCode BinaryCellCode::operator&(const BinaryCellCode &other) const {
  return BinaryCellCode((this->value & other.value).to_ulong());
}

// To string
std::string BinaryCellCode::to_string() const { return value.to_string(); }

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

bool rasterize_polygons(RasterGrid &grid, std::vector<Polygon> &lhs_polygons,
                        std::vector<Polygon> &rhs_polygons,
                        std::vector<RasterPolygonInfo> &lhs_i_j_to_rpoly_info,
                        std::vector<RasterPolygonInfo> &rhs_i_j_to_rpoly_info,
                        std::string err_poly_f_name, bool debug) {

  progressbar bar(lhs_polygons.size() + rhs_polygons.size());

  for (Polygon &polygon : lhs_polygons) {
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        i_j_to_rcell_info;
    if (!grid.weiler_rasterize_poly(polygon, i_j_to_rcell_info, debug)) {
      if (err_poly_f_name.size()) {
        polygon.save_poly(err_poly_f_name.c_str());
        return false;
      }
    }
    bar.update();
    printf("\n%d polygon_id\n", polygon.polygon_id);
    // printf polygon
    print_vec(polygon.vertices);
    for (auto &entry : i_j_to_rcell_info) {
      printf("i: %d, j: %d\n", entry.first.first, entry.first.second);
      printf("cell type: %s\n", entry.second.cell_type.to_type());
    }
    lhs_i_j_to_rpoly_info.emplace_back(polygon.polygon_id, i_j_to_rcell_info);
  }

  printf("done with lhs\n");
  for (Polygon &polygon : rhs_polygons) {
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
        i_j_to_rcell_info;
    if (!grid.weiler_rasterize_poly(polygon, i_j_to_rcell_info, debug)) {
      if (err_poly_f_name.size()) {
        polygon.save_poly(err_poly_f_name.c_str());
      }
      return false;
    }
    rhs_i_j_to_rpoly_info.emplace_back(polygon.polygon_id, i_j_to_rcell_info);
    bar.update();
  }
  return true;
}
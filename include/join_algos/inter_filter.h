#ifndef INTER_FILTER_H
#define INTER_FILTER_H
#include "../../include/utils/geometry_types.h"
#include <bitset>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

enum Direction {
  RIGHT = 1,
  LEFT = 2,
  NONE = 0,
};

struct PointInfo {
  const Point *point;
  Direction dir;
  int inter_point_idx;

  PointInfo(const Point *p, Direction d, int inter_p_idx = -1)
      : point(p), dir(d), inter_point_idx(inter_p_idx) {}
  void print_(bool line_change = true) const {
    std::cout << "p: (" << point->x << ", " << point->y << "), dir: " << dir
              << ", inter_point_idx: " << inter_point_idx;
    if (line_change) {
      std::cout << std::endl;
    }
  }
};

struct InterPointInfo {
  const Point *inter_point;
  Direction dir;
  int nxt_point_idx;
  bool visited;

  InterPointInfo(const Point *p, Direction d, int nxt_p_idx = -1,
                 bool visited_ = false)
      : inter_point(p), dir(d), nxt_point_idx(nxt_p_idx), visited{(visited_)} {}
  void print_(bool line_change = true) const {
    std::cout << "inter_point: (" << inter_point->x << ", " << inter_point->y
              << "), dir: " << dir << ", nxt_point_idx: " << nxt_point_idx;
    if (line_change) {
      std::cout << std::endl;
    }
  }
};

struct BinaryCellCode {
  std::bitset<3> value;
  bool type_R;

  enum Code {
    FULL_R = 0b011,
    STRONG_R = 0b101,
    WEAK_R = 0b100,
    FULL_S = 0b101,
    STRONG_S = 0b011,
    WEAK_S = 0b010,
    NULL_CODE = 0b000,
  };

  BinaryCellCode(unsigned char val = 0, bool type_R_ = true);

  BinaryCellCode(Code code, bool type_R_ = true);

  // Bitwise AND
  BinaryCellCode operator&(const BinaryCellCode &other) const;

  // To string
  std::string to_string() const;

  // Get value
  unsigned long to_ulong() const;

  const char *to_type(bool debug = false) const;
};

struct RasterCellInfo {
  std::vector<std::vector<const Point *>> cell_polygons;
  BinaryCellCode cell_type;

  RasterCellInfo(std::vector<std::vector<const Point *>> cell_polygons_,
                 BinaryCellCode cell_type_)
      : cell_polygons(std::move(cell_polygons_)),
        cell_type(std::move(cell_type_)) {}

  RasterCellInfo() = default;
};

struct RasterPolygonInfo {
  int idx;
  std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
      i_j_to_rcell_info;

  RasterPolygonInfo(int idx_, std::map<std::pair<unsigned int, unsigned int>,
                                       RasterCellInfo> &i_j_to_rcell_info_)
      : idx(std::move(idx_)),
        i_j_to_rcell_info(std::move(i_j_to_rcell_info_)) {};
};

struct RasterGrid {
private:
  // Calculate the smallest number greater than or equal to k in a
  // sequence that starts at xmin ends to xmax and increments by step
  double min_above_or_equal_k(double xmin, double xmax, double k);

  // Calculate the largest number less than or equal to k in a
  // sequence that starts at xmin ends to xmax and increments by step
  double max_below_or_equal_k(double xmin, double xmax, double k);

  // Compute the index of a given value x in a sequence
  // that starts at xmin and increments by step
  unsigned int sequence_idx(double x, double xmin);

  double idx_to_coord(int idx, double xmin);

  void split_points(double i_pos, double nxt_i_pos, const Point *p_inter,
                    const Point *p_nxt_i,
                    std::vector<const Point *> &left_vert_clipped_points,
                    std::vector<const Point *> &right_vert_clipped_points,
                    bool debug);

  bool is_certain_full_cell(const std::vector<const Point *> &vertices,
                            const Point &p1, const Point &p2,
                            bool debug = false);

  BinaryCellCode
  encode(std::vector<std::vector<const Point *>> &vertices_vectors,
         const Point &p1, const Point &p2, 
         bool debug = false);

  // clip the polygon with vertices elements (should be in cw order)
  // w.r.t. clipper line with points x1, y1, x2, y2
  // Insert points on the left of the line in left_vert_clipped_points
  // and points on the right of the line in right_vert_clipped_points
  bool clip(const std::vector<const Point *> &vertices,
            std::vector<const Point *> &left_vert_clipped_points,
            std::vector<const Point *> &right_vert_clipped_points,
            const Point &p1, const Point &p2, bool debug);

  int clip_check_full(const std::vector<const Point *> &vertices,
                      std::vector<const Point *> &left_vert_clipped_points,
                      std::vector<const Point *> &right_vert_clipped_points,
                      const Point &p1, const Point &p2, bool debug = false,
                      bool debug2 = false);

  int weiler_scan_polygon(const std::vector<const Point *> &vertices,
                          const Point &p1, const Point &p2,
                          std::vector<PointInfo> &vertices_with_info,
                          std::vector<InterPointInfo> &inter_pts_info,
                          bool debug);

  int walk_on_polygon_vertices(std::vector<const Point *> &polygon_to_fill,
                               std::vector<PointInfo> &vertices_with_info,
                               Direction exit_dir, int start_idx,
                               unsigned int max_iters, bool debug = false);

  int walk_on_inter_vertices(std::vector<const Point *> &polygon_to_fill,
                             std::vector<InterPointInfo> &inter_pts_info,
                             Direction exit_dir, int start_idx, int step,
                             const Point **exit_point, unsigned int max_iters,
                             bool debug = false);

  bool
  weiler_clip(const std::vector<const Point *> &vertices,
              std::vector<std::vector<const Point *>> &left_vertices_vectors,
              std::vector<std::vector<const Point *>> &right_vertices_vectors,
              const Point &p1, const Point &p2, unsigned int max_iters,
              bool debug = false);

  bool weiler_clip_vertices_vectors(
      const std::vector<std::vector<const Point *>> &vertices_vectors,
      std::vector<std::vector<const Point *>> &left_vertices_vectors,
      std::vector<std::vector<const Point *>> &right_vertices_vectors,
      const Point &p1, const Point &p2, bool debug = false);

public:
  std::vector<Point *> allocated_points;
  double n, step, cell_area;
  Point min_corner, max_corner;

  RasterGrid(double n_, Point min_corner_, Point max_corner_);

  bool weiler_rasterize_poly(Polygon &polygon,
                             std::map<std::pair<unsigned int, unsigned int>,
                                      RasterCellInfo> &i_j_to_rcell_info,
                             bool debug = false);
  bool save_mbr_grid(const char* output_file, Point& mbr_minCorner, Point& mbr_maxCorner, const char* mode="w");


  // bool weiler_rasterize_poly(
  // const std::vector<const Point *> &vertices,
  // Point &mbr_min_corner, Point &mbr_max_corner,
  // std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
  // &i_j_to_rcell_info, bool debug = true);
};

void save_rasterization(const char *output_file,
                        const std::map<std::pair<unsigned int, unsigned int>,
                                       RasterCellInfo> &i_j_to_rcell_info,
                        const char *mode = "w");

void save_rasterization(
    const char *output_file,
    const std::map<std::pair<unsigned int, unsigned int>,
                   const std::vector<const Point *>> &i_j_to_clipped_vertices,
    const char *mode);

void print_vec(const std::vector<const Point *> &vec);

bool are_equal(double a, double b, double epsilon = 1e-20);
bool is_greater_or_equal(double a, double b, double epsilon = 1e-20);
bool is_less_or_equal(double a, double b, double epsilon = 1e-20);

void save_vertices(std::vector<const Point *> vertices, const char *output_file,
                   const char *mode = "w");

void save_vertices_vectors(
    const char *output_file,
    std::map<std::pair<unsigned int, unsigned int>,
             std::vector<std::vector<const Point *>>> &i_j_to_vertices_vectors,
    const char *mode = "w");

void save_vertices_vectors(
    const char *output_file,
    std::map<unsigned int, std::vector<std::vector<const Point *>>>
        &i_j_to_vertices_vectors,
    const char *mode = "w");

void save_vertices_vectors(
    const char *output_file,
    std::vector<std::vector<const Point *>> &vertices_vectors,
    const char *mode = "w");

int get_nxt_idx(int idx, int size, int step);

double
get_polygons_area(std::vector<std::vector<const Point *>> &vertices_vectors);

bool rasterize_polygons(RasterGrid &grid, std::vector<Polygon> &lhs_polygons,
                        std::vector<Polygon> &rhs_polygons,
                        std::vector<RasterPolygonInfo> &lhs_i_j_to_rpoly_info,
                        std::vector<RasterPolygonInfo> &rhs_i_j_to_rpoly_info,
                        std::string err_poly_f_name, bool debug);

#endif
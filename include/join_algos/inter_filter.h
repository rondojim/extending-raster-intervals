#ifndef INTER_FILTER_H
#define INTER_FILTER_H
#include "../../include/utils/geometry_types.h"
#include <bitset>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
#include <set>

// expresses the direction of interesctio 
// of a segment wrt another segment
enum Direction {
  RIGHT = 1,
  LEFT = 2,
  NONE = 0,
};

struct PointInfo {
  const Point *point; // the polygon point
  Direction dir; // the direction in case of intersection of the polygon segment wrt to the clipping segment
  int inter_point_idx; // the next intersection point index

  PointInfo(const Point *p, Direction d, int inter_p_idx = -1)
      : point(p), dir(d), inter_point_idx(inter_p_idx) {}
  void print_(bool line_change = true) const {
    if (line_change) {
      std::cout << std::endl;
    }
  }
};

struct InterPointInfo {
  const Point *inter_point; // the interesction point
  Direction dir; // the direction of the polygon segment wrt to the clipping segment for which we have interesction
  int nxt_point_idx; // the next point index on the polygon vertices
  bool visited; // wether the interesction point is visited or not

  InterPointInfo(const Point *p, Direction d, int nxt_p_idx = -1,
                 bool visited_ = false)
      : inter_point(p), dir(d), nxt_point_idx(nxt_p_idx), visited{(visited_)} {}
  void print_(bool line_change = true) const {
    if (line_change) {
      std::cout << std::endl;
    }
  }
};

struct BinaryCellCode {
  std::bitset<3> value; // the cell type is represented in binary form of 3 bits
  bool type_R; // whether the cell is of dataset of type R or S 

  // the possible codes of a cell
  // NULL_CODE indicates error
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

  // Get value in ulong
  unsigned long to_ulong() const;

  // convert binary to char* FULL, STRONG, WEAK, NULL
  const char *to_type(bool debug = false) const;
};

struct RasterCellInfo {
  std::vector<std::vector<const Point *>> cell_polygons; // the cell polygons after clipping
  BinaryCellCode cell_type; // the cell type

  RasterCellInfo(std::vector<std::vector<const Point *>> cell_polygons_={},
                 BinaryCellCode cell_type_=BinaryCellCode())
      : cell_polygons(cell_polygons_),
        cell_type(cell_type_) {}
};

struct RasterPolygonInfo {
  int idx; // the index of a polygon
  std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo>
      i_j_to_rcell_info; // a map from i, j to RasterCellInfo

  RasterPolygonInfo(int idx_, std::map<std::pair<unsigned int, unsigned int>,
                                       RasterCellInfo> &i_j_to_rcell_info_)
      : idx(idx_),
        i_j_to_rcell_info(i_j_to_rcell_info_) {};
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

  // returns true if k is part of the sequence start from xmin
  // to xmax with the step of the grid
  int k_belongs_in_sequence(double xmin, double xmax, double k, double epsilon = 1e-6);

  // double idx_to_coord(int idx, double xmin);

  // i_pos, nxt_i_pos describe the side of a point and each next: 
  // wrt to another segment: < 0 on the left, > 0 on the right else on
  // the other segment. Based on that insert points in left_vert_clipped_points,
  // right_vert_clipped_points according to hodgman algo. For the first vector 
  // the visible side is the left while fot the second the right
  void hodgman_clip_segment(int i_pos_sigh, int nxt_i_pos_sigh, const Point *p_inter,
                    const Point *p_nxt_i,
                    std::vector<const Point *> &left_vert_clipped_points,
                    std::vector<const Point *> &right_vert_clipped_points,
                    bool debug);

  // if the vertices form the square with bottom corners p1,p2
  // returns true, else false 
  bool is_certain_full_cell(const std::vector<const Point *> &vertices,
                            const Point &p1, const Point &p2,
                            bool debug = false);

  // clip the polygon with vertices elements (should be in cw order)
  // w.r.t. clipper line with points x1, y1, x2, y2
  // Insert points on the left of the line in left_vert_clipped_points
  // and points on the right of the line in right_vert_clipped_points
  bool hodgman_clip(const std::vector<const Point *> &vertices,
            std::vector<const Point *> &left_vert_clipped_points,
            std::vector<const Point *> &right_vert_clipped_points,
            const Point &p1, const Point &p2, bool debug);


  // we find intersections between vertices and the segment p1,p2
  // and split the vertices, interescions in vertices_with_info, 
  // inter_pts_info according to weiler clip algorithm
  int weiler_scan_polygon(const std::vector<const Point *> &vertices,
                          const Point &p1, const Point &p2,
                          std::vector<PointInfo> &vertices_with_info,
                          std::vector<InterPointInfo> &inter_pts_info,
                          bool debug);

  // walk on the polygon vertices (vertices_with_info) starting from start_idx,
  // adding points until meet exit_dir 
  // max_iters are added for case of error (infinite loop)
  // Returns the entering (on intersection points) idx on success,
  // else -1
  int walk_on_polygon_vertices(std::vector<const Point *> &polygon_to_fill,
                               std::vector<PointInfo> &vertices_with_info,
                               Direction exit_dir, int start_idx,
                               unsigned int max_iters, bool debug = false);

  // walk on the intersection vertices (inter_pts_info) starting 
  // from start_idx, go on the direction based on step (-1/1) up 
  // or down adding points in polygon_to_fill until meet exit_dir 
  // saving the exiting point in exit_point
  // max_iters are added for case of error (infinite loop)
  // Returns the entering (on polygon) point index on success, else -1
  int walk_on_inter_vertices(std::vector<const Point *> &polygon_to_fill,
                             std::vector<InterPointInfo> &inter_pts_info,
                             Direction exit_dir, int start_idx, int step,
                             const Point **exit_point, unsigned int max_iters,
                             bool debug = false);

  // It clips the vertices wrt to the segment p1, p2 and 
  // the result is save in left_vertices_vectors, right_vertices_vectors
  // The weiler clip algo is used and adjusted to clip wrt to a segment
  // Returns true on success, else false
  bool weiler_clip(const std::vector<const Point *> &vertices,
              std::vector<std::vector<const Point *>> &left_vertices_vectors,
              std::vector<std::vector<const Point *>> &right_vertices_vectors,
              const Point &p1, const Point &p2, unsigned int max_iters,
              bool debug = false);

  // It clips each of the vertices included in vertices_vectors 
  // using the weiler_clip function 
  // Returns true on success, else false
  bool weiler_clip_vertices_vectors(
      const std::vector<std::vector<const Point *>> &vertices_vectors,
      std::vector<std::vector<const Point *>> &left_vertices_vectors,
      std::vector<std::vector<const Point *>> &right_vertices_vectors,
      const Point &p1, const Point &p2, bool debug = false);

  // for each segment in vertices we call set_segment_borders_types
  // seting weak cell types for edge cases that are saved in i_j_to_rcell_info
  bool set_polygon_borders_cell_types(std::vector<const Point*> &vertices,
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo> &i_j_to_rcell_info);

  // we check which cell sides of the grid the segment of p1,p2  
  // touches. For each such side we set weak cell types 
  // neighbooring cells of the side which are saved in i_j_to_rcell_info
  bool set_segment_borders_types(const Point &p1, const Point &p2, 
    std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo> &i_j_to_rcell_info);

  // return the R encoding of the vertices_vectors which is
  // in the cell with the bottom corners p1,p2
  // the clipped polygon wrt to the cell area 
  BinaryCellCode
  encode(std::vector<std::vector<const Point *>> &vertices_vectors,
         const Point &p1, const Point &p2, 
         bool debug = false);


public:
  std::vector<Point *> allocated_points;
  //n: grid of size nXn
  // min_corner,max_corner the cornerd of the grid
  // step: the cells' side  
  double n, step, cell_area;
  Point min_corner, max_corner;

  // The side of the grid is set to 2^n_.
  // Given the min_corner_, max_corner_
  // we find the biggest side between x or y
  // and expand the smallest side for the corners
  // of the grid to form a square. Cell side
  // of the grid is grid_side/2^n_
  RasterGrid(double n_, Point min_corner_, Point max_corner_);

  // rasterixe the polygon and save result for each cell 
  // in map with key the cell pos i,j and value RasterCellInfo
  // using the weiler clipping algorithm
  // return 2 in case of null cell code,
  // 0 for any other error
  // 1, on success
  int weiler_rasterize_poly(Polygon &polygon,
                             std::map<std::pair<unsigned int, unsigned int>,
                                      RasterCellInfo> &i_j_to_rcell_info,
                             bool debug = false);

  // rasterixe the polygon and save result for each cell 
  // in map with key the cell pos i,j and value RasterCellInfo
  // using the hodgman clipping algorithm
  bool hodgman_rasterize_poly(Polygon &polygon,
                             std::map<std::pair<unsigned int, unsigned int>,
                                      RasterCellInfo> &i_j_to_rcell_info,
                             bool debug = false);

  // save the rasterization of a poly in output_file in given mode
  // format:
  // xmin_idx, ymin_idx, xmax_idx, ymax_idx // idx of the polygon mbr on the grid
  // min_corner.x, min_corner.y, max_corner.x, max_corner.y, step // of the grid
  // i j cell_type // eg 0 0 FULL
  // i' j' cell_type'
  // ...
  // x0 y0 // the write the polygon vertices
  // x1 y1
  // xn yn
  bool save_poly_raster(const char *output_file, Polygon& poly,
      std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo> &i_j_to_rcell_info, const char *mode="w");

  
  // TO BE IMPLEMENTED
  // IS NOT COMPLETED
  bool save_polygons_grid(const char *output_file, std::vector<Polygon>& polygons, 
   std::map<std::pair<unsigned int, unsigned int>, RasterCellInfo> &i_j_to_rcell_info, const char *mode);
};

// rasterize the lhs_polygons, rhs_polygons and save their 
// raster representation in the lhs_i_j_to_rpoly_info, rhs_i_j_to_rpoly_info
bool rasterize_polygons(RasterGrid &grid, std::vector<Polygon> &lhs_polygons,
                        std::vector<Polygon> &rhs_polygons,
                        std::vector<RasterPolygonInfo> &lhs_i_j_to_rpoly_info,
                        std::vector<RasterPolygonInfo> &rhs_i_j_to_rpoly_info,
                        std::set<int> &null_cell_code_poly_idxs,
                        std::string err_poly_f_name="", bool debug=false);

int get_double_sigh(double x, double epsilon = 1e-24);
bool are_equal(double a, double b, double epsilon = 1e-20);
bool is_greater_or_equal(double a, double b, double epsilon = 1e-20);
bool is_less_or_equal(double a, double b, double epsilon = 1e-20);


// write vertices in output_file in given mode
// in format 
// x0 y0
// x1 y1
// xn yn
void save_vertices(std::vector<const Point *> vertices, const char *output_file,
                   const char *mode = "w");
// write vertices_vectors in output_file in given mode
// in format 
// x0 y0
// x1 y1
// xn yn
// poly
// ...
// x0 y0
// x1 y1
// xn yn
// poly
void save_vertices_vectors(
    const char *output_file,
    std::vector<std::vector<const Point *>> &vertices_vectors,
    const char *mode = "w");

// save segment p1, p2 and vertices vectors in output file in given mode
void save_vertices_vectors_seg(std::vector<std::vector<const Point *>> vertices_vectors,
                     const Point &p1, const Point &p2,
                     const char *output_file="clip_error.txt");

//  save clipped vertices in output_file in given mode
// in format:
// x0 y0 
// x1 y1
// xn yn
// POLY i j
// x0 y0 
// x1 y1
// xn yn
// POLY i j
void save_clipped_vertices(
    const char *output_file,
    const std::map<std::pair<unsigned int, unsigned int>,
                   const std::vector<const Point *>> &i_j_to_clipped_vertices,
    const char *mode= "w");

// save clipped vertices and cell types in output_file in given mode
// in format:
// x0 y0 
// x1 y1
// xn yn
// poly 
// x0 y0 
// x1 y1
// xn yn
// poly 
// cell i j
// cell_type
// ...
void save_clipped_vertices_cell_type(const char *output_file,
                        const std::map<std::pair<unsigned int, unsigned int>,
                                       RasterCellInfo> &i_j_to_rcell_info,
                        const char *mode = "w");


#endif
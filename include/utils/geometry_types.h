#ifndef GEOMETRY_TYPES_H
#define GEOMETRY_TYPES_H

#include <iostream>
#include <limits>
#include <vector>

struct Point {
  double x, y;

  // Parameterized constructor
  Point(double x_, double y_);

  // Default constructor
  Point();

  // Copy assignment operator
  Point &operator=(const Point &other);

  // Return true if points are equal, else false
  bool equal_points(const Point &other, double epsilon) const;

  // returns point in str format
  // (x y)
  std::string to_str() const;

  // Destructor
  ~Point();
};

struct Polygon {
  std::vector<const Point *> vertices;
  int polygon_id;

  Point minCorner; // for mbr
  Point maxCorner; // for mbr

  // find the mbr of the vertices and set it to minCorner, maxCorner
  void findMBR();

  // save polygon vertices in output_file in given mode
  // x0 y0
  // x1 y1
  // xn yn
  void save_poly(const char *output_file, const char *mode = "w",
                 bool append_poly = false);

  // prints the vertices of the polygon
  void print_();

  bool intersects(const Polygon &other) const;
  bool point_inside(const Point &p) const;

  // return string containing the WKT representation of the polygon.
  std::string to_wkt();
  void save_vertices_to_csv(const char *output_file);

  // The signed area of the polygon is calculated using the formula
  // If the area is positive, the vertices are in a counterclockwise order
  // (true) if negative, they are in a clockwise order (false)
  bool is_ccw() const;

  // reverse the vertices order
  void make_cw();
};

// returns polygon area using shoelace algo
double polygon_area(const std::vector<const Point *> &vertices);

/* The result of this cross product can be used to determine the relative
   orientation of the points:
- If the result is positive, then p1p2 to p1p3 makes a counterclockwise turn.
- If the result is negative, then p1p2 to p1p3 makes a clockwise turn.
- If the result is zero, then the points are collinear.
return The cross product of vectors p1p2 and p1p3. */
double cross_product(const Point &p1, const Point &p2, const Point &p3);

// returns true of p1, p2, p3 are collinear, else false
bool are_collinear(const Point &p1, const Point &p2, const Point &p3,
                   double epsilon);

// calculates the x-coordinate where two lines, defined by
// the points (p1, p2) and (p3, p4), intersect. If the lines are parallel
// or collinear, the function returns `NaN` to indicate no valid intersection.
double x_intersect_(const Point &p1, const Point &p2, const Point &p3,
                    const Point &p4, double epsilon);

// calculates the y-coordinate where two lines, defined by
// the points (p1, p2) and (p3, p4), intersect. If the lines are parallel
// or collinear, the function returns `NaN` to indicate no valid intersection.
double y_intersect_(const Point &p1, const Point &p2, const Point &p3,
                    const Point &p4, double epsilon);

// calculates the point where two lines, defined by
// the points (p1, p2) and (p3, p4), intersect. If the lines are parallel
// or collinear, the function returns nullptr to indicate no valid intersection.
Point *p_intersect(const Point &p1, const Point &p2, const Point &p3,
                   const Point &p4, double epsilon = 1e-10);

// Function to save a list of polygons to a CSV file
// in wkt
void save_polygons_to_csv(std::vector<Polygon> &polygons,
                          const char *output_file);

// get total area of polygons in vertices_vectors
double
get_polygons_area(std::vector<std::vector<const Point *>> &vertices_vectors);

// print vertices of vec
void print_vec(const std::vector<const Point *> &vec);

// if p1, p2 ara parallel to x axis return 1
// else if on y axis return 2
// else 0
unsigned int checkParallelism(const Point &p1, const Point &p2);

#endif

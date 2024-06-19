#ifndef GEOMETRY_TYPES_H
#define GEOMETRY_TYPES_H

#include <iostream>
#include <limits>
#include <vector>

struct Point {
  static constexpr double EPSILON = 1e-7;
  double x, y;

  // Parameterized constructor
  Point(double x_, double y_);

  // Default constructor
  Point();

  // Copy assignment operator
  Point &operator=(const Point &other);

  // Equality operator
  bool operator==(const Point &other) const;

  std::string to_str() const;

  // Destructor
  ~Point();
};

struct Polygon {
  std::vector<const Point *> vertices;
  int polygon_id;

  Point minCorner;
  Point maxCorner;

  void findMBR();
  void save_poly(const char *output_file, const char *mode = "w");
  void print_();
  void save_poly_to_wkt(const char *output_file, const char *mode = "w");

  bool intersects(const Polygon &other) const;
  bool point_inside(const Point &p) const;
};

double polygon_area(const std::vector<const Point *> &vertices);

double cross_product_(const Point &p1, const Point &p2, const Point &p3);

bool are_collinear(const Point &p1, const Point &p2, const Point &p3);

double x_intersect_(const Point &p1, const Point &p2, const Point &p3,
                    const Point &p4);

double y_intersect_(const Point &p1, const Point &p2, const Point &p3,
                    const Point &p4);

Point *p_intersect(const Point &p1, const Point &p2, const Point &p3,
                   const Point &p4);

unsigned int equal_with_corner(const Point &p, const Point &p1, const Point &p2,
                               const Point &p3, const Point &p4,
                               bool debug = false);

bool collinear_with_any_edge(const Point &p, const Point &p1, const Point &p2,
                             const Point &p3, const Point &p4);

#endif

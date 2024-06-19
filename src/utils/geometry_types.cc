#include "../../include/utils/geometry_types.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

Point::Point(double x_, double y_) : x(x_), y(y_) {}

// Default constructor
Point::Point() : x(0), y(0) {}

// Copy assignment operator
Point &Point::operator=(const Point &other) {
  if (this == &other)
    return *this; // handle self-assignment
  x = other.x;
  y = other.y;
  return *this;
}

// Equality operator
bool Point::operator==(const Point &other) const {
  return (std::abs(x - other.x) < EPSILON && std::abs(y - other.y) < EPSILON);
}

std::string Point::to_str() const {
  std::string pt_to_str("(" + std::to_string(x) + ", " + std::to_string(y) +
                        "), ");
  return pt_to_str;
}

// // Destructor
Point::~Point() {}

void Polygon::findMBR() {

  if (vertices.empty()) {
    std::cerr << "Vertice vector in findMbr are empty\n";
    return;
  }

  minCorner = maxCorner = *vertices[0];
  for (const Point *vertex : vertices) {
    minCorner.x = std::min(minCorner.x, vertex->x);
    minCorner.y = std::min(minCorner.y, vertex->y);
    maxCorner.x = std::max(maxCorner.x, vertex->x);
    maxCorner.y = std::max(maxCorner.y, vertex->y);
  }
}

void Polygon::save_poly(const char *output_file, const char *mode) {
  FILE *fp = fopen(output_file, mode);
  for (int i = 0; i < vertices.size(); ++i) {
    fprintf(fp, "%lf %lf\n", vertices[i]->x, vertices[i]->y);
  }
  fclose(fp);
}

void Polygon::print_() {
  std::cout << std::fixed << std::setprecision(6);
  for (const Point *p : vertices) {
    std::cout << "(" << p->x << ", " << p->y << "), ";
  }
  std::cout << std::endl;
}

void Polygon::save_poly_to_wkt(const char *output_file, const char *mode) {
  FILE *fp = fopen(output_file, mode);
  fprintf(fp, "\"POLYGON ((");
  if (vertices.size()) {
    fprintf(fp, "%lf %lf", vertices[0]->x, vertices[0]->y);
  }
  for (int i = 1; i < vertices.size(); ++i) {
    fprintf(fp, ",%lf %lf", vertices[i]->x, vertices[i]->y);
  }
  fprintf(fp, "))\"\n");
  fclose(fp);
}

double cross_product_(const Point &p1, const Point &p2, const Point &p3) {
  return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

bool are_collinear(const Point &p1, const Point &p2, const Point &p3) {
  return std::abs(cross_product_(p1, p2, p3)) <
         1e-24; // Consider floating-point precision issues
}

double x_intersect_(const Point &p1, const Point &p2, const Point &p3,
                    const Point &p4) {
  const double EPSILON = 1e-24;
  double num = (p1.x * p2.y - p1.y * p2.x) * (p3.x - p4.x) -
               (p1.x - p2.x) * (p3.x * p4.y - p3.y * p4.x);
  double den = (p1.x - p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x - p4.x);
  if (std::fabs(den) < EPSILON) {
    std::cerr
        << "x_intersect_: Lines are parallel or collinear; no intersection.\n";
    std::cout << "new\n";
    std::cout << p1.x << " " << p1.y << " " << p2.x << " " << p2.y << " "
              << p3.x << " " << p3.y << " " << p4.x << " " << p4.y << std::endl;
    return std::numeric_limits<double>::
        quiet_NaN(); // Return Not-A-Number to indicate no valid intersection
  }
  return num / den;
}

double y_intersect_(const Point &p1, const Point &p2, const Point &p3,
                    const Point &p4) {
  const double EPSILON = 1e-24;
  double num = (p1.x * p2.y - p1.y * p2.x) * (p3.y - p4.y) -
               (p1.y - p2.y) * (p3.x * p4.y - p3.y * p4.x);
  double den = (p1.x - p2.x) * (p3.y - p4.y) - (p1.y - p2.y) * (p3.x - p4.x);
  if (std::fabs(den) < EPSILON) {
    std::cerr
        << "y_intersect_: Lines are parallel or collinear; no intersection.\n";
    return std::numeric_limits<double>::quiet_NaN();
  }
  return num / den;
}

Point *p_intersect(const Point &p1, const Point &p2, const Point &p3,
                   const Point &p4) {
  double x_inter = x_intersect_(p1, p2, p3, p4);
  double y_inter = y_intersect_(p1, p2, p3, p4);
  if (x_inter == std::numeric_limits<double>::quiet_NaN() ||
      y_inter == std::numeric_limits<double>::quiet_NaN()) {
    std::cerr << "Nan intersection\n";
    return nullptr;
  }
  Point *p_inter = new Point(x_inter, y_inter);
  return p_inter;
}

unsigned int equal_with_corner(const Point &p, const Point &p1, const Point &p2,
                               const Point &p3, const Point &p4, bool debug) {
  if (debug) {
    std::cout << "Point: (" << p.x << ", " << p.y << ") == Corner: (";
  }
  if (p == p1) {
    if (debug) {
      std::cout << p1.x << ", " << p1.y << ") -> 1\n";
    }
    return 1;
  } else if (p == p2) {
    if (debug) {
      std::cout << p2.x << ", " << p2.y << ") -> 2\n";
    }
    return 2;
  } else if (p == p3) {
    if (debug) {
      std::cout << p3.x << ", " << p3.y << ") -> 3\n";
    }
    return 3;
  } else if (p == p4) {
    if (debug) {
      std::cout << p4.x << ", " << p4.y << ") -> 4\n";
    }
    return 4;
  }
  return 0;
}

bool collinear_with_any_edge(const Point &p, const Point &p1, const Point &p2,
                             const Point &p3, const Point &p4) {
  return (are_collinear(p, p1, p2) || are_collinear(p, p2, p3) ||
          are_collinear(p, p3, p4) || are_collinear(p, p4, p1));
}

double polygon_area(const std::vector<const Point *> &points) {
  double area = 0.0;
  double c1 = 0.0; // Compensation for the first sum
  double c2 = 0.0; // Compensation for the second sum
  int n = points.size();

  for (int i = 0; i < n; i++) {
    int j = (i + 1) % n;

    // Apply Kahan Summation Algorithm for sum1
    double product1 = points[i]->x * points[j]->y;
    double y1 = product1 - c1;
    double t1 = area + y1;
    c1 = (t1 - area) - y1;
    area = t1;

    // Apply Kahan Summation Algorithm for sum2
    double product2 = points[j]->x * points[i]->y;
    double y2 = product2 - c2;
    double t2 =
        area - y2; // Subtracting here as we need to subtract sum2 from sum1
    c2 = (t2 - area) + y2; // Note the change in the compensation formula
    area = t2;
  }

  area = std::abs(area) / 2.0;
  return area;
}

int orientation(const Point &p, const Point &q, const Point &r) {
  double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
  if (val == 0)
    return 0;               // collinear
  return (val > 0) ? 1 : 2; // clock or counterclockwise
}

bool onSegment(const Point &p, const Point &q, const Point &r) {
  if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) &&
      q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y))
    return true;
  return false;
}

bool doIntersect(const Point &p1, const Point &q1, const Point &p2,
                 const Point &q2) {
  int o1 = orientation(p1, q1, p2);
  int o2 = orientation(p1, q1, q2);
  int o3 = orientation(p2, q2, p1);
  int o4 = orientation(p2, q2, q1);

  if (o1 != o2 && o3 != o4)
    return true;

  if (o1 == 0 && onSegment(p1, p2, q1))
    return true;
  if (o2 == 0 && onSegment(p1, q2, q1))
    return true;
  if (o3 == 0 && onSegment(p2, p1, q2))
    return true;
  if (o4 == 0 && onSegment(p2, q1, q2))
    return true;

  return false;
}

bool Polygon::point_inside(const Point &p) const {
  int n = vertices.size();
  if (n < 3)
    return false;

  Point extreme = {1e9, p.y};

  int count = 0, i = 0;
  do {
    int next = (i + 1) % n;

    if (doIntersect(*vertices[i], *vertices[next], p, extreme)) {
      if (orientation(*vertices[i], p, *vertices[next]) == 0)
        return onSegment(*vertices[i], p, *vertices[next]);

      count++;
    }
    i = next;
  } while (i != 0);

  return count & 1;
}

bool Polygon::intersects(const Polygon &poly2) const {
  for (size_t i = 0; i < vertices.size(); ++i) {
    for (size_t j = 0; j < vertices.size(); ++j) {
      size_t next_i = (i + 1) % vertices.size();
      size_t next_j = (j + 1) % poly2.vertices.size();
      if (doIntersect(*vertices[i], *vertices[next_i], *poly2.vertices[j],
                      *poly2.vertices[next_j])) {
        return true;
      }
    }
  }

  for (const Point *vertex : vertices) {
    if (poly2.point_inside(*vertex)) {
      return true;
    }
  }

  for (const Point *vertex : poly2.vertices) {
    if (point_inside(*vertex)) {
      return true;
    }
  }

  return false;
}
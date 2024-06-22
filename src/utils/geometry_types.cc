#include "../../include/utils/geometry_types.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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

// Destructor
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

std::string Polygon::to_wkt() {
    std::ostringstream oss;
    // Set fixed format and precision for output
    oss << std::fixed << std::setprecision(7);
    oss << "\"POLYGON ((";
    for (size_t i = 0; i < vertices.size(); ++i) {
        if (i != 0) {
            oss << ", ";
        }
        oss << vertices[i]->x << " " << vertices[i]->y;
    }
    oss << "))\"";
    return oss.str();
}

void Polygon::save_vertices_to_csv(const char* output_file) {
    std::ofstream file(output_file);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << output_file << std::endl;
        return;
    }

    // Set fixed format and precision for output
    file << std::fixed << std::setprecision(7);

    // Write vertices to CSV
    for (const auto& vertex : vertices) {
        file << vertex->x << ", " << vertex->y << "\n";
    }

    file.close();
    std::cout << "Vertices saved to " << output_file << std::endl;
}

void save_polygons_to_csv(std::vector<Polygon>& polygons, const char* output_file) {
    std::ofstream file(output_file);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << output_file << std::endl;
        return;
    }

    // Set fixed format and precision for output
    file << std::fixed << std::setprecision(7);

    for (auto& polygon : polygons) {
        file << polygon.to_wkt() << "\n";
    }

  file.close();
  std::cout << "Polygons saved to " << output_file << std::endl;
}

double cross_product(const Point &p1, const Point &p2, const Point &p3) {
  return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

bool are_collinear(const Point &p1, const Point &p2, const Point &p3) {
  return std::abs(cross_product(p1, p2, p3)) <
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

// shoelace algo
double polygon_area(const std::vector<const Point *> &points) {
  double area = 0.0;
  int n = points.size();

  for (int i = 0; i < n; i++) {
    int j = (i + 1) % n;
    area += points[i]->x * points[j]->y;
    area -= points[j]->x * points[i]->y;
  }

  // to handle any given polygon order
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
    for (size_t j = 0; j < poly2.vertices.size(); ++j) {
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

double
get_polygons_area(std::vector<std::vector<const Point *>> &vertices_vectors) {
  double polys_area = 0.0;
  for (std::vector<const Point *> &vertices : vertices_vectors) {
    polys_area += polygon_area(vertices);
    ;
  }
  return polys_area;
}

void print_vec(const std::vector<const Point *> &vec) {
  for (const Point *p : vec) {
    std::cout << std::fixed << std::setprecision(15);
    std::cout << "(" << p->x << ", " << p->y << "), ";
  }
  std::cout << std::endl;
}

unsigned int checkParallelism(const Point &p1, const Point &p2) {
  if (p1.y == p2.y) {
    return 1;
  } else if (p1.x == p2.x) {
    return 2;
  } else {
    return 0;
  }
}
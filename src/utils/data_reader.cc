#include "../../include/utils/data_reader.h"
#include <fstream>
#include <iostream>
#include <sstream>

std::vector<const Point *> parse_wkt(const std::string &wkt) {
  std::vector<const Point *> vertices;
  std::string clean_wkt = wkt.substr(11); // Skip the "POLYGON ((" part
  clean_wkt =
      clean_wkt.substr(0, clean_wkt.size() - 3); // Remove the trailing "))"
  std::istringstream iss(clean_wkt);
  std::string temp;
  while (std::getline(iss, temp, ',')) {
    std::istringstream coord_stream(temp);
    double x, y;
    coord_stream >> x >> y;
    const Point *p = new Point(x, y);
    vertices.push_back(p);
  }
  return vertices;
}

std::vector<Polygon> read_data_find_MBR(const std::string &filename,
                                        Point &gridMinCorner,
                                        Point &gridMaxCorner, int n) {
  std::vector<Polygon> polygons;
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return polygons;
  }

  int cnt = 0;
  std::string line;

  while (std::getline(file, line) && n--) {

    Polygon polygon;
    polygon.vertices = parse_wkt(line);
    // for (const Point* p : polygon.vertices) {
    //     std::cout << "(" << p->x << ", " << p->y << "), ";
    // }
    // std::cout << std::endl;
    polygon.findMBR();

    if (cnt == 0) {
      gridMinCorner.x = polygon.minCorner.x,
      gridMinCorner.y = polygon.minCorner.y,
      gridMaxCorner.x = polygon.maxCorner.x,
      gridMaxCorner.y = polygon.maxCorner.y;
    } else {
      if (polygon.minCorner.x < gridMinCorner.x) {
        gridMinCorner.x = polygon.minCorner.x;
      }
      if (polygon.minCorner.y < gridMinCorner.y) {
        gridMinCorner.y = polygon.minCorner.y;
      }
      if (polygon.maxCorner.x > gridMaxCorner.x) {
        gridMaxCorner.x = polygon.maxCorner.x;
      }
      if (polygon.maxCorner.y > gridMaxCorner.y) {
        gridMaxCorner.y = polygon.maxCorner.y;
      }
    }
    cnt++;
    polygons.push_back(polygon);
  }

  file.close();
  return polygons;
}
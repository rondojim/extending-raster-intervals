#include "../../include/utils/data_reader.h"
#include "../../include/utils/geometry_types.h"
#include <iostream>

int main() {
  Polygon polygon;
  polygon.vertices = parse_wkt(
      "POLYGON ((172.5775765 -43.553195, 172.5775242 -43.5533022, 172.57755 "
      "-43.5534639, 172.5777216 -43.5535759, 172.577919 -43.5536381, "
      "172.5780981 -43.5536499, 172.5782109 -43.5535634, 172.5782566 "
      "-43.5534857, 172.5782341 -43.5533595, 172.5781422 -43.5532773, "
      "172.5780564 -43.5532026, 172.5778161 -43.553128, 172.5776701 "
      "-43.553128, 172.5775765 -43.553195))");

  Point minCorner, maxCorner;
  polygon.findMBR();

  std::cout << "MBR of the polygon: ";
  std::cout << "(" << polygon.minCorner.x << ", " << polygon.minCorner.y
            << "), ";
  std::cout << "(" << polygon.maxCorner.x << ", " << polygon.maxCorner.y << ")"
            << std::endl;

  return 0;
}
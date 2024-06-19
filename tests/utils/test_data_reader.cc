#include "../../include/utils/data_reader.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main() {
  std::string filename = "../../../dataset_files/OSM_by_continent/O5OC";
  std::vector<Polygon> polygons = read_data(filename);
  std::cout << "Total polygons: " << polygons.size() << std::endl;
  int i = 0;
  for (const auto &polygon : polygons) {

    std::cout << "Polygon with " << polygon.vertices.size() << " vertices:\n";
    for (const auto &vertex : polygon.vertices) {
      std::cout << "(" << vertex.x << ", " << vertex.y << "), ";
    }
    std::cout << std::endl;
    break;
  }

  return 0;
}

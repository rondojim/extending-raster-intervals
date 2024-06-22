#include "../include/utils/pipeline.h"

int main() {
  unsigned int N = 12;
  std::string rhs_f_name("../../dataset_files/OSM_by_continent/O5OC");
  std::string lhs_f_name("../../dataset_files/OSM_by_continent/O6OC");
  // std::string lhs_f_name("../datasets/O5AF.txt");
  // std::string rhs_f_name("../datasets/O6AF.txt");
  std::pair<std::vector<Polygon>, std::vector<Polygon>> intersections;

  find_interesctions(N, lhs_f_name, rhs_f_name, intersections);

  return 0;
}
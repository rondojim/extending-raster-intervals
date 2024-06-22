#include "../include/utils/pipeline.h"

int main() {
  unsigned int N = 12;
  std::string rhs_f_name("../../dataset_files/archive/T1.csv");
  std::string lhs_f_name("../../dataset_files/archive/T2.csv");
  // std::string lhs_f_name("../datasets/O5AF.txt");
  // std::string rhs_f_name("../datasets/O6AF.txt");
  std::pair<std::vector<Polygon>, std::vector<Polygon>> intersections;

  find_interesctions(N, lhs_f_name, rhs_f_name, intersections);

  return 0;
}
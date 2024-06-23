#include "geometry_types.h"
#include <algorithm>
#include <chrono>
#include <set>
#include <vector>

void id_polygons(std::vector<Polygon> &polygons, bool lhs);

void get_preprocessed_polygons(std::vector<Polygon> &lhs,
                               std::vector<Polygon> &rhs,
                               const std::string &filename_lhs,
                               const std::string &filename_rhs,
                               Point &gridMinCorner, Point &gridMaxCorner,
                               int n);

void find_interesctions(
    unsigned int N, std::string lhs_f_name, std::string rhs_f_name,
    std::pair<std::vector<Polygon>, std::vector<Polygon>> &intersections,
    int n = -1);
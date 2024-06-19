#include "../../include/join_algos/combined_sweep.h"
#include "../../include/utils/data_reader.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/progressbar.hpp"
#include <algorithm>
#include <chrono>
#include <map>
#include <random>
#include <set>
#include <time.h>
#include <vector>

static void
id_polygons(std::vector<Polygon> &polygons, bool lhs,
            std::vector<std::pair<double, std::pair<int, int>>> &all_points) {
  int id = 1;
  for (auto &it : polygons) {
    it.findMBR();
    it.polygon_id = (lhs ? id : -id);
    id++;

    all_points.push_back({it.minCorner.x, {0, it.polygon_id}});
    all_points.push_back({it.minCorner.y, {1, it.polygon_id}});
    all_points.push_back({it.maxCorner.x, {2, it.polygon_id}});
    all_points.push_back({it.maxCorner.y, {3, it.polygon_id}});
  }
}

static std::pair<int, int>
coordinate_compression(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                       std::map<std::pair<int, int>, int> &compress) {
  std::vector<std::pair<double, std::pair<int, int>>> all_points;
  id_polygons(lhs, true, all_points);
  id_polygons(rhs, false, all_points);

  std::sort(all_points.begin(), all_points.end());
  printf("All points size: %d\n", all_points.size());

  double prev = all_points[0].first - 1;
  int idx = 0;
  int MAXX = 0, MAXY = 0;

  for (auto &it : all_points) {

    if (it.first > prev) {
      idx++;
    }
    if (it.second.first == 0 || it.second.first == 2) {
      MAXX = idx;
    } else {
      MAXY = idx;
    }
    compress[{it.second}] = idx;
    prev = it.first;
  }
  return {MAXX, MAXY};
}

static bool check_intersection(Polygon &a, Polygon &b) {
  if ((a.maxCorner.x < b.minCorner.x) || (a.minCorner.x > b.maxCorner.x) ||
      (a.maxCorner.y < b.minCorner.y) || (a.minCorner.y > b.maxCorner.y)) {
    return false;
  }
  return true;
}

static bool
test_results(std::pair<std::vector<Polygon>, std::vector<Polygon>> &result) {
  progressbar bar(result.first.size());

  for (auto &it : result.first) {
    bool ok = false;
    for (auto &it2 : result.second) {

      if (check_intersection(it, it2)) {
        ok = true;
        break;
      }
    }
    if (!ok) {
      return false;
    }
    bar.update();
  }
  return true;
}

// void generate_polygons(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs)
// {
//     Polygon l1, l2, l3, l4;
//     Polygon r1, r2, r3;
//     l1.vertices = {{1.0, 3.0}, {5.0, 3.0}, {5.0, 7.0}, {1.0, 7.0}};
//     l2.vertices = {{2.0, 4.0}, {4.0, 4.0}, {4.0, 6.0}, {2.0, 6.0}};
//     l3.vertices = {{4.0, 2.0}, {6.0, 2.0}, {6.0, 4.0}, {4.0, 4.0}};
//     l4.vertices = {{6.0, 7.0}, {9.0, 7.0}, {9.0, 11.0}, {6.0, 11.0}};

//     r1.vertices = {{2.0, 7.0}, {4.0, 7.0}, {4.0, 8.0}, {2.0, 8.0}};
//     r2.vertices = {{7.0, 9.0}, {8.0, 9.0}, {8.0, 10.0}, {7.0, 10.0}};
//     r3.vertices = {{4.0, 9.0}, {5.0, 9.0}, {5.0, 10.0}, {4.0, 10.0}};

//     lhs = {l1, l2, l3, l4};
//     rhs = {r1, r2, r3};
// }

int main() {
  auto start = std::chrono::high_resolution_clock::now();
  std::vector<Polygon> lhs_polygons =
      read_data_find_MBR("/Users/dimronto/Desktop/DSIT/Spring/Database "
                         "Systems/project/dataset_files/OSM_by_continent/O5OC");
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;

  printf("Read first set in %.1f seconds\n", duration.count());

  start = std::chrono::high_resolution_clock::now();
  std::vector<Polygon> rhs_polygons =
      read_data("/Users/dimronto/Desktop/DSIT/Spring/Database "
                "Systems/project/dataset_files/OSM_by_continent/O6OC");
  end = std::chrono::high_resolution_clock::now();
  duration = end - start;

  printf("Read second set in %.1f seconds\n", duration.count());

  std::pair<std::vector<Polygon>, std::vector<Polygon>> final_result;

  mbr_combined_sweep(lhs_polygons, rhs_polygons, final_result);

  printf("%d\n", final_result.first.size() + final_result.second.size());

  return 0;
}
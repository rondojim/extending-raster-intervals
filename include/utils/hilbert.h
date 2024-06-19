#ifndef HILBERT_H
#define HILBERT_H

#include "../join_algos/inter_filter.h"
#include <map>
#include <vector>

void hilbert_sort(int N, std::map<std::pair<int, int>, int> &result);
void saveMapToFile(const std::map<std::pair<int, int>, int> &myMap,
                   const std::string &filename);
void loadMapFromFile(std::map<std::pair<int, int>, int> &myMap,
                     const std::string &filename);

#endif
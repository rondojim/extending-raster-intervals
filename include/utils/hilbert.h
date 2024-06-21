#ifndef HILBERT_H
#define HILBERT_H

#include "../join_algos/inter_filter.h"
#include <map>
#include <string>
#include <vector>

// Function to perform Hilbert sort on a collection of coordinates
// @param N: Size of the grid
// @param result: Map to store the Hilbert curve index of each coordinate pair
void hilbert_sort(int N, std::map<std::pair<int, int>, int> &result);

// Function to save a map to a file
// @param myMap: The map to be saved
// @param filename: The name of the file to save the map to
void saveMapToFile(const std::map<std::pair<int, int>, int> &myMap,
                   const std::string &filename);

// Function to load a map from a file
// @param myMap: The map to be loaded
// @param filename: The name of the file to load the map from
void loadMapFromFile(std::map<std::pair<int, int>, int> &myMap,
                     const std::string &filename);

#endif // HILBERT_H

#ifndef DATA_READER_H
#define DATA_READER_H

#include "geometry_types.h"
#include <string>
#include <vector>

std::vector<const Point *> parse_wkt(const std::string &wkt);
std::vector<Polygon> read_data_find_MBR(const std::string &filename,
                                        Point &gridMinCorner,
                                        Point &gridMaxCorner, int n = -1);

#endif // DATA_READER_H
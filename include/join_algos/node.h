#ifndef NODE_H
#define NODE_H

#include <set>

struct NodeY {
  bool active_polygons;
  bool whole;
  NodeY *left, *right;

  NodeY();
};

struct NodeX {
  NodeY *rootY;
  NodeY *rootY_whole;
  NodeX *left, *right;

  NodeX();
};

#endif // NODE_H
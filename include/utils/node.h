#ifndef NODE_H
#define NODE_H

#include <set>

// Structure representing a node in the Y segment tree
struct NodeY {
  bool active_polygons; // Flag indicating if the node has active polygons
  bool whole; // Flag indicating if the node represents a whole interval
  NodeY *left, *right; // Pointers to the left and right children

  // Constructor for initializing a NodeY
  NodeY();
};

// Structure representing a node in the X segment tree
struct NodeX {
  NodeY *rootY;        // Root of the Y segment tree
  NodeY *rootY_whole;  // Root of the Y segment tree for whole intervals
  NodeX *left, *right; // Pointers to the left and right children

  // Constructor for initializing a NodeX
  NodeX();
};

#endif // NODE_H

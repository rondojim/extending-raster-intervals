#include "../../include/mbr_algos/mbr_2d_segt.h"
#include "../../include/utils/compressing.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/node.h"
#include "../../include/utils/progressbar.hpp"
#include <algorithm>
#include <chrono>
#include <limits>
#include <map>
#include <set>
#include <vector>

// Function to update the Y segment tree with an interval
static void updateY(NodeY *&rootY, int lo, int hi, int min_y, int max_y) {
  rootY->active_polygons = true;
  if (lo == min_y && hi == max_y) {
    rootY->whole = true;
    return;
  } else {
    int mid = (lo + hi) / 2;

    if (max_y <= mid) {
      // Update the left child if necessary
      if (!rootY->left) {
        rootY->left = new NodeY();
      }
      updateY(rootY->left, lo, mid, min_y, max_y);
    } else if (min_y > mid) {
      // Update the right child if necessary
      if (!rootY->right) {
        rootY->right = new NodeY();
      }
      updateY(rootY->right, mid + 1, hi, min_y, max_y);
    } else {
      // Update both children if necessary
      if (!rootY->left) {
        rootY->left = new NodeY();
      }
      if (!rootY->right) {
        rootY->right = new NodeY();
      }
      updateY(rootY->left, lo, mid, min_y, mid);
      updateY(rootY->right, mid + 1, hi, mid + 1, max_y);
    }
  }
}

// Function to query the Y segment tree for an interval
static bool queryY(NodeY *&rootY, int lo, int hi, int min_y, int max_y) {
  if (rootY->whole) {
    return true;
  }
  if (lo == min_y && hi == max_y) {
    return rootY->active_polygons;
  } else {
    int mid = (lo + hi) / 2;

    if (max_y <= mid) {
      if (!rootY->left) {
        return false;
      } else {
        return queryY(rootY->left, lo, mid, min_y, max_y);
      }
    } else if (min_y > mid) {
      if (!rootY->right) {
        return false;
      } else {
        return queryY(rootY->right, mid + 1, hi, min_y, max_y);
      }
    } else {
      bool result = false;
      if (rootY->left) {
        result |= queryY(rootY->left, lo, mid, min_y, mid);
      }
      if (!result && rootY->right) {
        result |= queryY(rootY->right, mid + 1, hi, mid + 1, max_y);
      }
      return result;
    }
  }
}

// Function to update the X segment tree with an interval, which involves
// updating the Y segment tree as well
static void updateX(NodeX *&rootX, int lo, int hi, int min_x, int max_x,
                    int min_y, int max_y, int MAXC) {
  // Update the Y segment tree associated with the current X node
  updateY(rootX->rootY, 1, MAXC, min_y, max_y);
  if (lo == min_x && hi == max_x) {
    if (!rootX->rootY_whole) {
      rootX->rootY_whole = new NodeY();
    }
    updateY(rootX->rootY_whole, 1, MAXC, min_y, max_y);
    return;
  } else {
    int mid = (lo + hi) / 2;
    if (max_x <= mid) {
      // Update the left child if necessary
      if (!rootX->left) {
        rootX->left = new NodeX();
      }
      updateX(rootX->left, lo, mid, min_x, max_x, min_y, max_y, MAXC);
    } else if (min_x > mid) {
      // Update the right child if necessary
      if (!rootX->right) {
        rootX->right = new NodeX();
      }
      updateX(rootX->right, mid + 1, hi, min_x, max_x, min_y, max_y, MAXC);
    } else {
      // Update both children if necessary
      if (!rootX->left) {
        rootX->left = new NodeX();
      }
      if (!rootX->right) {
        rootX->right = new NodeX();
      }
      updateX(rootX->left, lo, mid, min_x, mid, min_y, max_y, MAXC);
      updateX(rootX->right, mid + 1, hi, mid + 1, max_x, min_y, max_y, MAXC);
    }
  }
}

// Function to query the X segment tree for an interval, which involves querying
// the Y segment tree as well
static bool queryX(NodeX *&rootX, int lo, int hi, int min_x, int max_x,
                   int min_y, int max_y, int MAXC) {

  if (rootX->rootY_whole) {
    if (queryY(rootX->rootY_whole, 1, MAXC, min_y, max_y)) {
      return true;
    }
  }

  if (lo == min_x && hi == max_x) {
    return queryY(rootX->rootY, 1, MAXC, min_y, max_y);
  } else {
    int mid = (lo + hi) / 2;

    if (max_x <= mid) {
      if (!rootX->left) {
        return false;
      } else {
        return queryX(rootX->left, lo, mid, min_x, max_x, min_y, max_y, MAXC);
      }
    } else if (min_x > mid) {
      if (!rootX->right) {
        return false;
      } else {
        return queryX(rootX->right, mid + 1, hi, min_x, max_x, min_y, max_y,
                      MAXC);
      }
    } else {
      bool result = false;
      if (rootX->left) {
        result |= queryX(rootX->left, lo, mid, min_x, mid, min_y, max_y, MAXC);
      }
      if (!result && rootX->right) {
        result |= queryX(rootX->right, mid + 1, hi, mid + 1, max_x, min_y,
                         max_y, MAXC);
      }
      return result;
    }
  }
}

// Function to filter polygons using the 2D segment tree
static void
filter_polygons(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                std::map<std::pair<int, int>, int> &compress, int MAXX,
                int MAXY,
                std::pair<std::vector<Polygon>, std::vector<Polygon>> &result) {
  // Initialize the root of the X segment tree
  NodeX *rootX = new NodeX();

  // Update the segment tree with the left-hand side (lhs) polygons
  for (auto &it : lhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];

    // Update the X segment tree with the MBR coordinates of the polygon
    updateX(rootX, 1, MAXX, min_x, max_x, min_y, max_y, MAXY);
  }

  // Query the segment tree with the right-hand side (rhs) polygons
  for (auto &it : rhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];

    // Check if the polygon intersects with any polygons in the segment tree
    if (queryX(rootX, 1, MAXX, min_x, max_x, min_y, max_y, MAXY)) {
      // If intersection is found, add the polygon to the result set
      result.second.push_back(it);
    }
  }

  // Initialize a new root for the X segment tree
  NodeX *rootX2 = new NodeX();

  // Update the new segment tree with the right-hand side (rhs) polygons
  for (auto &it : rhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];

    // Update the X segment tree with the MBR coordinates of the polygon
    updateX(rootX2, 1, MAXX, min_x, max_x, min_y, max_y, MAXY);
  }

  // Query the new segment tree with the left-hand side (lhs) polygons
  for (auto &it : lhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];

    // Check if the polygon intersects with any polygons in the segment tree
    if (queryX(rootX2, 1, MAXX, min_x, max_x, min_y, max_y, MAXY)) {
      // If intersection is found, add the polygon to the result set
      result.first.push_back(it);
    }
  }
}

// Function to perform MBR filtering using a 2D segment tree
void mbr_2D_segt(
    std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
    std::pair<std::vector<Polygon>, std::vector<Polygon>> &result) {
  std::map<std::pair<int, int>, int> compress; // Map for coordinate compression

  // Perform coordinate compression and get the maximum compressed coordinates
  std::pair<int, int> coords = coordinate_compression(lhs, rhs, compress);
  int MAXX = coords.first, MAXY = coords.second;
  filter_polygons(lhs, rhs, compress, MAXX, MAXY, result);
}

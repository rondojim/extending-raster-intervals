#include "../../include/join_algos/mbr_suboptimal_2d_segt.h"
#include "../../include/join_algos/node.h"
#include "../../include/utils/geometry_types.h"
#include "../../include/utils/progressbar.hpp"
#include <algorithm>
#include <chrono>
#include <limits>
#include <map>
#include <set>
#include <vector>

static void updateY(NodeY *&rootY, int lo, int hi, int min_y, int max_y) {
  rootY->active_polygons = true;
  if (lo == min_y && hi == max_y) {
    rootY->whole = true;
    return;
  } else {
    int mid = (lo + hi) / 2;

    if (max_y <= mid) {
      if (!rootY->left) {
        rootY->left = new NodeY();
      }
      updateY(rootY->left, lo, mid, min_y, max_y);
    } else if (min_y > mid) {
      if (!rootY->right) {
        rootY->right = new NodeY();
      }
      updateY(rootY->right, mid + 1, hi, min_y, max_y);
    } else {
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

static void updateX(NodeX *&rootX, int lo, int hi, int min_x, int max_x,
                    int min_y, int max_y, int MAXC) {
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
      if (!rootX->left) {
        rootX->left = new NodeX();
      }
      updateX(rootX->left, lo, mid, min_x, max_x, min_y, max_y, MAXC);
    } else if (min_x > mid) {
      if (!rootX->right) {
        rootX->right = new NodeX();
      }
      updateX(rootX->right, mid + 1, hi, min_x, max_x, min_y, max_y, MAXC);
    } else {
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

static void
filter_polygons(std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
                std::map<std::pair<int, int>, int> &compress, int MAXX,
                int MAXY,
                std::pair<std::vector<Polygon>, std::vector<Polygon>> &result) {
  NodeX *rootX = new NodeX();

  printf("\nUpdating lhs into segment tree\n");
  int iters = lhs.size();
  progressbar bar(iters);
  for (auto &it : lhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];

    updateX(rootX, 1, MAXX, min_x, max_x, min_y, max_y, MAXY);
    bar.update();
  }

  printf("\nQuerying segment tree for rhs\n");
  iters = rhs.size();

  progressbar bar2(iters);
  for (auto &it : rhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];

    if (queryX(rootX, 1, MAXX, min_x, max_x, min_y, max_y, MAXY)) {
      result.second.push_back(it);
    }
    bar2.update();
  }

  NodeX *rootX2 = new NodeX();

  printf("\nUpdating rhs into new segment tree\n");
  iters = rhs.size();
  progressbar bar3(iters);
  for (auto &it : rhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];

    updateX(rootX2, 1, MAXX, min_x, max_x, min_y, max_y, MAXY);
    bar3.update();
  }

  printf("\nQuerying new segment tree for lhs\n");
  iters = lhs.size();
  progressbar bar4(iters);
  for (auto &it : lhs) {
    int min_x = compress[{0, it.polygon_id}];
    int max_x = compress[{2, it.polygon_id}];
    int min_y = compress[{1, it.polygon_id}];
    int max_y = compress[{3, it.polygon_id}];

    if (queryX(rootX2, 1, MAXX, min_x, max_x, min_y, max_y, MAXY)) {
      result.first.push_back(it);
    }
    bar4.update();
  }
}

void mbr_suboptimal_2D(
    std::vector<Polygon> &lhs, std::vector<Polygon> &rhs,
    std::map<std::pair<int, int>, int> &compress, int MAXX, int MAXY,
    std::pair<std::vector<Polygon>, std::vector<Polygon>> &result) {
  filter_polygons(lhs, rhs, compress, MAXX, MAXY, result);
}
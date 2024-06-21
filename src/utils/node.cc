#include "../../include/utils/node.h"
#include <cassert>

NodeY::NodeY()
    : active_polygons(false), whole(false), left(nullptr), right(nullptr) {}

NodeX::NodeX() : left(nullptr), right(nullptr), rootY_whole(nullptr) {
  rootY = new NodeY();
}
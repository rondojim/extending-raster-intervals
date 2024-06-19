#ifndef INTERVAL_H
#define INTERVAL_H

struct Interval {
  int x;
  int min_y;
  int max_y;
  int id;
  bool is_start;

  Interval(int x, int min_y, int max_y, int id, bool is_start)
      : x(x), min_y(min_y), max_y(max_y), id(id), is_start(is_start) {}
};

struct IntervalComparePositiveIdPriority {
  bool operator()(const Interval &lhs, const Interval &rhs) const {
    if (lhs.x != rhs.x) {
      return lhs.x < rhs.x;
    }
    if ((lhs.id > 0) == (rhs.id > 0)) {
      return lhs.id < rhs.id;
    }
    if (lhs.id > 0) {
      if (!lhs.is_start && !rhs.is_start)
        return true;
      if (!lhs.is_start && rhs.is_start)
        return false;
      if (lhs.is_start && !rhs.is_start)
        return true;
      if (lhs.is_start && rhs.is_start)
        return false;
    }
    if (rhs.id > 0) {
      if (!lhs.is_start && !rhs.is_start)
        return false;
      if (!lhs.is_start && rhs.is_start)
        return true;
      if (lhs.is_start && !rhs.is_start)
        return false;
      if (lhs.is_start && rhs.is_start)
        return true;
    }
    return lhs.id < rhs.id;
  }
};

struct IntervalCompareNegativeIdPriority {
  bool operator()(const Interval &lhs, const Interval &rhs) const {
    if (lhs.x != rhs.x) {
      return lhs.x < rhs.x;
    }
    if ((lhs.id > 0) == (rhs.id > 0)) {
      return lhs.id < rhs.id;
    }
    if (lhs.id < 0) {
      if (!lhs.is_start && !rhs.is_start)
        return true;
      if (!lhs.is_start && rhs.is_start)
        return false;
      if (lhs.is_start && !rhs.is_start)
        return true;
      if (lhs.is_start && rhs.is_start)
        return false;
    }
    if (rhs.id < 0) {
      if (!lhs.is_start && !rhs.is_start)
        return false;
      if (!lhs.is_start && rhs.is_start)
        return true;
      if (lhs.is_start && !rhs.is_start)
        return false;
      if (lhs.is_start && rhs.is_start)
        return true;
    }
    return lhs.id < rhs.id;
  }
};

#endif // INTERVAL_H
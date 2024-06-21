#ifndef INTERVAL_H
#define INTERVAL_H

// Interval structure to represent the segments of Minimum Bounding Rectangles
// (MBRs)
struct Interval {
  int x;         // X coordinate of the segment
  int min_y;     // Minimum Y coordinate of the segment
  int max_y;     // Maximum Y coordinate of the segment
  int id;        // Identifier for the interval (polygon ID)
  bool is_start; // Flag to indicate if it is the starting segment of the
                 // interval

  // Constructor to initialize the Interval structure
  Interval(int x, int min_y, int max_y, int id, bool is_start)
      : x(x), min_y(min_y), max_y(max_y), id(id), is_start(is_start) {}
};

// Comparator to prioritize intervals with positive IDs
struct IntervalComparePositiveIdPriority {
  // Operator to compare two intervals
  bool operator()(const Interval &lhs, const Interval &rhs) const {
    // Compare the X coordinates of the intervals
    if (lhs.x != rhs.x) {
      return lhs.x < rhs.x;
    }
    // If they belong in the same collection (positive or negative IDs),
    // is_start does not matter
    if ((lhs.id > 0) == (rhs.id > 0)) {
      return lhs.id < rhs.id;
    }
    // If lhs has a positive ID and the other has a negative
    if (lhs.id > 0) {
      // If both intervals are not starting segments, prioritize lhs
      if (!lhs.is_start && !rhs.is_start)
        return true;
      // If the left interval is not a starting segment and the right one is,
      // prioritize the right one
      if (!lhs.is_start && rhs.is_start)
        return false;
      // If the left interval is a starting segment and the right one is not,
      // prioritize the left one
      if (lhs.is_start && !rhs.is_start)
        return true;
      // If both intervals are starting segments, prioritize the one with the
      // negative ID
      if (lhs.is_start && rhs.is_start)
        return false;
    }
    // If lhs has a negative ID and the other has a positive
    if (rhs.id > 0) {
      // If both intervals are not starting segments, prioritize lhs
      if (!lhs.is_start && !rhs.is_start)
        return false;
      // If the left interval is not a starting segment and the right one is,
      // prioritize the left one
      if (!lhs.is_start && rhs.is_start)
        return true;
      // If the left interval is a starting segment and the right one is not,
      // prioritize the right one
      if (lhs.is_start && !rhs.is_start)
        return false;
      // If both intervals are starting segments, prioritize left one
      if (lhs.is_start && rhs.is_start)
        return true;
    }
    return lhs.id < rhs.id;
  }
};

// Comparator to prioritize intervals with negative IDs
struct IntervalCompareNegativeIdPriority {
  // Operator to compare two intervals
  bool operator()(const Interval &lhs, const Interval &rhs) const {
    // Compare the X coordinates of the intervals
    if (lhs.x != rhs.x) {
      return lhs.x < rhs.x;
    }
    // If they belong in the same collection (positive or negative IDs),
    // is_start does not matter
    if ((lhs.id > 0) == (rhs.id > 0)) {
      return lhs.id < rhs.id;
    }
    // If lhs has a positive ID and the other has a negative
    if (lhs.id < 0) {
      // If both intervals are not starting segments, prioritize lhs
      if (!lhs.is_start && !rhs.is_start)
        return true;
      // If the left interval is not a starting segment and the right one is,
      // prioritize the right one
      if (!lhs.is_start && rhs.is_start)
        return false;
      // If the left interval is a starting segment and the right one is not,
      // prioritize the left one
      if (lhs.is_start && !rhs.is_start)
        return true;
      // If both intervals are starting segments, prioritize the one with the
      // positive ID
      if (lhs.is_start && rhs.is_start)
        return false;
    }
    // If rhs has a negative ID and the other has a positive
    if (rhs.id < 0) {
      // If both intervals are not starting segments, prioritize rhs
      if (!lhs.is_start && !rhs.is_start)
        return false;
      // If the left interval is not a starting segment and the right one is,
      // prioritize the left one
      if (!lhs.is_start && rhs.is_start)
        return true;
      // If the left interval is a starting segment and the right one is not,
      // prioritize the right one
      if (lhs.is_start && !rhs.is_start)
        return false;
      // If both intervals are starting segments, prioritize right one
      if (lhs.is_start && rhs.is_start)
        return true;
    }
    return lhs.id < rhs.id;
  }
};

#endif // INTERVAL_H

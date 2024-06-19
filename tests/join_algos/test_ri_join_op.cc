#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <stdio.h>

int main() {
  // [5, 7]

  // create lhs bitmask with the bits 100101110
  boost::dynamic_bitset<> lhs_bitmask(9, 0b100101110);

  std::cout << lhs_bitmask << std::endl;

  int start = 7, end = 7;
  int ls = 5, le = 7;

  lhs_bitmask >>= (start - ls) * 3;
  lhs_bitmask.resize(lhs_bitmask.size() - (start - ls) * 3);

  std::cout << "After first truncation: " << lhs_bitmask << std::endl;
  std::cout << "Size: " << lhs_bitmask.size() << std::endl;

  lhs_bitmask <<= (le - end) * 3;
  lhs_bitmask >>= (le - end) * 3;
  lhs_bitmask.resize(lhs_bitmask.size() - (le - end) * 3);

  std::cout << "After second truncation: " << lhs_bitmask << std::endl;
  std::cout << "Size: " << lhs_bitmask.size() << std::endl;

  return 0;
}
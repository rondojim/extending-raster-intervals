#include "../include/utils/hilbert.h"
#include <algorithm>
#include <map>
#include <stdio.h>

int main() {
  for (int N = 2; N <= 2; ++N) {
    std::map<std::pair<int, int>, int> result;
    hilbert_sort(N, result);

    char filename[60];
    snprintf(filename, sizeof(filename), "../hilbert_maps/hilbert_%d.txt", N);

    printf("Saving to file %s\n", filename);

    saveMapToFile(result, filename);
  }

  return 0;
}
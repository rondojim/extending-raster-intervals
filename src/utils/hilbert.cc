#include <fstream>
#include <iostream>
#include <map>
#include <vector>

int dir[] = {0, 1, 0, -1};

static std::pair<std::pair<int, int>, std::pair<int, int>>
hilbert(int level, int dx, int dy, int angle, int x, int y,
        std::vector<std::pair<int, int>> &coords) {
  if (level == 0) {
    coords.push_back({x, y});
    return {{x, y}, {dx, dy}};
  }

  dx = (dx + angle + 4) % 4;
  dy = (dy + angle + 4) % 4;

  auto state = hilbert(level - 1, dx, dy, -angle, x, y, coords);
  x = state.first.first;
  y = state.first.second;
  dx = state.second.first;
  dy = state.second.second;

  x += dir[dx];
  y += dir[dy];

  dx = (dx - angle + 4) % 4;
  dy = (dy - angle + 4) % 4;

  state = hilbert(level - 1, dx, dy, angle, x, y, coords);

  x = state.first.first;
  y = state.first.second;
  dx = state.second.first;
  dy = state.second.second;

  x += dir[dx];
  y += dir[dy];

  state = hilbert(level - 1, dx, dy, angle, x, y, coords);
  x = state.first.first;
  y = state.first.second;
  dx = state.second.first;
  dy = state.second.second;

  dx = (dx - angle + 4) % 4;
  dy = (dy - angle + 4) % 4;
  x += dir[dx];
  y += dir[dy];

  state = hilbert(level - 1, dx, dy, -angle, x, y, coords);
  x = state.first.first;
  y = state.first.second;
  dx = state.second.first;
  dy = state.second.second;

  dx = (dx + angle + 4) % 4;
  dy = (dy + angle + 4) % 4;

  return {{x, y}, {dx, dy}};
}

void hilbert_sort(int N, std::map<std::pair<int, int>, int> &result) {
  std::vector<std::pair<int, int>> coords;
  hilbert(N, 1, 0, -1, 0, (1 << N) - 1, coords);

  int id = 1;
  for (auto it : coords) {
    result[it] = id++;
  }
}

void saveMapToFile(const std::map<std::pair<int, int>, int> &myMap,
                   const std::string &filename) {
  std::ofstream ofs(filename);
  if (ofs.is_open()) {
    for (const auto &item : myMap) {
      ofs << item.first.first << " " << item.first.second << " " << item.second
          << "\n";
    }
    ofs.close();
  } else {
    std::cerr << "Unable to open file for writing: " << filename << "\n";
  }
}

void loadMapFromFile(std::map<std::pair<int, int>, int> &myMap,
                     const std::string &filename) {

  std::ifstream ifs(filename);
  if (ifs.is_open()) {
    int first, second, value;
    while (ifs >> first >> second >> value) {
      myMap[std::make_pair(first, second)] = value;
    }
    ifs.close();
  } else {
    std::cerr << "Unable to open file for reading: " << filename << "\n";
  }
}
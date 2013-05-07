#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

#include <pari/pari.h>

#include "liboptarith/s128_c.h"

static inline uint64_t current_nanos(void) {
  struct timespec res;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &res);
  return (res.tv_sec * 1000000000ULL) + res.tv_nsec;
}

uint64_t avg(const vector<uint64_t>& xs) {
  if (xs.size() == 0)
    return 0;
  s128 sum = 0;
  for (uint64_t x : xs) {
    sum += x;
  }
  return (sum / xs.size()).to_u64();
}

vector<GEN> readfile(const string& filename) {
  vector<GEN> res;
  ifstream f(filename);
  while (!f.eof()) {
    string s;
    f >> s;
    if (s != "") {
      res.push_back(strtoi(s.c_str()));
    }
  }
  return res;
}

int main(int argc, char** argv) {
  pari_init(1<<30, 10000000);
  vector<GEN> xs = readfile(argv[1]);
  vector<uint64_t> times;
  for (auto x : xs) {
    uint64_t start = current_nanos();
    Z_factor(x);
    times.push_back(current_nanos() - start);
  }
  
  sort(times.begin(), times.end());
  cout << "min: " << times[0] << ' ';
  cout << "max: " << times[times.size()-1] << ' ';
  cout << "median: " << times[times.size()/2] << ' ';
  cout << "avg: " << avg(times) << endl;
  return 0;
}

/**
 * Batch process a file on composites for timing purposes.
 * Farms the factoring out to an implementation of SQUFOF taken from
 * Pari/GP, and another modified to use 128-bit arithmetic.
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

#include "liboptarith/s128_c.h"
#include "squfof-parigp/squfof64.h"
#include "squfof-parigp/squfof128.h"

using namespace std;

/// Gives the time from system on in nanoseconds
static inline uint64_t current_nanos(void) {
#ifdef __linux__
  struct timespec res;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &res);
  return (res.tv_sec * 1000000000ULL) + res.tv_nsec;
#else
  struct timeval tv;
  gettimeofday(&tv, 0);
  return ((uint64_t)tv.tv_sec * 1000000ULL + (uint64_t)tv.tv_usec) * 1000;
#endif
}

s128 squfof(s128 N) {
  if (N.msb() >= 59) return squfof128(N);
  else squfof64(N.to_u64());
}

void read_file_or_die(vector<s128>* xs, const string& filename) {
  xs->clear();
  ifstream f(filename.c_str());
  if (!f.is_open()) {
    cerr << "Couldn't find file " << filename << endl;
    exit(-1);
  }
  while (!f.eof()) {
    string tmp;
    f >> tmp;
    if (!f.fail()) {
      xs->push_back(s128(tmp.c_str()));
    }
  }
  f.close();
}

void time_factoring(vector<uint64_t>* times, const vector<s128>& xs) {
  times->clear();
  for (auto x : xs) {
    uint64_t start = current_nanos();
    s128 p = squfof(x);
    if (p != 0 && p != 1 && p != x) {
      uint64_t total_time = current_nanos() - start;
      times->push_back(total_time);
    }
  }
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


int main(int argc, char** argv) {
  if (argc != 2) {
    cout << "Usage: " << argv[0] << " <composites.txt>" << endl;
    cout << endl;
    return 0;
  }

  vector<s128> composites;
  read_file_or_die(&composites, argv[1]);
  cout << "Factoring " << composites.size() << " composites." << endl;
  vector<uint64_t> times;
  time_factoring(&times, composites);
  cout << "Factored " << times.size() << " composites." << endl;

  // Compute stats.
  sort(times.begin(), times.end());
  cout << "min: " << times[0] << ' ';
  cout << "max: " << times[times.size()-1] << ' ';
  cout << "median: " << times[times.size()/2] << ' ';
  cout << "avg: " << avg(times) << endl;

  return 0;
}




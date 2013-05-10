#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <stdint.h>
#include <time.h>

#include "liboptarith/s128_c.h"

using namespace std;

/// These must be provided.
string out_filename();
s128 factor(s128 N);

string in_filename(int bits) {
  return "/home/max/bitbucket/sspar/comparisons/composites/composites-" + to_string(bits) + ".txt";
}

/// Gives the time from system on in nanoseconds
static inline uint64_t current_nanos(void) {
  struct timespec res;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &res);
  return (res.tv_sec * 1000000000ULL) + res.tv_nsec;
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

void time_factoring(vector<double>* times, const vector<s128>& xs) {
  times->clear();
  for (auto x : xs) {
    uint64_t start = current_nanos();
    factor(x);
    uint64_t total_time = current_nanos() - start;
    times->push_back(static_cast<double>(total_time) / 1000);
  }
}

double avg(const vector<double>& xs) {
  if (xs.size() == 0)
    return 0;
  double sum = 0;
  for (double x : xs) {
    sum += x;
  }
  return sum / xs.size();
}

void dobits(const int bits) {
  cout << "Stats for " << bits << " bits" << endl;
  string filename = in_filename(bits);

  vector<s128> composites;
  read_file_or_die(&composites, filename);
  cout << "Factoring " << composites.size() << " semiprimes." << endl;
  vector<double> times;
  time_factoring(&times, composites);

  // Compute stats.
  sort(times.begin(), times.end());
  double m0, m1, m2, m3, m4;
  size_t s = times.size();
  m0 = times[0];
  m1 = times[s/4];
  m2 = times[s/2];
  m3 = times[3*s/4];
  m4 = times[s-1];
  cout << fixed << setprecision(5);
  cout << "m0: " << m0 << ' ';
  cout << "m1: " << m1 << ' ';
  cout << "m2: " << m2 << ' ';
  cout << "m3: " << m3 << ' ';
  cout << "m4: " << m4 << ' ';
  cout << "avg: " << avg(times) << endl;

  // append to output file
  ofstream out;
  out.open(out_filename().c_str(), ios::out | ios::app);
  out << fixed << setprecision(5);
  out << bits << ", "
      << m0 << ", "
      << m1 << ", "
      << m2 << ", "
      << m3 << ", "
      << m4 << ", "
      << avg(times) << endl;
  out.close();
}


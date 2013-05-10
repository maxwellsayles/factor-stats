#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <stdint.h>
#include <time.h>

#include <flint/flint.h>
#include <flint/qsieve.h>
#include <gmp.h>

#include "liboptarith/s128_c.h"
#include "factor-stats/factor-stats.h"

using namespace std;

string out_filename() {
  return "flint-timings.dat";
}

s128 factor(s128 x) {
  return qsieve_ll_factor(x.v1, x.v0);
}

int main(int argc, char** argv) {
  for (int i = 16; i <= 100; i += 2) {
    if (i != 24 && i != 28 && i != 42) {
      dobits(i);
    }
  }

  return 0;
}

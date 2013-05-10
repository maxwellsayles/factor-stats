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
#include "factor-stats/factor-stats.h"
#include "factor-stats/squfof-parigp/squfof64.h"
#include "factor-stats/squfof-parigp/squfof128.h"

using namespace std;

string out_filename() {
  return "squfof-timings.dat";
}

s128 factor(s128 N) {
  if (N.msb() >= 59) return squfof128(N);
  else squfof64(N.to_u64());
}

int main(int argc, char** argv) {
  for (int i = 16; i <= 76; i += 2) {
    dobits(i);
  }

  return 0;
}




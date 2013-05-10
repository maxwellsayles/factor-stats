#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

#include <pari/pari.h>

#include "factor-stats/factor-stats.h"
#include "liboptarith/s128_c.h"

// Convert s128 to GEN.
static inline GEN to_gen(const s128_t* x_) {
  long s = cmp_s128_s64(x_, 0);
  s128_t x = *x_;
  if (s < 0) {
    neg_s128_s128(&x, &x);
  }
  GEN r;
  if (x.v1 == 0) {
    r = stoi(x.v0);
  } else {
    r = mkintn(4,
               (uint64_t)(x.v1 >> 32) & 0xFFFFFFFF,
               (uint64_t)x.v1 & 0xFFFFFFFF,
               (x.v0 >> 32) & 0xFFFFFFFF,
               x.v0 & 0xFFFFFFFF);
  }
  if (s < 0) setsigne(r, -1);
  return r;
}

string out_filename() {
  return "pari-timings.dat";
}

s128 factor(s128 x) {
  pari_sp ltop = avma;

  GEN xg = to_gen(&x);
  Z_factor(xg);  // NOTE: This returns an array of factors.

  avma = ltop;
  return -1;  // NOTE: This tricks factor-stats into counting the result.
}

int main(int argc, char** argv) {
  pari_init(1<<30, 10000000);
  for (int i = 16; i <= 100; i += 2) {
    dobits(i);
  }
  return 0;
}

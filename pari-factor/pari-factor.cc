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

/// Convert a GEN into an s128_t.
static void to_s128(s128_t* x, GEN g) {
  long l = lgefint(g);
  if (l == 2) {
    setzero_s128(x);
  } else if (l == 3) {
    long* p = int_LSW(g);
    set_s128_u64(x, *p);
    if (signe(g) == -1) {
      neg_s128_s128(x, x);
    }
  } else if (l == 4) {
    long* p = int_LSW(g);
    x->v0 = *p;
    p = int_nextW(p);
    x->v1 = *p;
    if (signe(g) == -1) {
      neg_s128_s128(x, x);
    }
  } else {
    assert(false);
  }
}

string out_filename() {
  return "pari-timings.dat";
}

s128 factor(s128 x) {
  pari_sp ltop = avma;

  GEN xg = to_gen(&x);
  GEN yg = Z_factor(xg);
  s128_t y;
  to_s128(&y, yg);

  avma = ltop;
  return y;
}

int main(int argc, char** argv) {
  pari_init(1<<30, 10000000);
  for (int i = 16; i <= 100; i += 2) {
    dobits(i);
  }
  return 0;
}

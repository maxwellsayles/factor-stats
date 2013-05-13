#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <stdint.h>
#include <time.h>
#include <gmp.h>

#include "liboptarith/s128_c.h"
#include "spar/spar.h"
#include "../factor-stats.h"

using namespace std;

string out_filename() {
  return "spar-timings.dat";
}

Spar spar;
mpz_t N, d;

s128 factor(s128 x) {
  x.to_mpz(N);
  if (!spar.factor(d, N)) {
    return 1;
  }
  s128 r;
  r.from_mpz(d);
  return r;
}

int main(int argc, char** argv) {
  mpz_init(N);
  mpz_init(d);

  for (int i = 16; i <= 40; i += 2) {
    dobits(i);
  }

  mpz_clear(N);
  mpz_clear(d);

  return 0;
}

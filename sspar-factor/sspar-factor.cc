#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <stdint.h>
#include <time.h>
#include <gmp.h>

extern "C" {
#include "libsspar/sspar.h"
}

#include "liboptarith/s128_c.h"
#include "../factor-stats.h"

using namespace std;

string out_filename() {
  return "sspar-timings.dat";
}

sspar_t sspar;
mpz_t N, d;

s128 factor(s128 x) {
  x.to_mpz(N);
  sspar_factor(&sspar, d, N);
  s128 r;
  r.from_mpz(d);
  return r;
}

int main(int argc, char** argv) {
  mpz_init(N);
  mpz_init(d);
  sspar_init(&sspar);

  for (int i = 16; i <= 100; i += 2) {
    dobits(i);
  }

  mpz_clear(N);
  mpz_clear(d);
  sspar_clear(&sspar);

  return 0;
}

/**
 * An implementation of squfof taken from the Pari-GP source.
 */

#include <gmp.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <time.h>

#include "liboptarith/math64.h"
#include "liboptarith/gcd_binary_l2r.h"

static inline uint64_t current_nanos() {
  struct timespec res;
  clock_gettime(CLOCK_MONOTONIC, &res);
  return (res.tv_sec * 1000000000ULL) + res.tv_nsec;
}

int is_square(uint64_t A, uint64_t* sqrtA) {
  *sqrtA = sqrt_u64(A);
  return *sqrtA * *sqrtA == A;
}

/***********************************************************************/
/**                                                                   **/
/**              FACTORIZATION (Shanks' SQUFOF) --GN2000Sep30-Oct01   **/
/**  squfof() returns a nontrivial factor of n, assuming n is odd,    **/
/**  composite, not a pure square, and has no small prime divisor,    **/
/**  or NULL if it fails to find one.  It works on two discriminants  **/
/**  simultaneously  (n and 5n for n=1 (mod 4), 3n                    **/
/**  and 4n for n=3 (mod 4)).                                         **/
/**  Present implementation is limited to input <2^59, and works most **/
/**  of the time in signed arithmetic on integers <2^31 in absolute   **/
/**  size. (Cf. Algo 8.7.2 in ACiCNT)                                 **/
/**                                                                   **/
/***********************************************************************/

/* The following is invoked to walk back along the ambiguous cycle* until we
 * hit an ambiguous form and thus the desired factor, which it returns.  If it
 * fails for any reason, it returns 0.  It doesn't interfere with timing and
 * diagnostics, which it leaves to squfof().
 *
 * Before we invoke this, we've found a form (A, B, -C) with A = a^2, where a
 * isn't blacklisted and where gcd(a, B) = 1.  According to ACiCANT, we should
 * now proceed reducing the form (a, -B, -aC), but it is easy to show that the
 * first reduction step always sends this to (-aC, B, a), and the next one,
 * with q computed as usual from B and a (occupying the c position), gives a
 * reduced form, whose third member is easiest to recover by going back to D.
 * From this point onwards, we're once again working with single-word numbers.
 * No need to track signs, just work with the abs values of the coefficients. */
int64_t squfof_ambig(int64_t a, int64_t B, int64_t dd, int64_t D) {
  int64_t b, c, q, qc, qcb, a0, b0, b1, c0;

  q = (dd + (B >> 1)) / a;
  b = ((q * a) << 1) - B;

  c = ((D - b * b) >> 2) / a;

  a0 = a;
  b0 = b; // end of loop detection and safeguard
  b1 = b; 

  // reduced cycles are finite
  // this is the reduction step
  for (;;) {
    c0 = c;
    if (c0 > dd) {
      qcb = c0 - b;
      b = c0 + qcb;
      c = a - qcb;
    } else {
      q = (dd + (b >> 1)) / c0;
      if (q == 1) {
	qcb = c0 - b;
	b = c0 + qcb;
	c = a - qcb;
      } else {
	qc = q*c0;
	qcb = qc - b;
	b = qc + qcb;
	c = a - q*qcb;
      }
    }
    a = c0;

    if (b == b1) break;

    /* safeguard against infinite loop: recognize when we've walked the entire
     * cycle in vain. (I don't think this can actually happen -- exercise.) */
    if (b == b0 && a == a0) {
      return 0;
    }

    b1 = b;
  }
  q = a & 1 ? a : a >> 1;

  return q;
}

#define SQUFOF_BLACKLIST_SZ 256

/* assume 2,3,5 do not divide n */
uint64_t squfof(uint64_t n) {
  int64_t d1, d2;
  uint64_t nm4; // n (mod 4)
  int cnt = 0;
  int64_t a1, b1, c1, dd1, L1;
  int64_t a2, b2, c2, dd2, L2;
  int64_t a, q, c, qc, qcb;
  int64_t D1, D2;
  int64_t blacklist1[SQUFOF_BLACKLIST_SZ];
  int64_t blacklist2[SQUFOF_BLACKLIST_SZ];
  int blp1 = 0; // black list pointer
  int blp2 = 0;
  int act1 = 1; // is this multiple of N active
  int act2 = 1;
  int j;

  // now we have 5 < n < 2^59
  nm4 = n & 3;
  if (nm4 == 1) { // n = 1 (mod4):  run one iteration on D1 = n, another on D2 = 5n
    D1 = n;
    D2 = 5 * n;
    d2 = (int64_t)sqrt_u64(D2);
    dd2 = (d2 >> 1) + (d2 & 1);
    b2 = (d2 - 1) | 1; // b1, b2 will always stay odd
  } else { // n = 3 (mod4):  run one iteration on D1 = 3n, another on D2 = 4n
    D1 = 3 * n;
    D2 = 4 * n;
    dd2 = (int64_t)sqrt_u64(D2);
    d2 = dd2 << 1;
    b2 = d2;// & (~1UL); // largest even below d2, will stay even
  }
  d1 = sqrt_u64(D1);
  b1 = (d1 - 1) | 1; // largest odd number not exceeding d1

  c1 = (D1 - b1 * b1) >> 2;
  c2 = (D2 - b2 * b2) >> 2;
  L1 = sqrt_u64(d1);
  L2 = sqrt_u64(d2);

  /* dd1 used to compute floor((d1+b1)/2) as dd1+floor(b1/2), without
   * overflowing the 31bit signed integer size limit. Same for dd2. */
  dd1 = (d1 >> 1) + (d1 & 1);
  a1 = 1;
  a2 = 1;

  /* The two (identity) forms (a1,b1,-c1) and (a2,b2,-c2) are now set up.
   *
   * a1 and c1 represent the absolute values of the a,c coefficients; we keep
   * track of the sign separately, via the iteration counter cnt: when cnt is
   * even, c is understood to be negative, else c is positive and a < 0.
   *
   * L1, L2 are the limits for blacklisting small leading coefficients
   * on the principal cycle, to guarantee that when we find a square form,
   * its square root will belong to an ambiguous cycle  (i.e. won't be an
   * earlier form on the principal cycle).
   *
   * When n = 3(mod 4), D2 = 12(mod 16), and b^2 is always 0 or 4 mod 16.
   * It follows that 4*a*c must be 4 or 8 mod 16, respectively, so at most
   * one of a,c can be divisible by 2 at most to the first power.  This fact
   * is used a couple of times below.
   *
   * The flags act1, act2 remain true while the respective cycle is still
   * active;  we drop them to false when we return to the identity form with-
   * out having found a square form  (or when the blacklist overflows, which
   * shouldn't happen). */

  /* MAIN LOOP: walk around the principal cycle looking for a square form.
   * Blacklist small leading coefficients.
   *
   * The reduction operator can be computed entirely in 32-bit arithmetic:
   * Let q = floor(floor((d1+b1)/2)/c1)  (when c1>dd1, q=1, which happens
   * often enough to special-case it).  Then the new b1 = (q*c1-b1) + q*c1,
   * which does not overflow, and the new c1 = a1 - q*(q*c1-b1), which is
   * bounded by d1 in abs size since both the old and the new a1 are positive
   * and bounded by d1. */
  while (act1 || act2) {
    if (act1) { // send first form through reduction operator if active
      c = c1;
      q = (c > dd1) ? 1 : (dd1 + (b1 >> 1)) / c;
      if (q == 1) {
	qcb = c - b1;
	b1 = c + qcb;
	c1 = a1 - qcb;
      } else {
	qc = q*c;
	qcb = qc - b1;
	b1 = qc + qcb;
	c1 = a1 - q*qcb;
      }
      a1 = c;

      if (a1 <= L1) { // blacklist this
	if (blp1 >= SQUFOF_BLACKLIST_SZ) {
	  // overflows: shouldn't happen
	  printf("Black list for second discriminant overflowed.  Increase blacklist.\n");
	  exit(-1);
	} else {
	  blacklist1[blp1++] = a1;
	}
      }
    }
    if (act2) { // send second form through reduction operator if active
      c = c2;
      q = (c > dd2) ? 1 : (dd2 + (b2 >> 1)) / c;
      if (q == 1) {
	qcb = c - b2;
	b2 = c + qcb;
	c2 = a2 - qcb;
      } else {
	qc = q*c;
	qcb = qc - b2;
	b2 = qc + qcb;
	c2 = a2 - q*qcb;
      }
      a2 = c;

      if (a2 <= L2) { // blacklist this
	if (blp2 >= SQUFOF_BLACKLIST_SZ) {
	  // overflows: shouldn't happen
	  printf("Black list for second discriminant overflowed.  Increase blacklist.\n");
	  exit(-1);
	} else {
	  blacklist2[blp2++] = a2;
	}
      }
    }

    // bump counter, loop if this is an odd iteration (i.e. if the real
    // leading coefficients are negative)
    cnt ++;
    if (cnt & 1)
      continue;

    // second half of main loop entered only when the leading coefficients
    // are positive (i.e., during even-numbered iterations)

    // examine first form if active
    if (act1 && a1 == 1) {
      // back to identity
      // drop this discriminant
      act1 = 0;
    }
    if (act1) {
      if (is_square((uint64_t) a1, (uint64_t*)&a)) { // square form
	// check if this is blacklisted
	if (a <= L1) {
	  for (j = 0; j < blp1; j++)
	    if (a == blacklist1[j]) {
	      a = 0;
	      break;
	    }
	}
	// not blacklisted
	if (a > 0) {
	  q = gcd_binary_l2r_u64(a, b1); // imprimitive form?
	  if (q > 1) { /* q^2 divides D1 hence n [ assuming n % 3 != 0 ] */
	    return q*q;
	  }
	  /* chase the inverse root form back along the ambiguous cycle */
	  q = squfof_ambig(a, b1, dd1, D1);
	  if (nm4 == 3 && q % 3 == 0) q /= 3;
	  if (q > 1) {
	    return q;
	  }
	}
      }
    }

    // examine second form if active
    if (act2 && a2 == 1) {
      // back to identity
      // drop this discriminant
      act2 = 0;
    }
    if (act2) {
      if (is_square((uint64_t)a2, (uint64_t*)&a)) { // square form
	// check if this is blacklisted
	if (a <= L2) {
	  for (j = 0; j < blp2; j++)
	    if (a == blacklist2[j]) {
	      a = 0;
	      break;
	    }
	}
	// not blacklisted
	if (a > 0) {
	  q = gcd_binary_l2r_u64(a, b2); // imprimitive form?
	  /* NB if b2 is even, a is odd, so the gcd is always odd */
	  if (q > 1) { /* q^2 divides D2 hence n [ assuming n % 5 != 0 ] */
	    return q*q;
	  }
	  /* chase the inverse root form along the ambiguous cycle */
	  q = squfof_ambig(a, b2, dd2, D2);
	  if (nm4 == 1 && q % 5 == 0) q /= 5;
	  if (q > 1) {
	    return q;
	  }
	}
      }
    }
  }

  // both discriminants turned out to be useless.
  return 0;
}

int main(int argc, char** argv) {
  int i = 0;
  int count = 0;
  int t = 0;
  int failed = 0;
  uint64_t tmp = 0;    
  uint64_t start_time = 0;
  uint64_t finish_time = 0;

  if (argc != 2) {
    printf("Usage: %s <composites.txt>\n", argv[0]);
    printf("\n");
    exit(0);
  }
  
  // read numbers in two passes
  // first pass determines how many numbers, second pass reads them
  FILE* f = fopen(argv[1], "r");
  if (!f) {
    printf("Couldn't find file '%s'\n", argv[1]);
    exit(-1);
  }
  count = 0;
  while (!feof(f)) {
    t = fscanf(f, "%"PRIu64, &tmp);
    count ++;
  }
  count --;
  fseek(f, 0, SEEK_SET);
  printf("Factoring %d composites.\n", count);

  // allocate memory and read in the numbers
  uint64_t* composites = (uint64_t*)malloc(sizeof(uint64_t)*count);

  for (i = 0; i < count; i++) {
    t = fscanf(f, "%"PRIu64, &composites[i]);
  }

  fclose(f);

  start_time = current_nanos();
  for (i = 0; i < count; i++) {
    uint64_t N = composites[i];
    uint64_t d = squfof(N);
        
    if (d == 0 || d == N || d == 1) {
      //printf("Failed to factor %"PRIu64"\n", N);
      failed ++;
    } else {
      //printf("%"PRIu64" = %"PRIu64" * %"PRIu64"\n", N, d, N / d);
      if (N % d) {
	printf("Insane N=%"PRIu64", d=%"PRIu64"\n", N, d);
	exit(-1);
      }
    }
  }
  finish_time = current_nanos();

  free(composites);

  printf("Failed to factor %d integers\n", failed);
  
  // Extrapolate time based on successes.
  uint64_t guessed_time = (finish_time - start_time) * count / (count - failed);

  printf("Took %0.2f milliseconds for %d composites.\n", (finish_time - start_time)/1000000.0, count);
  printf("Extrapolated time = %0.5f seconds.\n", (double)guessed_time / 1000000000.0);

  return 0;
}


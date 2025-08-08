// -*- coding: utf-8 -*-

#pragma once

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#define FE_PAIR_IMPLEMENTATION
#include "../f64_pair.h"

#define REPORT_TABLE_IMPLEMENTATION
#include "report_table.h"

#define  PRNG_IMPLEMENTATION
#include "prng_small_global.h"

#include <mpfr.h>

#define LENGTHOF(X) (sizeof(X)/sizeof(X[0]))

#define DEF_FE(N,M) .fe=N,.name=#N,.mp=&M
#define DEF_FR(N,M) .fr=N,.name=#N,.mp=&M
#define OP_U(U,X,R) .max_ulp=U,.max_x=X,.max_r=R
#define OP_B(U,A,B,R) .max_ulp=U,.max_a=A,.max_b=B,.max_r=R

#define UP(H,L) {.hi=H,.lo=L}

static const uint64_t f64_one_bits_k = UINT64_C(0x3ff0000000000000);


// random on [1,2)
// ensures hi and lo bits are fully populated
static inline fe_pair_t prng_fe_12(void)
{
  uint64_t u = prng_u64() >> 12;         // uniform 64-bit draw (reduce to 53)
  uint64_t b = u ^ f64_one_bits_k;       // [1,2) {bit pattern}
  double   h = type_pun(b,double);       // [1,2) {binary64}
  double   l;

  u = prng_u64() >> 12;                  // uniform 64-bit draw (reduce to 53)
  b = u ^ f64_one_bits_k;                // [1,2) {bit pattern}
  l = type_pun(b,double);                // [1,2) {binary64}

  l *= 0x1.0p-54;
  
  return fe_pair(h,l);
}

// on [0,1)
// hi bits are fully populated, lo is "almost" (is odd 25% of the time)
static inline fe_pair_t prng_fe(void)
{
  double   hi = prng_odd_f64();
  uint64_t u  = prng_u64() >> 12;
  double   lo = hi * 0x1.0p-54;

  lo = type_pun(type_pun(lo,uint64_t) ^ u, double);

  return fe_fast_sum(hi,lo);
}


static inline void mp_set(mpfr_t r, fe_pair_t x)
{
  mpfr_set_d(r,  x.hi, MPFR_RNDN);
  mpfr_add_d(r,r,x.lo, MPFR_RNDN);
};

extern mpfr_t mp_t;

static fe_pair_t mp2fe(mpfr_t r)
{
  double hi = mpfr_get_d(r,MPFR_RNDN);

  mpfr_set_d(mp_t,hi,   MPFR_RNDN);
  mpfr_sub(mp_t,r,mp_t, MPFR_RNDN);
  
  double lo = mpfr_get_d(mp_t,MPFR_RNDN);
  
  return fe_pair(hi,lo);
};

// WARNING: modifies global 'mp_e' and passed in 'mp_t'
// assumes finite and normal
extern mpfr_t mp_e;

double ulp_dist(mpfr_t mp_r, fe_pair_t r)
{
  mpfr_exp_t rexp = mpfr_get_exp(mp_r);

  mp_set(mp_e,r);

  mpfr_set(mp_t,mp_r,MPFR_RNDN);
  mpfr_sub(mp_t,mp_t,mp_e,      MPFR_RNDN);
  mpfr_abs(mp_t,mp_t,           MPFR_RNDN);
  mpfr_mul_2si(mp_t,mp_t,53+54-rexp, MPFR_RNDN);
  
  return mpfr_get_d(mp_t, MPFR_RNDU);
}

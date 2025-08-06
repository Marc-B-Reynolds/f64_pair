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

#define LENGTHOF(X) (sizeof(X)/sizeof(X[0]))


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

// -*- coding: utf-8 -*-

#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>

#include "../f64_pair.h"

#define REPORT_TABLE_IMPLEMENTATION
#include "report_table.h"

#define  PRNG_IMPLEMENTATION
#include "prng_small_global.h"

#include <mpfr.h>

#define TRIALS 0x7fffff

#ifndef TRIALS
#define TRIALS 0x1fff
#endif


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

// globals 
mpfr_t mp_r0;
mpfr_t mp_r1;
mpfr_t mp_a;
mpfr_t mp_b;
mpfr_t mp_e;
mpfr_t mp_t;


static inline int mpfr_inv(mpfr_t r ,const mpfr_t x, mpfr_rnd_t rnd)
{
  return mpfr_ui_div(r,1,x,rnd);
}

static inline void mp_set(mpfr_t r, fe_pair_t x)
{
  mpfr_set_d(r,  x.hi, MPFR_RNDN);
  mpfr_add_d(r,r,x.lo, MPFR_RNDN);
};

// test this
static fe_pair_t mp2fe(mpfr_t r)
{
  double hi = mpfr_get_d(r,MPFR_RNDN);

  mpfr_set_d(mp_t,hi,   MPFR_RNDN);
  mpfr_sub(mp_t,r,mp_t, MPFR_RNDN);
  
  double lo = mpfr_get_d(mp_t,MPFR_RNDN);
  
  return fe_pair(hi,lo);
};


// validate round-trip to MPFR value is bit exact
void test_mp2fe(void)
{
  for(int i=0; i<0xfffff; i++) {
    fe_pair_t a = prng_fe();
    mp_set(mp_r0,a);
    fe_pair_t r = mp2fe(mp_r0);

    if (a.hi == r.hi && a.lo == r.lo) continue;

    printf("{%a,%a} {%a,%a} {%a,%a}\n",a.hi,a.lo,r.hi,r.lo,a.hi-r.hi,a.lo-r.lo);
  }
}


// WARNING: modifies global 'mp_e' and passed in 'mp_t'
// assumes finite and normal
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

void test_test(void)
{
  double max_ulp = 0.0;

  for(int i=0; i<0xfffff; i++) {
    fe_pair_t a = prng_fe();

    mp_set(mp_r0,a);

    a.lo = type_pun(type_pun(a.lo,uint64_t)^1,double);

    double ulp = ulp_dist(mp_r0,a);
    
    if (max_ulp >= ulp) continue;
    
    max_ulp = ulp;
    printf("%f ", ulp);
  }

  printf("\n%f\n\n", max_ulp);
}


//**********************************************************
// unary ops (double input)

typedef struct {
  union {
    fe_pair_t (*fe)(double);
    fr_pair_t (*fr)(double);
  };

  int (*mp)(mpfr_t,const mpfr_t,mpfr_rnd_t);
  char* name;

  double    max_ulp;
  double    max_x;
  fe_pair_t max_r;
  fe_pair_t max_e;
  
} uop_d_table_t;


#define DEF_FE(N,M) .fe=N,.name=#N,.mp=&M
#define DEF_FR(N,M) .fr=N,.name=#N,.mp=&M

uop_d_table_t op_d[] =
{
  { DEF_FE(fe_rsqrt_d,  mpfr_rec_sqrt)  },
  { DEF_FE(fe_rsqrt_dh, mpfr_rec_sqrt)  },
  { DEF_FE(fe_inv_d,    mpfr_inv)  },
  { DEF_FE(fe_inv_dn,   mpfr_inv)  },
  { DEF_FE(fe_inv_dh,   mpfr_inv)  },
  { DEF_FR(fr_sqrt_d,   mpfr_sqrt) },
};


// test unary f(double) → pair
void op_d_tests(void)
{
  printf(SGR_BOLD SGR_RGB(200,200,255) "\nf(double) → pair\n" SGR_RESET);
  
  report_table_t table = {
    .col = {
      { REPORT_TABLE_STR("f(x)",14), .just=report_table_justify_left },
      { REPORT_TABLE_F("ulp",1,4) },
      { REPORT_TABLE_A64("x") },
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
      { REPORT_TABLE_A64("e.hi") },      
      { REPORT_TABLE_A64("e.lo") },      
    }
  };
  
  report_table_header(stdout, &table);
  
  for(size_t i=0; i<LENGTHOF(op_d); i++) {
    double    max_ulp = 0.0;
    double    max_a;
    fe_pair_t max_r;
    fe_pair_t max_e;
    
    fe_pair_t (*f)(double) = op_d[i].fe;
    int       (*mp)(mpfr_t,const mpfr_t,mpfr_rnd_t) = op_d[i].mp;
    
    for(int i=0; i<TRIALS; i++) {
      double    a = prng_f64();
      fe_pair_t r = f(a);
      
      // perform with MPFR
      mpfr_set_d(mp_a,a,MPFR_RNDN);
      mp(mp_r0,mp_a,MPFR_RNDN);
      fe_pair_t e = mp2fe(mp_r0);
      
      double ulp = ulp_dist(mp_r0,r);
      
      if (max_ulp >= ulp) continue;
      
      max_ulp = ulp;
      max_a   = a;
      max_r   = r;
      max_e   = e;
    }
    
    if (max_ulp > op_d[i].max_ulp) {
      op_d[i].max_ulp = max_ulp;
      op_d[i].max_x   = max_a;
      op_d[i].max_r   = max_r;
      op_d[i].max_e   = max_e;
    }
    
    report_table_row(stdout,&table,
                     op_d[i].name,
                     op_d[i].max_ulp,
                     op_d[i].max_x,
                     op_d[i].max_r.hi,
                     op_d[i].max_r.lo,
                     op_d[i].max_e.hi,
                     op_d[i].max_e.lo
                     );
  }
  
  report_table_end(stdout, &table);
}


//**********************************************************
// unary ops (pair input)

typedef enum {
  UOP_TYPE_FE,
  UOP_TYPE_FR,
} uop_type_t;


typedef struct {
  uop_type_t type;
  
  union {
    fe_pair_t (*fe)(fe_pair_t);
    fr_pair_t (*fr)(fr_pair_t);
  };

  int (*mp)(mpfr_t,const mpfr_t,mpfr_rnd_t);
  char* name;

  //
  double    max_ulp;
  fe_pair_t max_x;
  fe_pair_t max_r;
  fe_pair_t max_e;
  
} uop_table_t;


static inline fe_pair_t foo_sqrt(fe_pair_t x)
{
  (void)x;
  return fe_fast_sum(0,0);
}


uop_table_t uops[] =
{
  { DEF_FR(fr_inv,   mpfr_inv)  },
  { DEF_FR(fr_inv_a, mpfr_inv)  },
  { DEF_FE(fe_inv,   mpfr_inv)  },
  { DEF_FE(fe_sqrt,  mpfr_sqrt) },
  { DEF_FE(fe_sq,    mpfr_sqr)  },
};

// test unary
void op_p_tests(void)
{
  printf(SGR_BOLD SGR_RGB(200,200,255) "\nf(pair) → pair\n" SGR_RESET);

  report_table_t table = {
    .col = {
      { REPORT_TABLE_STR("f(x)",10), .just=report_table_justify_left },
      { REPORT_TABLE_F("ulp",1,4) },
      { REPORT_TABLE_A64("x.hi") },      
      { REPORT_TABLE_A64("x.lo") },
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
      { REPORT_TABLE_A64("e.hi") },      
      { REPORT_TABLE_A64("e.lo") },      
    }
  };
  
  report_table_header(stdout, &table);
  
  for(size_t i=0; i<LENGTHOF(uops); i++) {
    double    max_ulp = -2.0;
    fe_pair_t max_a;
    fe_pair_t max_r;
    fe_pair_t max_e;
    
    fe_pair_t (*f)(fe_pair_t) = uops[i].fe;
    int       (*mp)(mpfr_t,const mpfr_t,mpfr_rnd_t) = uops[i].mp;
    
    for(int i=0; i<TRIALS; i++) {
      fe_pair_t a = prng_fe_12();
      fe_pair_t r = f(a);
      
      // perform with MPFR
      mp_set(mp_a, a);
      mp(mp_r0,mp_a,MPFR_RNDN);
      fe_pair_t e = mp2fe(mp_r0);
      
      double ulp = ulp_dist(mp_r0,r);

      if (max_ulp >= ulp) continue;
      
      max_ulp = ulp;
      max_a   = a;
      max_r   = r;
      max_e   = e;
    }
    
    if (max_ulp > op_d[i].max_ulp) {
      uops[i].max_ulp = max_ulp;
      uops[i].max_x   = max_a;
      uops[i].max_r   = max_r;
      uops[i].max_e   = max_e;
    }
    
    report_table_row(stdout,&table,
                     uops[i].name,
                     uops[i].max_ulp,
                     uops[i].max_x,
                     uops[i].max_r.hi,
                     uops[i].max_r.lo,
                     uops[i].max_e.hi,
                     uops[i].max_e.lo
                     );
  }
  
  report_table_end(stdout, &table);
}




//**********************************************************
// binary ops
//**********************************************************

//**********************************************************
// f(double,double) 

typedef struct {
  union {
    fe_pair_t (*fe)(double,double);
    fr_pair_t (*fr)(double,double);
  };

  int (*mp)(mpfr_t,const mpfr_t,const mpfr_t,mpfr_rnd_t);
  char* name;

  double    max_ulp;
  double    max_a;
  double    max_b;
  fe_pair_t max_r;
  fe_pair_t max_e;
  
} uop_dd_table_t;

uop_dd_table_t op_dd[] =
{
  { DEF_FE(fe_div_dd,  mpfr_div)  },
  { DEF_FR(fr_div_dd,  mpfr_div)  },
};

void op_dd_tests(void)
{
  printf(SGR_BOLD SGR_RGB(200,200,255) "\nf(double,double) → pair\n" SGR_RESET);
  
  report_table_t table = {
    .col = {
      { REPORT_TABLE_STR("f(a,b)",14), .just=report_table_justify_left },
      { REPORT_TABLE_F("ulp",1,4) },
      { REPORT_TABLE_A64("a") },
      { REPORT_TABLE_A64("b") },
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
      { REPORT_TABLE_A64("e.hi") },      
      { REPORT_TABLE_A64("e.lo") },      
    }
  };
  
  report_table_header(stdout, &table);
  
  for(size_t i=0; i<LENGTHOF(op_dd); i++) {
    double    max_ulp = 0.0;
    double    max_a;
    double    max_b;
    fe_pair_t max_r;
    fe_pair_t max_e;
    
    fe_pair_t (*f)(double,double) = op_dd[i].fe;
    int       (*mp)(mpfr_t,const mpfr_t,const mpfr_t,mpfr_rnd_t) = op_dd[i].mp;
    
    for(int i=0; i<TRIALS; i++) {
      double    a = prng_f64();
      double    b = prng_f64();
      fe_pair_t r = f(a,b);
      
      // perform with MPFR
      mpfr_set_d(mp_a,a,MPFR_RNDN);
      mpfr_set_d(mp_b,b,MPFR_RNDN);
      mp(mp_r0,mp_a,mp_b,MPFR_RNDN);
      fe_pair_t e = mp2fe(mp_r0);
      
      double ulp = ulp_dist(mp_r0,r);
      
      if (max_ulp >= ulp) continue;
      
      max_ulp = ulp;
      max_a   = a;
      max_b   = b;
      max_r   = r;
      max_e   = e;
    }
    
    if (max_ulp > op_d[i].max_ulp) {
      op_dd[i].max_ulp = max_ulp;
      op_dd[i].max_a   = max_a;
      op_dd[i].max_b   = max_b;
      op_dd[i].max_r   = max_r;
      op_dd[i].max_e   = max_e;
    }
    
    report_table_row(stdout,&table,
                     op_dd[i].name,
                     op_dd[i].max_ulp,
                     op_dd[i].max_a,
                     op_dd[i].max_b,
                     op_dd[i].max_r.hi,
                     op_dd[i].max_r.lo,
                     op_dd[i].max_e.hi,
                     op_dd[i].max_e.lo
                     );
  }
  
  report_table_end(stdout, &table);
}


//**********************************************************
// f(double,pair) 

typedef struct {
  union {
    fe_pair_t (*fe)(double,fe_pair_t);
    fr_pair_t (*fr)(double,fe_pair_t);
  };

  int (*mp)(mpfr_t,const mpfr_t,const mpfr_t,mpfr_rnd_t);
  char* name;

  double    max_ulp;
  double    max_a;
  fe_pair_t max_b;
  fe_pair_t max_r;
  fe_pair_t max_e;
  
} uop_dp_table_t;

uop_dp_table_t op_dp[] =
{
  { DEF_FE(fe_d_add,  mpfr_add)  },
  { DEF_FE(fe_d_div,  mpfr_div)  },
};

void op_dp_tests(void)
{
  printf(SGR_BOLD SGR_RGB(200,200,255) "\nf(double,pair) → pair\n" SGR_RESET);
  
  report_table_t table = {
    .col = {
      { REPORT_TABLE_STR("f(a,b)",14), .just=report_table_justify_left },
      { REPORT_TABLE_F("ulp",1,4) },
      { REPORT_TABLE_A64("a") },
      { REPORT_TABLE_A64("b.hi") },      
      { REPORT_TABLE_A64("b.lo") },      
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
      { REPORT_TABLE_A64("e.hi") },      
      { REPORT_TABLE_A64("e.lo") },      
    }
  };
  
  report_table_header(stdout, &table);
  
  for(size_t i=0; i<LENGTHOF(op_dp); i++) {
    double    max_ulp = 0.0;
    double    max_a;
    fe_pair_t max_b;
    fe_pair_t max_r;
    fe_pair_t max_e;
    
    fe_pair_t (*f)(double,fe_pair_t) = op_dp[i].fe;
    int       (*mp)(mpfr_t,const mpfr_t,const mpfr_t,mpfr_rnd_t) = op_dp[i].mp;
    
    for(int i=0; i<TRIALS; i++) {
      double    a = prng_f64();
      fe_pair_t b = prng_fe();
      fe_pair_t r = f(a,b);
      
      // perform with MPFR
      mpfr_set_d(mp_a,a,MPFR_RNDN);
      mp_set(mp_b,b);
      mp(mp_r0,mp_a,mp_b,MPFR_RNDN);
      fe_pair_t e = mp2fe(mp_r0);
      
      double ulp = ulp_dist(mp_r0,r);
      
      if (max_ulp >= ulp) continue;
      
      max_ulp = ulp;
      max_a   = a;
      max_b   = b;
      max_r   = r;
      max_e   = e;
    }
    
    if (max_ulp > op_dp[i].max_ulp) {
      op_dp[i].max_ulp = max_ulp;
      op_dp[i].max_a   = max_a;
      op_dp[i].max_b   = max_b;
      op_dp[i].max_r   = max_r;
      op_dp[i].max_e   = max_e;
    }
    
    report_table_row(stdout,&table,
                     op_dp[i].name,
                     op_dp[i].max_ulp,
                     op_dp[i].max_a,
                     op_dp[i].max_b.hi,
                     op_dp[i].max_b.lo,
                     op_dp[i].max_r.hi,
                     op_dp[i].max_r.lo,
                     op_dp[i].max_e.hi,
                     op_dp[i].max_e.lo
                     );
  }
  
  report_table_end(stdout, &table);
}


//**********************************************************
// f(pair,double) 

typedef struct {
  union {
    fe_pair_t (*fe)(fe_pair_t,double);
    fr_pair_t (*fr)(fr_pair_t,double);
  };

  int (*mp)(mpfr_t,const mpfr_t,const mpfr_t,mpfr_rnd_t);
  char* name;

  double    max_ulp;
  fe_pair_t max_a;
  double    max_b;
  fe_pair_t max_r;
  fe_pair_t max_e;
  
} uop_pd_table_t;

uop_pd_table_t op_pd[] =
{
  { DEF_FE(fe_mul_d,  mpfr_mul)  },
  { DEF_FR(fr_mul_d,  mpfr_mul)  },
  { DEF_FE(fe_div_d,  mpfr_div)  },
};

void op_pd_tests(void)
{
  printf(SGR_BOLD SGR_RGB(200,200,255) "\nf(pair,double) → pair\n" SGR_RESET);
  
  report_table_t table = {
    .col = {
      { REPORT_TABLE_STR("f(a,b)",14), .just=report_table_justify_left },
      { REPORT_TABLE_F("ulp",1,4) },
      { REPORT_TABLE_A64("a.hi") },
      { REPORT_TABLE_A64("a.lo") },      
      { REPORT_TABLE_A64("b") },
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
      { REPORT_TABLE_A64("e.hi") },      
      { REPORT_TABLE_A64("e.lo") },      
    }
  };
  
  report_table_header(stdout, &table);
  
  for(size_t i=0; i<LENGTHOF(op_pd); i++) {
    double    max_ulp = 0.0;
    fe_pair_t max_a;
    double    max_b;
    fe_pair_t max_r;
    fe_pair_t max_e;
    
    fe_pair_t (*f)(fe_pair_t,double) = op_pd[i].fe;
    int       (*mp)(mpfr_t,const mpfr_t,const mpfr_t,mpfr_rnd_t) = op_pd[i].mp;
    
    for(int i=0; i<TRIALS; i++) {
      fe_pair_t a = prng_fe();
      double    b = prng_f64();
      fe_pair_t r = f(a,b);
      
      // perform with MPFR
      mp_set(mp_a,a);
      mpfr_set_d(mp_b,b,MPFR_RNDN);
      mp(mp_r0,mp_a,mp_b,MPFR_RNDN);
      fe_pair_t e = mp2fe(mp_r0);
      
      double ulp = ulp_dist(mp_r0,r);
      
      if (max_ulp >= ulp) continue;
      
      max_ulp = ulp;
      max_a   = a;
      max_b   = b;
      max_r   = r;
      max_e   = e;
    }
    
    if (max_ulp > op_pd[i].max_ulp) {
      op_pd[i].max_ulp = max_ulp;
      op_pd[i].max_a   = max_a;
      op_pd[i].max_b   = max_b;
      op_pd[i].max_r   = max_r;
      op_pd[i].max_e   = max_e;
    }
    
    report_table_row(stdout,&table,
                     op_pd[i].name,
                     op_pd[i].max_ulp,
                     op_pd[i].max_a.hi,
                     op_pd[i].max_a.lo,
                     op_pd[i].max_b,
                     op_pd[i].max_r.hi,
                     op_pd[i].max_r.lo,
                     op_pd[i].max_e.hi,
                     op_pd[i].max_e.lo
                     );
  }
  
  report_table_end(stdout, &table);
}


//**********************************************************
// f(pair,pair) 

typedef struct {
  union {
    fe_pair_t (*fe)(fe_pair_t,fe_pair_t);
    fr_pair_t (*fr)(fr_pair_t,fr_pair_t);
  };

  int (*mp)(mpfr_t,const mpfr_t,const mpfr_t,mpfr_rnd_t);
  char* name;

  double    max_ulp;
  fe_pair_t max_a;
  fe_pair_t max_b;
  fe_pair_t max_r;
  fe_pair_t max_e;
  
} uop_pp_table_t;

uop_pp_table_t op_pp[] =
{
  { DEF_FR(fr_div,   mpfr_div) },
  { DEF_FR(fr_div_a, mpfr_div) },

  { DEF_FE(fe_sub,   mpfr_sub) },
  { DEF_FR(fr_sub,   mpfr_sub) },
  { DEF_FE(fe_osub,  mpfr_sub) },
  { DEF_FR(fr_osub,  mpfr_sub) },
  { DEF_FE(fe_add,   mpfr_add) },
  { DEF_FR(fr_add,   mpfr_add) },
  { DEF_FE(fe_add_s, mpfr_add) },
  { DEF_FE(fe_mul,   mpfr_mul) },
  { DEF_FR(fr_mul,   mpfr_mul) },
  { DEF_FE(fe_div,   mpfr_div) },
  { DEF_FE(fe_sub_s, mpfr_sub) },
};

void op_pp_tests(void)
{
  printf(SGR_BOLD SGR_RGB(200,200,255) "\nf(pair,pair) → pair : a∈[1,2) b∈[0,1) \n" SGR_RESET);
  
  report_table_t table = {
    .col = {
      { REPORT_TABLE_STR("f(a,b)",14), .just=report_table_justify_left },
      { REPORT_TABLE_F("ulp",8,4) },
      { REPORT_TABLE_A64("a.hi") },
      { REPORT_TABLE_A64("a.lo") },      
      { REPORT_TABLE_A64("b.hi") },
      { REPORT_TABLE_A64("b.lo") },      
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
      { REPORT_TABLE_A64("e.hi") },      
      { REPORT_TABLE_A64("e.lo") },      
    }
  };
  
  report_table_header(stdout, &table);
  
  for(size_t i=0; i<LENGTHOF(op_pp); i++) {
    double    max_ulp = 0.0;
    fe_pair_t max_a;
    fe_pair_t max_b;
    fe_pair_t max_r;
    fe_pair_t max_e;
    
    fe_pair_t (*f)(fe_pair_t,fe_pair_t) = op_pp[i].fe;
    int       (*mp)(mpfr_t,const mpfr_t,const mpfr_t,mpfr_rnd_t) = op_pp[i].mp;
    
    for(int i=0; i<TRIALS; i++) {
      fe_pair_t a = prng_fe_12();
      fe_pair_t b = prng_fe();
      fe_pair_t r = f(a,b);
      
      // perform with MPFR
      mp_set(mp_a,a);
      mp_set(mp_b,b);
      mp(mp_r0,mp_a,mp_b,MPFR_RNDN);
      fe_pair_t e = mp2fe(mp_r0);
      
      double ulp = ulp_dist(mp_r0,r);
      
      if (max_ulp >= ulp) continue;
      
      max_ulp = ulp;
      max_a   = a;
      max_b   = b;
      max_r   = r;
      max_e   = e;
    }
    
    if (max_ulp > op_pp[i].max_ulp) {
      op_pp[i].max_ulp = max_ulp;
      op_pp[i].max_a   = max_a;
      op_pp[i].max_b   = max_b;
      op_pp[i].max_r   = max_r;
      op_pp[i].max_e   = max_e;
    }
    
    report_table_row(stdout,&table,
                     op_pp[i].name,
                     op_pp[i].max_ulp,
                     op_pp[i].max_a.hi,
                     op_pp[i].max_a.lo,
                     op_pp[i].max_b,
                     op_pp[i].max_r.hi,
                     op_pp[i].max_r.lo,
                     op_pp[i].max_e.hi,
                     op_pp[i].max_e.lo
                     );
  }
  
  report_table_end(stdout, &table);
}


//**********************************************************

int main(void)
{
  mpfr_set_emin(-1074);
  mpfr_set_emax( 1024);

  mpfr_init2(mp_r0, 128);
  mpfr_init2(mp_r1, 128);
  mpfr_init2(mp_a,  128);
  mpfr_init2(mp_b,  128);

  mpfr_init2(mp_e,  128);
  mpfr_init2(mp_t,  128);

  test_mp2fe();

#if 0
//op_d_tests();
//op_p_tests();
//op_dd_tests();
//op_dp_tests();
//op_pd_tests();
  op_pp_tests();
#else  
  op_d_tests();
  op_p_tests();

  op_dd_tests();
  op_dp_tests();
  op_pd_tests();
  op_pp_tests();
#endif  

  return 0;
}



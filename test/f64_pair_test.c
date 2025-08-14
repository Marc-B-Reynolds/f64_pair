// -*- coding: utf-8 -*-
// Marc B. Reynolds, 2020-2025
// Public Domain under http://unlicense.org, see link for details.

// The tests for all functions that involve add/sub is about worthless
// ATM. In general none the testing is making effort to examine any
// hard cases. So it's really just weak spot checking for the moment.

#include "common.h"

#define TRIALS 0x40000

#ifndef TRIALS
#define TRIALS 0x02800000
#endif

// Since "fr" routines don't normalize the result "exact" and "computed" pairs can appear quite different.
// (specifically displayed values aren't normalized either)



// globals 
mpfr_t mp_r0;
mpfr_t mp_r1;
mpfr_t mp_a;
mpfr_t mp_b;
mpfr_t mp_c;
mpfr_t mp_d;
mpfr_t mp_e;
mpfr_t mp_t;


void report_ulp(report_tbuild_t* b, void* data)
{
  report_table_col_t* c = b->col;

  double d = ((double*)data)[0];

  if (d < 0x1.0p40)
    report_table_format(b,"%*.*f",c->digits,c->prec,d);
  else {
    report_table_format(b,SGR_RGB(255,150,150) "%*.*fu" SGR_RESET,c->digits,c->prec,d*0x1.0p-53);
  }
}

#define REPORT_ULP(N,D,P)          \
  .name   = N,                     \
  .type   = report_ttype_ptr,      \
  .pformat= report_ulp,            \
  .digits = D,                     \
  .prec   = P,                     \
  .extra  = D+P+2



static inline int mpfr_inv(mpfr_t r ,const mpfr_t x, mpfr_rnd_t rnd)
{
  return mpfr_ui_div(r,1,x,rnd);
}


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
  bool      updated;
  
} uop_d_table_t;



uop_d_table_t op_d[] =
{
//{ DEF_FE(fe_rsqrt_d,  mpfr_rec_sqrt) }  <- the OP_U(...) part can be deleted like this ("tracked" peak error info across runs)
  { DEF_FE(fe_rsqrt_d,  mpfr_rec_sqrt), OP_U(0x1.573b74p+2, 0x1.faea8d04b1668p-2, UP(0x1.6bd97ca99753cp+0, 0x1.acd1ede5c7de7p-53)) },
  { DEF_FE(fe_rsqrt_dh, mpfr_rec_sqrt), OP_U(0x1.a56438p+1, 0x1.9c1b95807da1ep-2, UP(0x1.9389c166b0628p+0, 0x1.96f00fc997be8p-53)) },
  { DEF_FE(fe_inv_d,    mpfr_inv),      OP_U(0x1p-1,        0x1.56f581ba03a0ap-2, UP(0x1.7e2e0633c775ap+1, 0x1.d6872e3e1af22p-53)) },
  { DEF_FE(fe_inv_dn,   mpfr_inv),      OP_U(0x1.7f683p+0,  0x1.ff5b5fd4f00d4p-1, UP(0x1.00526a957415dp+0,-0x1.ffcdfa3432c6ep-54)) },
  { DEF_FE(fe_inv_dh,   mpfr_inv),      OP_U(0x1.34c4ep+0,  0x1.6a324624fbf52p-1, UP(0x1.69e18b2b0b124p+0, 0x1.fff6dd32b13f2p-54)) },
  { DEF_FR(fr_sqrt_d,   mpfr_sqrt),     OP_U(0x1.ff159p-1,  0x1.001b080141b3p-2,  UP(0x1.000d83a54f9b3p-1,-0x1.ff523f8237456p-55)) }
};


// test unary f(double) → pair
void op_d_tests(void)
{
  printf(SGR_BOLD SGR_RGB(200,200,255) "\nf(double) → pair : x ∈ [0,1) \n" SGR_RESET);
  
  report_table_t table = {
    .col = {
      { REPORT_TABLE_STR("f(x)",14), .just=report_table_justify_left },
      { REPORT_TABLE_POS_F("ulp",1,4) },
      { REPORT_TABLE_A64("x") },
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
      { REPORT_TABLE_A64("e.hi") },      
      { REPORT_TABLE_A64("e.lo") },      
    }
  };

  report_table_t htable = {
    .col = {
      { REPORT_TABLE_STR("f(x)",14), .just=report_table_justify_left },
      { REPORT_TABLE_POS_F("ulp",1,4) },
      { REPORT_TABLE_A64("x") },
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") }
    }
  };


  // historic data
  report_table_header(stdout, &htable);

  for(size_t i=0; i<LENGTHOF(op_d); i++) {
    uop_d_table_t* t = op_d + i;

    if (t->max_ulp != 0.0) {
      fe_pair_t (*f)(double) = t->fe;
      fe_pair_t r = f(t->max_x);
      
      if (!fe_eq(r, t->max_r)) {
        htable.eol = " ← disagrees : bug, change of implementation, etc.\n";
      }

      report_table_row(stdout,&htable, t->name, t->max_ulp, t->max_x,
                       t->max_r.hi,
                       t->max_r.lo
                       );

      htable.eol = "\n";
    }
  }
  report_table_end(stdout, &htable);


  // current run
  report_table_header(stdout, &table);
  
  for(size_t i=0; i<LENGTHOF(op_d); i++) {
    uop_d_table_t* t = op_d + i;
    
    double    max_ulp = 0.0;
    double    max_x;
    fe_pair_t max_r;
    fe_pair_t max_e;
    
    fe_pair_t (*f)(double) = t->fe;
    int       (*mp)(mpfr_t,const mpfr_t,mpfr_rnd_t) = t->mp;
    
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
      max_x   = a;
      max_r   = r;
      max_e   = e;
    }


    // print-out current
    report_table_row(stdout,&table,
                     t->name,
                     max_ulp,
                     max_x,
                     max_r.hi,max_r.lo,
                     max_e.hi,max_e.lo
                     );

    // update table
    if (max_ulp > t->max_ulp) {
      t->max_ulp = max_ulp;
      t->max_x   = max_x;
      t->max_r   = max_r;
      t->max_e   = max_e;
      t->updated = true;
    }
  }

  
  
  report_table_end(stdout, &table);
}


//**********************************************************
// unary ops (pair input)

typedef struct {
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
  bool      updated;
  
} uop_table_t;


// for inspecting just forwarding to mul.
fe_pair_t fe_sq_x(fe_pair_t x) { return fe_mul(x,x); }
fr_pair_t fr_sq_x(fr_pair_t x) { return fr_mul(x,x); }

// slightly better than default
fe_pair_t fe_rsqrt_x2(fe_pair_t x) {  return fe_sqrt(fe_inv_a(x)); }


uop_table_t uops[] =
{
//{ DEF_FE(fe_sq_x,  mpfr_sqr), OP_U(0x1.bfeac4p+1, UP(0x1.a38d4bd4d9cfcp+0, 0x1.ffcec18bf0e71p-54), UP(0x1.57cbea1baff64p+1, 0x1.48180fd1708p-59))   },
//{ DEF_FR(fr_sq_x,  mpfr_sqr), OP_U(0x1.ff88ccp+1, UP(0x1.00675d6df320ap+0, 0x1.ff12c8585be9cp-54), UP(0x1.00cee49836d3fp+0, 0x1.61e6e0a963d54p-52)) },

  { DEF_FE(fe_rsqrt,   mpfr_rec_sqrt), OP_U(0x1.393f6cp+2, UP(0x1.07225df4049fp+0, 0x1.f8ee93c787c06p-54), UP(0x1.f902f0cff5919p-1,-0x1.19d7b8458b996p-54)) },
  { DEF_FE(fe_rsqrt_s, mpfr_rec_sqrt), OP_U(0x1.ce7fcep+3, UP(0x1.039ecc881c261p+0,0x1.f36e6a19afc5ep-54), UP(0x1.fc6aeaa5de476p-1, 0x1.978a88841f2bp-56)) },
//{ DEF_FE(fe_rsqrt_x2,mpfr_rec_sqrt), OP_U(0x1.167444p+2, UP(0x1.12de2a9524ebp+0, 0x1.fab11e58faad7p-54), UP(0x1.ee1d6c6191ddcp-1,-0x1.0d22e17fa9d1bp-54))},

  { DEF_FE(fe_sq,    mpfr_sqr), OP_U(0x1.341ee8p+1, UP(0x1.69c9710dacf63p+0, 0x1.49456cb73be15p-54), UP(0x1.ff49bf5b4b351p+0,-0x1.193c992e56cp-60))   },
  { DEF_FE(fe_sq_hq, mpfr_sqr), OP_U(0x1.8p+0,      UP(0x1.8005230146627p+0, 0x1.467ae867a8928p-54), UP(0x1.2007b48f1afecp+1,-0x1.e6477cfc4fb7p-54))  }, 
  { DEF_FR(fr_sq,    mpfr_sqr), OP_U(0x1.7ff614p+1, UP(0x1.229486b38d6b2p+0, 0x1.fffd8c93f7beap-54), UP(0x1.49d4d75ad2e2bp+0, 0x1.0816997c3baeap-52)) },

  { DEF_FE(fe_inv,   mpfr_inv), OP_U(0x1.b3e456p+2, UP(0x1.03f98822d44b1p+0, 0x1.fdf06a56693fbp-54), UP(0x1.f82c0ce4c8549p-1,-0x1.dacf304ef2868p-55)) },
  { DEF_FE(fe_inv_n, mpfr_inv), OP_U(0x1.b5e106p+2, UP(0x1.00d0047006f1dp+0, 0x1.fb9d588e7e10ap-54), UP(0x1.fe61481c8b467p-1,-0x1.e8df70bf85325p-55)) },
  { DEF_FE(fe_inv,   mpfr_inv), OP_U(0x1.b3e456p+2, UP(0x1.03f98822d44b1p+0, 0x1.fdf06a56693fbp-54), UP(0x1.f82c0ce4c8549p-1,-0x1.dacf304ef2868p-55)) },
  { DEF_FE(fe_inv_n, mpfr_inv), OP_U(0x1.b5e106p+2, UP(0x1.00d0047006f1dp+0, 0x1.fb9d588e7e10ap-54), UP(0x1.fe61481c8b467p-1,-0x1.e8df70bf85325p-55)) },
  { DEF_FR(fr_inv,   mpfr_inv), OP_U(0x1.eebddap+2, UP(0x1.01674824684e7p+0, 0x1.f7291b1be27fp-54), UP(0x1.fd355aae68ef6p-1,-0x1.71c6f1a803bc1p-53)) },
  { DEF_FR(fr_inv_n, mpfr_inv), OP_U(0x1.085e49p+3, UP(0x1.00048b09bc88cp+0, 0x1.fe98c86e04422p-54), UP(0x1.fff6ea15cdd5fp-1,-0x1.7e0c5705437e4p-53)) },
  { DEF_FR(fr_inv_a, mpfr_inv), OP_U(0x1.b2cf6p+2,  UP(0x1.01c7416a5f0d3p+0, 0x1.fddd0b6ac709ep-54), UP(0x1.fc77c532d3367p-1,-0x1.7af51f66ed10bp-53)) },
  { DEF_FE(fe_sqrt,  mpfr_sqrt),OP_U(0x1.8da008p+1, UP(0x1.00ec6c9d80937p+0, 0x1.f9b76feee2a06p-54), UP(0x1.00761b10455b8p+0, 0x1.7d7d284f7655ap-53)) },
};

// test unary
void op_p_tests(void)
{
  printf(SGR_BOLD SGR_RGB(200,200,255) "\nf(pair) → pair : x ∈ [1,2)\n" SGR_RESET);

  report_table_t table = {
    .col = {
      { REPORT_TABLE_STR("f(x)",12), .just=report_table_justify_left },
      { REPORT_ULP("ulp",4,4) },
      { REPORT_TABLE_A64("x.hi") },      
      { REPORT_TABLE_A64("x.lo") },
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
      { REPORT_TABLE_A64("e.hi") },      
      { REPORT_TABLE_A64("e.lo") },      
    }
  };


  report_table_t htable = {
    .col = {
      { REPORT_TABLE_STR("f(x)",12), .just=report_table_justify_left },
      { REPORT_ULP("ulp",4,4) },
      { REPORT_TABLE_A64("x.hi") },      
      { REPORT_TABLE_A64("x.lo") },
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") }
    }
  };


  // historic data
  report_table_header(stdout, &htable);

  for(size_t i=0; i<LENGTHOF(uops); i++) {
    uop_table_t* t = uops + i;

    if (t->max_ulp != 0.0) {
      fe_pair_t (*f)(fe_pair_t) = t->fe;
      fe_pair_t r = f(t->max_x);
      
      if (!fe_eq(r, t->max_r)) {
        htable.eol = " ← disagrees : bug, change of implementation, etc.\n";
      }

      report_table_row(stdout,&htable, t->name, &(t->max_ulp),
                       t->max_x.hi, t->max_x.lo, 
                       t->max_r.hi, t->max_r.lo
                       );

      htable.eol = "\n";
    }
  }
  report_table_end(stdout, &htable);


  // current scan
  report_table_header(stdout, &table);
  
  for(size_t i=0; i<LENGTHOF(uops); i++) {
    uop_table_t* t = uops + i;
    double    max_ulp = -2.0;
    fe_pair_t max_a;
    fe_pair_t max_r;
    fe_pair_t max_e;
    
    fe_pair_t (*f)(fe_pair_t) = t->fe;
    int       (*mp)(mpfr_t,const mpfr_t,mpfr_rnd_t) = t->mp;
    
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

    report_table_row(stdout,&table,
                     t->name,
                     &max_ulp,
                     max_a,
                     max_r.hi,max_r.lo,
                     max_e.hi,max_e.lo
                     );
    
    if (max_ulp > t->max_ulp) {
      t->max_ulp = max_ulp;
      t->max_x   = max_a;
      t->max_r   = max_r;
      t->max_e   = max_e;
      t->updated = true;
    }
  }
  
  report_table_end(stdout, &table);

  // dump any update peak values
  for(size_t i=0; i<LENGTHOF(uops); i++) {
    uop_table_t* t = uops + i;
    if (t->updated) {
      printf("%13s : OP_U(%a, UP(%a,% a), UP(%a,% a))\n", t->name,
             t->max_ulp,
             t->max_x.hi,t->max_x.lo,
             t->max_r.hi,t->max_r.lo);
    }
  }
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
  bool      updated;
  
} uop_dd_table_t;

uop_dd_table_t op_dd[] =
{
  { DEF_FE(fe_div_dd,  mpfr_div), OP_B(0x1p-1, 0x1.f53551678cfa2p-1,0x1.b634d5858bc6p-2,  UP(0x1.24ce410baddbbp+1,-0x1.117aebad2f645p-53)) },
  { DEF_FR(fr_div_dd,  mpfr_div), OP_B(0x1p-1, 0x1.e1c5056b83d64p-1,0x1.cdc3a0984a9bdp-1, UP(0x1.0b174727e9f41p+0, 0x1.e5a07f110338cp-54)) },
};


void op_dd_tests(void)
{
  printf(SGR_BOLD SGR_RGB(200,200,255) "\nf(double,double) → pair\n" SGR_RESET);
  
  report_table_t table = {
    .col = {
      { REPORT_TABLE_STR("f(a,b)",14), .just=report_table_justify_left },
      { REPORT_TABLE_POS_F("ulp",1,4) },
      { REPORT_TABLE_A64("a") },
      { REPORT_TABLE_A64("b") },
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
      { REPORT_TABLE_A64("e.hi") },      
      { REPORT_TABLE_A64("e.lo") },      
    }
  };

  report_table_t htable = {
    .col = {
      { REPORT_TABLE_STR("f(a,b)",14), .just=report_table_justify_left },
      { REPORT_TABLE_POS_F("ulp",1,4) },
      { REPORT_TABLE_A64("a") },
      { REPORT_TABLE_A64("b") },
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
    }
  };


  // historic data
  report_table_header(stdout, &htable);

  for(size_t i=0; i<LENGTHOF(op_dd); i++) {
    uop_dd_table_t* t = op_dd + i;

    if (t->max_ulp != 0.0) {
      fe_pair_t (*f)(double,double) = t->fe;
      fe_pair_t r = f(t->max_a,t->max_b);
      
      if (!fe_eq(r, t->max_r)) {
        htable.eol = " ← disagrees : bug, change of implementation, etc.\n";
      }

      report_table_row(stdout,&htable, t->name, t->max_ulp,
                       t->max_a, t->max_b, 
                       t->max_r.hi, t->max_r.lo
                       );

      htable.eol = "\n";
    }
  }
  report_table_end(stdout, &htable);

  // current scan
  report_table_header(stdout, &table);
  
  for(size_t i=0; i<LENGTHOF(op_dd); i++) {
    uop_dd_table_t* t = op_dd + i;
    double    max_ulp = 0.0;
    double    max_a;
    double    max_b;
    fe_pair_t max_r;
    fe_pair_t max_e;
    
    fe_pair_t (*f)(double,double) = t->fe;
    int       (*mp)(mpfr_t,const mpfr_t,const mpfr_t,mpfr_rnd_t) = t->mp;
    
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
    
    report_table_row(stdout,&table,
                     t->name,
                     max_ulp,
                     max_a,
                     max_b,
                     max_r.hi,max_r.lo,
                     max_e.hi,max_e.lo
                     );

    if (max_ulp > t->max_ulp) {
      t->max_ulp = max_ulp;
      t->max_a   = max_a;
      t->max_b   = max_b;
      t->max_r   = max_r;
      t->max_e   = max_e;
      t->updated = true;
    }    
  }
  
  report_table_end(stdout, &table);


  // dump any update peak values
  for(size_t i=0; i<LENGTHOF(op_dd); i++) {
    uop_dd_table_t* t = op_dd + i;
    if (t->updated) {
      printf("%13s : OP_B(%a, %a,%a, UP(%a,% a))\n", t->name,
             t->max_ulp,
             t->max_a, t->max_b,
             t->max_r.hi,t->max_r.lo);
    }
  }
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
  bool      updated;
  
} uop_dp_table_t;

uop_dp_table_t op_dp[] =
{
  { DEF_FE(fe_d_add, mpfr_add), OP_B(0x1p+0,        0x1.c2f8982c34d9p-3,  UP(0x1.a70c5cf1fde5dp-3, 0x1.5519c72ad6036p-57), UP(0x1.b5027a8f195f7p-2,-0x1.55731c6a94fe4p-56)) },
  { DEF_FE(fe_d_div, mpfr_div), OP_B(0x1.bc3698p+2, 0x1.f135dea7a280ap-1, UP(0x1.04fa52e0932f6p-2, 0x1.f416a52268bdep-56), UP(0x1.e7ba0cde25a4cp+1,-0x1.a1dc3508d491cp-53)) }
};

void op_dp_tests(void)
{
  printf(SGR_BOLD SGR_RGB(200,200,255) "\nf(double,pair) → pair\n" SGR_RESET);
  
  report_table_t table = {
    .col = {
      { REPORT_TABLE_STR("f(a,b)",14), .just=report_table_justify_left },
      { REPORT_TABLE_POS_F("ulp",1,4) },
      { REPORT_TABLE_A64("a") },
      { REPORT_TABLE_A64("b.hi") },      
      { REPORT_TABLE_A64("b.lo") },      
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
      { REPORT_TABLE_A64("e.hi") },      
      { REPORT_TABLE_A64("e.lo") },      
    }
  };

  report_table_t htable = {
    .col = {
      { REPORT_TABLE_STR("f(a,b)",14), .just=report_table_justify_left },
      { REPORT_TABLE_POS_F("ulp",1,4) },
      { REPORT_TABLE_A64("a") },
      { REPORT_TABLE_A64("b.hi") },      
      { REPORT_TABLE_A64("b.lo") },      
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
    }
  };

  // historic data
  report_table_header(stdout, &htable);

  for(size_t i=0; i<LENGTHOF(op_dp); i++) {
    uop_dp_table_t* t = op_dp + i;

    if (t->max_ulp != 0.0) {
      fe_pair_t (*f)(double,fe_pair_t) = t->fe;
      fe_pair_t r = f(t->max_a,t->max_b);
      
      if (!fe_eq(r, t->max_r)) {
        htable.eol = " ← disagrees : bug, change of implementation, etc.\n";
      }

      report_table_row(stdout,&htable, t->name, t->max_ulp,
                       t->max_a,
                       t->max_b.hi, t->max_b.lo, 
                       t->max_r.hi, t->max_r.lo
                       );

      htable.eol = "\n";
    }
  }

  report_table_end(stdout, &htable);

  // current scan
  report_table_header(stdout, &table);
  
  for(size_t i=0; i<LENGTHOF(op_dp); i++) {
    uop_dp_table_t* t = op_dp + i;
    double    max_ulp = 0.0;
    double    max_a;
    fe_pair_t max_b;
    fe_pair_t max_r;
    fe_pair_t max_e;
    
    fe_pair_t (*f)(double,fe_pair_t) = t->fe;
    int       (*mp)(mpfr_t,const mpfr_t,const mpfr_t,mpfr_rnd_t) = t->mp;
    
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

    report_table_row(stdout,&table,
                     t->name,
                     max_ulp,
                     max_a,
                     max_b.hi,max_b.lo,
                     max_r.hi,max_r.lo,
                     max_e.hi,max_e.lo);
    
    if (max_ulp > t->max_ulp) {
      t->max_ulp = max_ulp;
      t->max_a   = max_a;
      t->max_b   = max_b;
      t->max_r   = max_r;
      t->max_e   = max_e;
      t->updated = true;
    }
  }
  
  report_table_end(stdout, &table);

  // dump any update peak values
  for(size_t i=0; i<LENGTHOF(op_dp); i++) {
    uop_dp_table_t* t = op_dp + i;
    if (t->updated) {
      printf("%13s : OP_B(%a, %a, UP(%a,% a), UP(%a,% a))\n", t->name,
             t->max_ulp,
             t->max_a,
             t->max_b.hi,t->max_b.lo,
             t->max_r.hi,t->max_r.lo);
    }
  }
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
  bool      updated;
  
} uop_pd_table_t;

uop_pd_table_t op_pd[] =
{
  { DEF_FE(fe_mul_d, mpfr_mul), OP_B(0x1.fffffp+0,  UP(0x1.33e294d2806dbp-3, 0x1.c0a8d58c8151bp-57), 0x1.7b59decb9a5bp-3,  UP(0x1.c83c880cbcaf6p-6, 0x1.720b8574d91dp-62)) },
  { DEF_FR(fr_mul_d, mpfr_mul), OP_B(0x1.fffff8p+0, UP(0x1.3fd7c01f52d96p+0, 0x1.cfc2a62bae943p-54), 0x1.67da020fc674ap-1, UP(0x1.c197eeb8ecf27p-1, 0x1.1a9016ff1f84dp-53)) },
  { DEF_FE(fe_div_d, mpfr_div), OP_B(0x1.fee34p+1,  UP(0x1.002fc2a81154cp-1, 0x1.bdb64b52e4b97p-55), 0x1.010e69edf6804p-2, UP(0x1.fe4485e67ddfdp+0, 0x1.5bfc9502acp-56))  },
  { DEF_FR(fr_add_d, mpfr_add), OP_B(0x1p+0,        UP(0x1.1d8dd71bd7e35p+0, 0x1.a50834c49ddd1p-54), 0x1.91bf667159818p-4, UP(0x1.36a9cd82ed7b6p+0, 0x1.d2841a624eee8p-53)) },
  { DEF_FE(fe_add_d, mpfr_add), OP_B(0x1p+0,        UP(0x1.9fa24e02a90c9p-1, 0x1.1beddffa97ee6p-55), 0x1.c05242834fcfp-1,  UP(0x1.affa4842fc6ddp+0,-0x1.72091002b408cp-54)) },
  { DEF_FE(fe_add_d_cr, mpfr_add), OP_B(0x1p-1, UP(0x1.897cf19a37e37p-2, 0x1.7697d03ace3f1p-56), 0x1.6530d2ec293fp-2,  UP(0x1.7756e24330914p-1,-0x1.44b417e298e08p-55)) },
};

void op_pd_tests(void)
{
  printf(SGR_BOLD SGR_RGB(200,200,255) "\nf(pair,double) → pair\n" SGR_RESET);
  
  report_table_t table = {
    .col = {
      { REPORT_TABLE_STR("f(a,b)",14), .just=report_table_justify_left },
      { REPORT_TABLE_POS_F("ulp",1,4) },
      { REPORT_TABLE_A64("a.hi") },
      { REPORT_TABLE_A64("a.lo") },      
      { REPORT_TABLE_A64("b") },
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
      { REPORT_TABLE_A64("e.hi") },      
      { REPORT_TABLE_A64("e.lo") },      
    }
  };

  report_table_t htable = {
    .col = {
      { REPORT_TABLE_STR("f(a,b)",14), .just=report_table_justify_left },
      { REPORT_TABLE_POS_F("ulp",1,4) },
      { REPORT_TABLE_A64("a.hi") },
      { REPORT_TABLE_A64("a.lo") },      
      { REPORT_TABLE_A64("b") },
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
    }
  };

  // historic data
  report_table_header(stdout, &htable);

  for(size_t i=0; i<LENGTHOF(op_pd); i++) {
    uop_pd_table_t* t = op_pd + i;

    if (t->max_ulp != 0.0) {
      fe_pair_t (*f)(fe_pair_t,double) = t->fe;
      fe_pair_t r = f(t->max_a,t->max_b);
      
      if (!fe_eq(r, t->max_r)) {
        htable.eol = " ← disagrees : bug, change of implementation, etc.\n";
      }

      report_table_row(stdout,&htable, t->name, t->max_ulp,
                       t->max_a.hi, t->max_a.lo,
                       t->max_b, 
                       t->max_r.hi, t->max_r.lo
                       );

      htable.eol = "\n";
    }

  }
  report_table_end(stdout, &htable);

  
  // current scan
  report_table_header(stdout, &table);
  
  for(size_t i=0; i<LENGTHOF(op_pd); i++) {
    uop_pd_table_t* t = op_pd + i;
    double    max_ulp = 0.0;
    fe_pair_t max_a;
    double    max_b;
    fe_pair_t max_r;
    fe_pair_t max_e;
    
    fe_pair_t (*f)(fe_pair_t,double) = t->fe;
    int       (*mp)(mpfr_t,const mpfr_t,const mpfr_t,mpfr_rnd_t) = t->mp;
    
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
    
    report_table_row(stdout,&table,
                     t->name,
                     max_ulp,
                     max_a.hi,max_a.lo,
                     max_b,
                     max_r.hi,max_r.lo,
                     max_e.hi,max_e.lo
                     );

    if (max_ulp > t->max_ulp) {
      t->max_ulp = max_ulp;
      t->max_a   = max_a;
      t->max_b   = max_b;
      t->max_r   = max_r;
      t->max_e   = max_e;
      t->updated = true;
    }
  }
  
  report_table_end(stdout, &table);

  // dump any update peak values
  for(size_t i=0; i<LENGTHOF(op_pd); i++) {
    uop_pd_table_t* t = op_pd + i;
    if (t->updated) {
      printf("%13s : OP_B(%a, UP(%a,% a), %a, UP(%a,% a))\n", t->name,
             t->max_ulp,
             t->max_a.hi,t->max_a.lo,
             t->max_b,
             t->max_r.hi,t->max_r.lo);
    }
  }
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

  bool      updated;
  
} uop_pp_table_t;

uop_pp_table_t op_pp[] =
{
  { DEF_FE(fe_mul,   mpfr_mul), OP_B(0x1.ffdadcp+1, UP(0x1.33a2475df0007p+0, 0x1.4113fc8ecc3f8p-54), UP(0x1.58aaa545d661p+0,  0x1.e51907ceacb49p-54), UP(0x1.9e2f2045233dp+0,  0x1.b3b3c8fe40e08p-54)) },
  { DEF_FR(fr_mul,   mpfr_mul), OP_B(0x1.3e82e6p+2, UP(0x1.2ebcf4375c469p+0, 0x1.ffd4e75eccd61p-54), UP(0x1.771d5c40a1d91p+0, 0x1.fdd05c6ffb02p-54),  UP(0x1.bb99823d18d75p+0, 0x1.f7fbe914166d2p-53)) },
  { DEF_FE(fe_div,   mpfr_div), OP_B(0x1.f6445ap+2, UP(0x1.f46ff175bf788p+0, 0x1.255fb8eee3b9ap-54), UP(0x1.03907994c26d3p+0, 0x1.f68b4e34ab4fap-54), UP(0x1.ed90b38de2a11p+0,-0x1.af5971ef3c02p-56))  },
  { DEF_FR(fr_div,   mpfr_div), OP_B(0x1.bbb4bcp+2, UP(0x1.e36632b496997p+0, 0x1.338b14743e68p-54),  UP(0x1.004144deee7dfp-7, 0x1.ed0f1804acb91p-61), UP(0x1.e2eb1319c23dap+7,-0x1.0911f15d0370bp-45)) },
  { DEF_FR(fr_div_a, mpfr_div), OP_B(0x1.8246d6p+2, UP(0x1.f67db1f73460bp+0, 0x1.3f6e78a6d70a8p-54), UP(0x1.015016a836f95p-2, 0x1.fb4801df60511p-56), UP(0x1.f3ed5e2cb47b4p+2,-0x1.0fbc6656c0e2ep-50)) },
  { DEF_FE(fe_sub,   mpfr_sub), OP_B(0x1.4p+1,      UP(0x1.70c95f88e9523p+0, 0x1.fda6293f1b28ep-54), UP(0x1.c72b0e7b45b27p-2, 0x1.4567c8464fa93p-56), UP(0x1.fdfd37d42fcb3p-1, 0x1.58986e5b0e7dp-55)) },
  { DEF_FR(fr_sub,   mpfr_sub), OP_B(0x1p+13,       UP(0x1.0004f6d1bdea2p+0, 0x1.e849f520c9091p-54), UP(0x1.fffee0564a72ap-1, 0x1.afef4ffc80111p-55), UP(0x1.61a9a62c34p-14,   0x1.10524d2289008p-54)) },
  { DEF_FE(fe_osub,  mpfr_sub), OP_B(0x1.4p+1,      UP(0x1.088f69877af96p+0, 0x1.b30dbba2b9ff7p-54), UP(0x1.f3aa597ee7f76p-3, 0x1.94aab92331f02p-57), UP(0x1.94343caf3bf4fp-1, 0x1.00f0c8fca783p-55)) },
  { DEF_FR(fr_osub,  mpfr_sub), OP_B(0x1p+12,       UP(0x1.0007e27c8cddbp+0, 0x1.b4e36aaac5a9cp-54), UP(0x1.fffc1612fefc3p-1, 0x1.48404b8fc8e1fp-55), UP(0x1.3aee61abf3p-13,   0x1.10c344e2e138cp-54)) },
  { DEF_FE(fe_add,   mpfr_add), OP_B(0x1p+1,        UP(0x1.1244ecc34f519p+0, 0x1.b92a3357d3067p-54), UP(0x1.3db617b2ecdcbp-1, 0x1.e980b7ce4e1dfp-55), UP(0x1.b11ff89cc5bffp+0, 0x1.5bd51e7df42b1p-55)) },
  { DEF_FR(fr_add,   mpfr_add), OP_B(0x1.7ep+1,     UP(0x1.2f63c832f7e45p+0, 0x1.ff63065548ca8p-54), UP(0x1.0c300bdba276p-6,  0x1.ed8f996e5a2bfp-60), UP(0x1.33948862666e2p+0, 0x1.01c6512ec08ccp-52)) },
  { DEF_FE(fe_add_s, mpfr_add), OP_B(0x1.7fp+1,     UP(0x1.acfcbd2d17b7ap+0, 0x1.fd1e2873d4f7fp-54), UP(0x1.ec4ef7bb0ec4p-7,  0x1.af3b010488701p-61), UP(0x1.b0d55b1c8dd53p+0, 0x1.f279d77824p-64)) },
  { DEF_FE(fe_sub_s, mpfr_sub), OP_B(0x1p+11,       UP(0x1.00115f257f459p+0, 0x1.b30314678e636p-54), UP(0x1.ffeb8f2bd131dp-1, 0x1.649846149b5c7p-55), UP(0x1.b978f96acac03p-12,-0x1.21d457e95cp-67)) },
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

  report_table_t htable = {
    .col = {
      { REPORT_TABLE_STR("f(a,b)",14), .just=report_table_justify_left },
      { REPORT_TABLE_F("ulp",8,4) },
      { REPORT_TABLE_A64("a.hi") },
      { REPORT_TABLE_A64("a.lo") },      
      { REPORT_TABLE_A64("b.hi") },
      { REPORT_TABLE_A64("b.lo") },      
      { REPORT_TABLE_A64("r.hi") },      
      { REPORT_TABLE_A64("r.lo") },      
    }
  };

  // historic data
  report_table_header(stdout, &htable);
  
  for(size_t i=0; i<LENGTHOF(op_pp); i++) {
    uop_pp_table_t* t = op_pp + i;
    
    if (t->max_ulp != 0.0) {
      fe_pair_t (*f)(fe_pair_t,fe_pair_t) = t->fe;
      fe_pair_t r = f(t->max_a,t->max_b);
      
      if (!fe_eq(r, t->max_r)) {
        htable.eol = " ← disagrees : bug, change of implementation, etc.\n";
      }
      
      report_table_row(stdout,&htable, t->name, t->max_ulp,
                       t->max_a.hi, t->max_a.lo, 
                       t->max_b.hi, t->max_b.lo, 
                       t->max_r.hi, t->max_r.lo
                       );

      htable.eol = "\n";
    }
  }
  report_table_end(stdout, &htable);
  
  // current scan
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
        
    report_table_row(stdout,&table,
                     op_pp[i].name,
                     max_ulp,
                     max_a.hi,max_a.lo,
                     //max_b,
                     max_b.hi,max_b.lo,
                     max_r.hi,max_r.lo,
                     max_e.hi,max_e.lo
                     );

    if (max_ulp > op_pp[i].max_ulp) {
      op_pp[i].max_ulp = max_ulp;
      op_pp[i].max_a   = max_a;
      op_pp[i].max_b   = max_b;
      op_pp[i].max_r   = max_r;
      op_pp[i].max_e   = max_e;
      op_pp[i].updated = true;
    }
  }
  
  report_table_end(stdout, &table);

  // update any updated peaks
  for(size_t i=0; i<LENGTHOF(op_pp); i++) {
    uop_pp_table_t* t = op_pp + i;
    if (t->updated) {
      printf("%13s : OP_B(%a, UP(%a,% a), UP(%a,% a), UP(%a,% a))\n",
             t->name,
             t->max_ulp,
             t->max_a.hi, t->max_a.lo,
             t->max_b.hi, t->max_b.lo,
             t->max_r.hi, t->max_r.lo);
    }
  }
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
  mpfr_init2(mp_c,  128);
  mpfr_init2(mp_d,  128);

  mpfr_init2(mp_e,  128);
  mpfr_init2(mp_t,  128);
  
#if 0
  //op_d_tests();
  op_p_tests();
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



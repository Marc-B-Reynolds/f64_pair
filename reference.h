// -*- coding: utf-8 -*-

// paper or other reference-y versions of routines that are modified
// in f64_pair.h


// Boldo & Melquiond (2008) round-to-odd
static inline double bm_sum_ro(double x, double y)
{
  fe_pair_t v = fe_two_sum(x,y);
  uint64_t  r = fe_to_bits(v.hi);

  // expectation: probability of being odd/even is
  // equal and exact is about zero.  sadface for
  // branchy code.
  if ((r & 1) == 0 && (v.lo != 0.0)) {
    if ((v.hi > 0.0) ^ (v.lo < 0.0))
      r++;
    else
      r--;
  }

  return fe_from_bits(r);
}


// Boldo & Melquiond (2008) round-to-odd
static inline double bm_add3(double a, double b, double c)
{
  fe_pair_t x = fe_two_sum(a,b);
  fe_pair_t s = fe_two_sum(c,x.hi);
  double    v = bm_sum_ro(x.lo,s.lo);

  return(s.hi+v);
}

// Graillat & Muller (2025)
static inline int gm_not_1or3_times_pot(double x)
{
  static const double P = 2251799813685249.0, Q = 2251799813685248.0;
  double Delta;
  Delta = (P*x)-(Q*x);
  return (Delta != x);
}

// Graillat & Muller (2025)
double gm_add3(double a, double b, double c)
{
  fe_pair_t x = fe_two_sum(a,b);
  fe_pair_t s = fe_two_sum(x.hi,c);
  fe_pair_t v = fe_two_sum(x.lo,s.lo);

  if ((gm_not_1or3_times_pot(v.hi)) || (v.lo == 0))
    return s.hi+v.hi;

  if ((v.lo < 0) ^ (v.hi < 0))
    return s.hi+(0.875*v.hi);
  else
    return s.hi+(1.125*v.hi);
}

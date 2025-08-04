// -*- coding: utf-8 -*-

// paper or other reference-y versions of routines that are modified
// in f64_pair.h

static inline int paper_not_1or3_time_pot(double x)
{
  static const double P = 2251799813685249.0, Q = 2251799813685248.0;
  double Delta;
  Delta = (P*x)-(Q*x);
  return (Delta != x);
}


static inline double sum_ro_f64_(double x, double y)
{
  fe_pair_t v = fe_two_sum(x,y);
  uint64_t  r = type_pun(v.hi,uint64_t);

  // expectation: probability of being odd/even is
  // equal and exact is about zero.  sadface for
  // branchy code.
  if ((r & 1) == 0 && (v.lo != 0.0)) {
    if ((v.hi > 0.0) ^ (v.lo < 0.0))
      r++;
    else
      r--;
  }

  return type_pun(r,double);
}


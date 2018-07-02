/*
 * MgSO4_EOS2_planetary_smaller.cpp
 *
 * Code generation for function 'MgSO4_EOS2_planetary_smaller'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MgSO4_EOS2_planetary_smaller.h"
#include "interp3.h"

/* Function Declarations */
static int div_s32(int numerator, int denominator);

/* Function Definitions */
static int div_s32(int numerator, int denominator)
{
  int quotient;
  unsigned int absNumerator;
  unsigned int absDenominator;
  boolean_T quotientNeedsNegation;
  if (denominator == 0) {
    if (numerator >= 0) {
      quotient = MAX_int32_T;
    } else {
      quotient = MIN_int32_T;
    }
  } else {
    if (numerator < 0) {
      absNumerator = ~(unsigned int)numerator + 1U;
    } else {
      absNumerator = (unsigned int)numerator;
    }

    if (denominator < 0) {
      absDenominator = ~(unsigned int)denominator + 1U;
    } else {
      absDenominator = (unsigned int)denominator;
    }

    quotientNeedsNegation = ((numerator < 0) != (denominator < 0));
    absNumerator /= absDenominator;
    if (quotientNeedsNegation) {
      quotient = -(int)absNumerator;
    } else {
      quotient = (int)absNumerator;
    }
  }

  return quotient;
}

void MgSO4_EOS2_planetary_smaller(double m, double P, double T, double *rho,
  double *Cp, double *alpha)
{
  double Pgg[4200];
  int jtilecol;
  double mgg[4200];
  int ibtile;
  int ia;
  double y[15];
  int ibtmp;
  static const double dv0[15] = { -20.0, -11.428571428571429,
    -2.8571428571428577, 5.7142857142857153, 14.285714285714285,
    22.857142857142854, 31.428571428571431, 40.0, 48.571428571428569,
    57.142857142857139, 65.714285714285708, 74.285714285714292,
    82.857142857142861, 91.428571428571431, 100.0 };

  static const double a[14] = { 0.0, 0.19230769230769232, 0.38461538461538464,
    0.57692307692307687, 0.76923076923076927, 0.96153846153846156,
    1.1538461538461537, 1.3461538461538463, 1.5384615384615385,
    1.7307692307692308, 1.9230769230769231, 2.1153846153846154,
    2.3076923076923075, 2.5 };

  int q;
  int db[3];
  static const double dv1[20] = { 0.0, 42.10526315789474, 84.21052631578948,
    126.31578947368421, 168.42105263157896, 210.52631578947367,
    252.63157894736841, 294.73684210526318, 336.84210526315792,
    378.94736842105266, 421.05263157894734, 463.15789473684208,
    505.26315789473682, 547.36842105263156, 589.47368421052636, 631.578947368421,
    673.68421052631584, 715.78947368421052, 757.89473684210532, 800.0 };

  double Tgg[4200];
  static const signed char iv0[3] = { 14, 20, 15 };

  static const signed char iv1[3] = { 1, 1, 15 };

  /* ,Videal,Vex,Vw,VDH */
  /* function [rho,Vs,vel,Cp,alpha,Videal,Vex]=MgSO4_EOS2(m,P,T) */
  /*  usage: */
  /*   [rho,Vs,Videal,Vex,Vw,vel,Cp,alpha]=MgSO4_EOS(m,P,T,parms) */
  /*    units: P in MPa, T in �C, volumes in cc/g Cp in J/kg/C, rho in gm/cc */
  /*    velocities in km/s */
  /*  set up the EOS on the parameter grid */
  /*  interpolate the parameters onto the grid of the input */
  for (jtilecol = 0; jtilecol < 15; jtilecol++) {
    ibtile = jtilecol * 280;
    for (ia = 0; ia < 20; ia++) {
      ibtmp = ibtile + ia * 14;
      for (q = 0; q < 14; q++) {
        Pgg[ibtmp + q] = dv1[ia];
      }
    }
  }

  for (jtilecol = 0; jtilecol < 300; jtilecol++) {
    ibtile = jtilecol * 14;
    memcpy(&mgg[ibtile], &a[0], 14U * sizeof(double));
  }

  memcpy(&y[0], &dv0[0], 15U * sizeof(double));
  for (jtilecol = 0; jtilecol < 3; jtilecol++) {
    db[jtilecol] = 1;
  }

  for (jtilecol = 0; jtilecol < 2; jtilecol++) {
    db[jtilecol + 1] = db[jtilecol] * iv0[jtilecol];
  }

  for (ibtile = 0; ibtile < 4200; ibtile++) {
    ia = 0;
    ibtmp = ibtile;
    for (jtilecol = 2; jtilecol >= 0; jtilecol += -1) {
      q = div_s32(div_s32(ibtmp, db[jtilecol]) * db[jtilecol], db[jtilecol]);
      ia = (ia + q) - div_s32(q, (int)iv1[jtilecol]) * iv1[jtilecol];
      ibtmp -= div_s32(ibtmp, db[jtilecol]) * db[jtilecol];
    }

    Tgg[ibtile] = y[ia];
  }

  /* VDH=interp3(Pgg,mgg,Tgg,VDHg,Pc,mc,Tc,'spline'); */
  /* Vex=interp3(Pgg,mgg,Tgg,Vexg,Pc,mc,Tc,'spline'); */
  /* Videal=interp3(Pgg,mgg,Tgg,Videalg,Pc,mc,Tc,'spline'); */
  /* Vs=interp3(Pgg,mgg,Tgg,Vsg,Pc,mc,Tc,'spline'); */
  *rho = interp3(Pgg, mgg, Tgg, P, m, T);

  /*  vel=interp3(Pgg,mgg,Tgg,velg,Pc,mc,Tc,'spline'); */
  *alpha = b_interp3(Pgg, mgg, Tgg, P, m, T);
  *Cp = c_interp3(Pgg, mgg, Tgg, P, m, T);
}

/* End of code generation (MgSO4_EOS2_planetary_smaller.cpp) */

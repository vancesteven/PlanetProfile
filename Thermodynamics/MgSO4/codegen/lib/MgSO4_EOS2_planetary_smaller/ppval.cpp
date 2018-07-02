/*
 * ppval.cpp
 *
 * Code generation for function 'ppval'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MgSO4_EOS2_planetary_smaller.h"
#include "ppval.h"

/* Function Definitions */
void b_ppval(const double pp_breaks[20], const double pp_coefs[1064], double x,
             double v[14])
{
  int low_i;
  int low_ip1;
  int mid_i;
  int high_i;
  int icp;
  double xloc;
  if (rtIsNaN(x)) {
    for (mid_i = 0; mid_i < 14; mid_i++) {
      v[mid_i] = x;
    }
  } else {
    low_i = 1;
    low_ip1 = 2;
    high_i = 20;
    while (high_i > low_ip1) {
      mid_i = (low_i + high_i) >> 1;
      if (x >= pp_breaks[mid_i - 1]) {
        low_i = mid_i;
        low_ip1 = mid_i + 1;
      } else {
        high_i = mid_i;
      }
    }

    icp = (low_i - 1) * 14;
    xloc = x - pp_breaks[low_i - 1];
    memcpy(&v[0], &pp_coefs[icp], 14U * sizeof(double));
    for (low_ip1 = 0; low_ip1 < 3; low_ip1++) {
      high_i = icp + (low_ip1 + 1) * 266;
      for (mid_i = 0; mid_i < 14; mid_i++) {
        v[mid_i] = xloc * v[mid_i] + pp_coefs[high_i + mid_i];
      }
    }
  }
}

double c_ppval(const double pp_breaks[14], const double pp_coefs[52], double x)
{
  double v;
  int low_i;
  int low_ip1;
  int high_i;
  double xloc;
  int mid_i;
  if (rtIsNaN(x)) {
    v = x;
  } else {
    low_i = 0;
    low_ip1 = 2;
    high_i = 14;
    while (high_i > low_ip1) {
      mid_i = ((low_i + high_i) + 1) >> 1;
      if (x >= pp_breaks[mid_i - 1]) {
        low_i = mid_i - 1;
        low_ip1 = mid_i + 1;
      } else {
        high_i = mid_i;
      }
    }

    xloc = x - pp_breaks[low_i];
    v = pp_coefs[low_i];
    for (low_ip1 = 0; low_ip1 < 3; low_ip1++) {
      v = xloc * v + pp_coefs[low_i + (low_ip1 + 1) * 13];
    }
  }

  return v;
}

void ppval(const double pp_breaks[15], const double pp_coefs[15680], double x,
           double v[280])
{
  int low_i;
  int low_ip1;
  int mid_i;
  int high_i;
  int icp;
  double xloc;
  if (rtIsNaN(x)) {
    for (mid_i = 0; mid_i < 280; mid_i++) {
      v[mid_i] = x;
    }
  } else {
    low_i = 1;
    low_ip1 = 2;
    high_i = 15;
    while (high_i > low_ip1) {
      mid_i = (low_i + high_i) >> 1;
      if (x >= pp_breaks[mid_i - 1]) {
        low_i = mid_i;
        low_ip1 = mid_i + 1;
      } else {
        high_i = mid_i;
      }
    }

    icp = (low_i - 1) * 280;
    xloc = x - pp_breaks[low_i - 1];
    memcpy(&v[0], &pp_coefs[icp], 280U * sizeof(double));
    for (low_ip1 = 0; low_ip1 < 3; low_ip1++) {
      high_i = icp + (low_ip1 + 1) * 3920;
      for (mid_i = 0; mid_i < 280; mid_i++) {
        v[mid_i] = xloc * v[mid_i] + pp_coefs[high_i + mid_i];
      }
    }
  }
}

/* End of code generation (ppval.cpp) */

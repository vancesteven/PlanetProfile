/*
 * spline.cpp
 *
 * Code generation for function 'spline'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MgSO4_EOS2_planetary_smaller.h"
#include "spline.h"
#include "pwchcore.h"

/* Function Definitions */
void b_spline(const double x[20], const double y[280], double output_breaks[20],
              double output_coefs[1064])
{
  double dx[19];
  double dvdf[266];
  double s[280];
  int k;
  int pgm1;
  double d31;
  int pg;
  int pgp1;
  double dnnm2;
  int j;
  double d1;
  double md[20];
  for (k = 0; k < 19; k++) {
    dx[k] = x[k + 1] - x[k];
    pgm1 = k * 14;
    pgp1 = (k + 1) * 14;
    for (j = 0; j < 14; j++) {
      dvdf[pgm1 + j] = (y[pgp1 + j] - y[pgm1 + j]) / dx[k];
    }
  }

  for (k = 0; k < 18; k++) {
    pg = (k + 1) * 14;
    pgm1 = k * 14;
    for (j = 0; j < 14; j++) {
      s[pg + j] = 3.0 * (dx[k + 1] * dvdf[pgm1 + j] + dx[k] * dvdf[pg + j]);
    }
  }

  d31 = x[2] - x[0];
  dnnm2 = x[19] - x[17];
  d1 = dx[0];
  for (j = 0; j < 14; j++) {
    s[j] = ((d1 + 2.0 * d31) * dx[1] * dvdf[j] + d1 * d1 * dvdf[j + 14]) / d31;
  }

  d1 = dx[18];
  for (j = 0; j < 14; j++) {
    s[j + 266] = ((d1 + 2.0 * dnnm2) * dx[17] * dvdf[j + 252] + d1 * d1 * dvdf[j
                  + 238]) / dnnm2;
  }

  md[0] = dx[1];
  md[19] = dx[17];
  for (k = 0; k < 18; k++) {
    md[k + 1] = 2.0 * (dx[k + 1] + dx[k]);
  }

  d1 = dx[1] / md[0];
  md[1] -= d1 * d31;
  for (j = 0; j < 14; j++) {
    s[j + 14] -= d1 * s[j];
  }

  for (k = 0; k < 17; k++) {
    d1 = dx[k + 2] / md[k + 1];
    md[k + 2] -= d1 * dx[k];
    pg = (k + 2) * 14;
    pgm1 = (k + 1) * 14;
    for (j = 0; j < 14; j++) {
      s[pg + j] -= d1 * s[pgm1 + j];
    }
  }

  d1 = dnnm2 / md[18];
  md[19] -= d1 * dx[17];
  for (j = 0; j < 14; j++) {
    s[j + 266] -= d1 * s[j + 252];
    s[j + 266] /= md[19];
  }

  for (k = 17; k >= 0; k += -1) {
    pg = (k + 1) * 14;
    pgp1 = (k + 2) * 14;
    for (j = 0; j < 14; j++) {
      s[pg + j] = (s[pg + j] - dx[k] * s[pgp1 + j]) / md[k + 1];
    }
  }

  for (j = 0; j < 14; j++) {
    s[j] = (s[j] - d31 * s[j + 14]) / md[0];
  }

  b_pwchcore(x, y, s, dx, dvdf, output_breaks, output_coefs);
}

void c_spline(const double x[14], const double y[14], double output_breaks[14],
              double output_coefs[52])
{
  double dvdf[13];
  double s[14];
  double dx[13];
  int k;
  double d31;
  double dnnm2;
  double md[14];
  double r;
  for (k = 0; k < 13; k++) {
    d31 = x[k + 1] - x[k];
    dvdf[k] = (y[k + 1] - y[k]) / d31;
    dx[k] = d31;
  }

  d31 = x[2] - x[0];
  dnnm2 = x[13] - x[11];
  s[0] = ((dx[0] + 2.0 * d31) * dx[1] * dvdf[0] + dx[0] * dx[0] * dvdf[1]) / d31;
  s[13] = ((dx[12] + 2.0 * dnnm2) * dx[11] * dvdf[12] + dx[12] * dx[12] * dvdf
           [11]) / dnnm2;
  md[0] = dx[1];
  md[13] = dx[11];
  for (k = 0; k < 12; k++) {
    s[k + 1] = 3.0 * (dx[k + 1] * dvdf[k] + dx[k] * dvdf[k + 1]);
    md[k + 1] = 2.0 * (dx[k + 1] + dx[k]);
  }

  r = dx[1] / md[0];
  md[1] -= r * d31;
  s[1] -= r * s[0];
  for (k = 0; k < 11; k++) {
    r = dx[k + 2] / md[k + 1];
    md[k + 2] -= r * dx[k];
    s[k + 2] -= r * s[k + 1];
  }

  r = dnnm2 / md[12];
  md[13] -= r * dx[11];
  s[13] -= r * s[12];
  s[13] /= md[13];
  for (k = 11; k >= 0; k += -1) {
    s[k + 1] = (s[k + 1] - dx[k] * s[k + 2]) / md[k + 1];
  }

  s[0] = (s[0] - d31 * s[1]) / md[0];
  memcpy(&output_breaks[0], &x[0], 14U * sizeof(double));
  for (k = 0; k < 13; k++) {
    d31 = (dvdf[k] - s[k]) / dx[k];
    dnnm2 = (s[k + 1] - dvdf[k]) / dx[k];
    output_coefs[k] = (dnnm2 - d31) / dx[k];
    output_coefs[k + 13] = 2.0 * d31 - dnnm2;
    output_coefs[k + 26] = s[k];
    output_coefs[k + 39] = y[k];
  }
}

void spline(const double x[15], const double y[4200], double output_breaks[15],
            double output_coefs[15680])
{
  double dx[14];
  double dvdf[3920];
  double s[4200];
  int k;
  int pgm1;
  double d31;
  int pg;
  int pgp1;
  double dnnm2;
  int j;
  double d1;
  double md[15];
  for (k = 0; k < 14; k++) {
    dx[k] = x[k + 1] - x[k];
    pgm1 = k * 280;
    pgp1 = (k + 1) * 280;
    for (j = 0; j < 280; j++) {
      dvdf[pgm1 + j] = (y[pgp1 + j] - y[pgm1 + j]) / dx[k];
    }
  }

  for (k = 0; k < 13; k++) {
    pg = (k + 1) * 280;
    pgm1 = k * 280;
    for (j = 0; j < 280; j++) {
      s[pg + j] = 3.0 * (dx[k + 1] * dvdf[pgm1 + j] + dx[k] * dvdf[pg + j]);
    }
  }

  d31 = x[2] - x[0];
  dnnm2 = x[14] - x[12];
  d1 = dx[0];
  for (j = 0; j < 280; j++) {
    s[j] = ((d1 + 2.0 * d31) * dx[1] * dvdf[j] + d1 * d1 * dvdf[j + 280]) / d31;
  }

  d1 = dx[13];
  for (j = 0; j < 280; j++) {
    s[j + 3920] = ((d1 + 2.0 * dnnm2) * dx[12] * dvdf[j + 3640] + d1 * d1 *
                   dvdf[j + 3360]) / dnnm2;
  }

  md[0] = dx[1];
  md[14] = dx[12];
  for (k = 0; k < 13; k++) {
    md[k + 1] = 2.0 * (dx[k + 1] + dx[k]);
  }

  d1 = dx[1] / md[0];
  md[1] -= d1 * d31;
  for (j = 0; j < 280; j++) {
    s[j + 280] -= d1 * s[j];
  }

  for (k = 0; k < 12; k++) {
    d1 = dx[k + 2] / md[k + 1];
    md[k + 2] -= d1 * dx[k];
    pg = (k + 2) * 280;
    pgm1 = (k + 1) * 280;
    for (j = 0; j < 280; j++) {
      s[pg + j] -= d1 * s[pgm1 + j];
    }
  }

  d1 = dnnm2 / md[13];
  md[14] -= d1 * dx[12];
  for (j = 0; j < 280; j++) {
    s[j + 3920] -= d1 * s[j + 3640];
    s[j + 3920] /= md[14];
  }

  for (k = 12; k >= 0; k += -1) {
    pg = (k + 1) * 280;
    pgp1 = (k + 2) * 280;
    for (j = 0; j < 280; j++) {
      s[pg + j] = (s[pg + j] - dx[k] * s[pgp1 + j]) / md[k + 1];
    }
  }

  for (j = 0; j < 280; j++) {
    s[j] = (s[j] - d31 * s[j + 280]) / md[0];
  }

  pwchcore(x, y, s, dx, dvdf, output_breaks, output_coefs);
}

/* End of code generation (spline.cpp) */

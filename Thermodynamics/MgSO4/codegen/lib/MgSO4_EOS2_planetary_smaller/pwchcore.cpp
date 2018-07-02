/*
 * pwchcore.cpp
 *
 * Code generation for function 'pwchcore'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "MgSO4_EOS2_planetary_smaller.h"
#include "pwchcore.h"

/* Function Definitions */
void b_pwchcore(const double x[20], const double y[280], const double s[280],
                const double dx[19], const double divdif[266], double pp_breaks
                [20], double pp_coefs[1064])
{
  int j;
  double dxj;
  int joffset;
  int i;
  double dzzdx;
  double dzdxdx;
  memcpy(&pp_breaks[0], &x[0], 20U * sizeof(double));
  for (j = 0; j < 19; j++) {
    dxj = dx[j];
    joffset = j * 14;
    for (i = 0; i < 14; i++) {
      dzzdx = (divdif[joffset + i] - s[joffset + i]) / dxj;
      dzdxdx = (s[(joffset + i) + 14] - divdif[joffset + i]) / dxj;
      pp_coefs[joffset + i] = (dzdxdx - dzzdx) / dxj;
      pp_coefs[(joffset + i) + 266] = 2.0 * dzzdx - dzdxdx;
      pp_coefs[(joffset + i) + 532] = s[joffset + i];
      pp_coefs[(joffset + i) + 798] = y[joffset + i];
    }
  }
}

void pwchcore(const double x[15], const double y[4200], const double s[4200],
              const double dx[14], const double divdif[3920], double pp_breaks
              [15], double pp_coefs[15680])
{
  int j;
  double dxj;
  int joffset;
  int i;
  double dzzdx;
  double dzdxdx;
  memcpy(&pp_breaks[0], &x[0], 15U * sizeof(double));
  for (j = 0; j < 14; j++) {
    dxj = dx[j];
    joffset = j * 280;
    for (i = 0; i < 280; i++) {
      dzzdx = (divdif[joffset + i] - s[joffset + i]) / dxj;
      dzdxdx = (s[(joffset + i) + 280] - divdif[joffset + i]) / dxj;
      pp_coefs[joffset + i] = (dzdxdx - dzzdx) / dxj;
      pp_coefs[(joffset + i) + 3920] = 2.0 * dzzdx - dzdxdx;
      pp_coefs[(joffset + i) + 7840] = s[joffset + i];
      pp_coefs[(joffset + i) + 11760] = y[joffset + i];
    }
  }
}

/* End of code generation (pwchcore.cpp) */

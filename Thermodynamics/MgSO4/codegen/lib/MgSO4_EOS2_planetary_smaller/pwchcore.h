/*
 * pwchcore.h
 *
 * Code generation for function 'pwchcore'
 *
 */

#ifndef PWCHCORE_H
#define PWCHCORE_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "MgSO4_EOS2_planetary_smaller_types.h"

/* Function Declarations */
extern void b_pwchcore(const double x[20], const double y[280], const double s
  [280], const double dx[19], const double divdif[266], double pp_breaks[20],
  double pp_coefs[1064]);
extern void pwchcore(const double x[15], const double y[4200], const double s
                     [4200], const double dx[14], const double divdif[3920],
                     double pp_breaks[15], double pp_coefs[15680]);

#endif

/* End of code generation (pwchcore.h) */

/*
 * ppval.h
 *
 * Code generation for function 'ppval'
 *
 */

#ifndef PPVAL_H
#define PPVAL_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "MgSO4_EOS2_planetary_smaller_types.h"

/* Function Declarations */
extern void b_ppval(const double pp_breaks[20], const double pp_coefs[1064],
                    double x, double v[14]);
extern double c_ppval(const double pp_breaks[14], const double pp_coefs[52],
                      double x);
extern void ppval(const double pp_breaks[15], const double pp_coefs[15680],
                  double x, double v[280]);

#endif

/* End of code generation (ppval.h) */

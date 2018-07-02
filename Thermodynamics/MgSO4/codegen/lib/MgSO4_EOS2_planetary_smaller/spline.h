/*
 * spline.h
 *
 * Code generation for function 'spline'
 *
 */

#ifndef SPLINE_H
#define SPLINE_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "MgSO4_EOS2_planetary_smaller_types.h"

/* Function Declarations */
extern void b_spline(const double x[20], const double y[280], double
                     output_breaks[20], double output_coefs[1064]);
extern void c_spline(const double x[14], const double y[14], double
                     output_breaks[14], double output_coefs[52]);
extern void spline(const double x[15], const double y[4200], double
                   output_breaks[15], double output_coefs[15680]);

#endif

/* End of code generation (spline.h) */

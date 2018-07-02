/*
 * interp3.h
 *
 * Code generation for function 'interp3'
 *
 */

#ifndef INTERP3_H
#define INTERP3_H

/* Include files */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "MgSO4_EOS2_planetary_smaller_types.h"

/* Function Declarations */
extern double b_interp3(const double varargin_1[4200], const double varargin_2
  [4200], const double varargin_3[4200], double varargin_5, double varargin_6,
  double varargin_7);
extern double c_interp3(const double varargin_1[4200], const double varargin_2
  [4200], const double varargin_3[4200], double varargin_5, double varargin_6,
  double varargin_7);
extern double interp3(const double varargin_1[4200], const double varargin_2
                      [4200], const double varargin_3[4200], double varargin_5,
                      double varargin_6, double varargin_7);

#endif

/* End of code generation (interp3.h) */

/*
 * main.cpp
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include files */
#include "rt_nonfinite.h"
#include "MgSO4_EOS2_planetary_smaller.h"
#include "main.h"
#include "MgSO4_EOS2_planetary_smaller_terminate.h"
#include "MgSO4_EOS2_planetary_smaller_initialize.h"

/* Function Declarations */
static double argInit_real_T();
static void c_main_MgSO4_EOS2_planetary_sma();

/* Function Definitions */
static double argInit_real_T()
{
  return 0.0;
}

static void c_main_MgSO4_EOS2_planetary_sma()
{
  double rho;
  double Cp;
  double alpha;

  /* Initialize function 'MgSO4_EOS2_planetary_smaller' input arguments. */
  /* Call the entry-point 'MgSO4_EOS2_planetary_smaller'. */
  MgSO4_EOS2_planetary_smaller(argInit_real_T(), argInit_real_T(),
    argInit_real_T(), &rho, &Cp, &alpha);
}

int main(int, const char * const [])
{
  /* Initialize the application.
     You do not need to do this more than one time. */
  MgSO4_EOS2_planetary_smaller_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  c_main_MgSO4_EOS2_planetary_sma();

  /* Terminate the application.
     You do not need to do this more than one time. */
  MgSO4_EOS2_planetary_smaller_terminate();
  return 0;
}

/* End of code generation (main.cpp) */

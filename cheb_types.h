/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * cheb_types.h
 *
 * Code generation for function 'cheb'
 *
 */

#ifndef CHEB_TYPES_H
#define CHEB_TYPES_H

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif

/* End of code generation (cheb_types.h) */

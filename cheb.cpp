/*
 *
 * Code generation for function 'cheb', This code used to a MATLAB code
 *
 * This files contains function that generates a Chebyshev mnatrix
 *
 *
 */

/* Include files */
#include <cmath>
#include <math.h>
#include "rt_nonfinite.h"
#include "cheb.h"
#include "cheb_emxutil.h"

/* Function Declarations */
static double rt_powd_snf(double u0, double u1);

/* Function Definitions */
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d0;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d0 = std::abs(u0);
    d1 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (d0 == 1.0) {
        y = 1.0;
      } else if (d0 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

void cheb(double N, emxArray_real_T *D, emxArray_real_T *x)
{
  emxArray_real_T *y;
  int k;
  int ibtile;
  int npages;
  int nx;
  emxArray_real_T *b_y;
  emxArray_real_T *r0;
  emxArray_real_T *c;
  emxArray_real_T *X;
  int xpageoffset;
  double t;
  emxArray_real_T *r1;

  /*  CHEB compute D = differentiation matrix, x = Chebyshev grid */
  if (N == 0.0) {
    k = D->size[0] * D->size[1];
    D->size[0] = 1;
    D->size[1] = 1;
    emxEnsureCapacity_real_T(D, k);
    D->data[0] = 0.0;
    k = x->size[0];
    x->size[0] = 1;
    emxEnsureCapacity_real_T(x, k);
    x->data[0] = 1.0;
  } else {
    emxInit_real_T(&y, 2);
    if (rtIsNaN(N)) {
      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity_real_T(y, k);
      y->data[0] = rtNaN;
    } else if (N < 0.0) {
      y->size[0] = 1;
      y->size[1] = 0;
    } else if (rtIsInf(N) && (0.0 == N)) {
      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity_real_T(y, k);
      y->data[0] = rtNaN;
    } else {
      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      ibtile = static_cast<int>(std::floor(N));
      y->size[1] = ibtile + 1;
      emxEnsureCapacity_real_T(y, k);
      for (k = 0; k <= ibtile; k++) {
        y->data[k] = k;
      }
    }

    k = y->size[0] * y->size[1];
    npages = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_real_T(y, npages);
    ibtile = k - 1;
    for (k = 0; k <= ibtile; k++) {
      y->data[k] = 3.1415926535897931 * y->data[k] / N;
    }

    nx = y->size[1];
    for (k = 0; k < nx; k++) {
      y->data[k] = std::cos(y->data[k]);
    }

    k = x->size[0];
    x->size[0] = y->size[1];
    emxEnsureCapacity_real_T(x, k);
    ibtile = y->size[1];
    for (k = 0; k < ibtile; k++) {
      x->data[k] = y->data[k];
    }

    if (rtIsNaN(N)) {
      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity_real_T(y, k);
      y->data[0] = rtNaN;
    } else if (N < 0.0) {
      y->size[0] = 1;
      y->size[1] = 0;
    } else if (rtIsInf(N) && (0.0 == N)) {
      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity_real_T(y, k);
      y->data[0] = rtNaN;
    } else {
      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = static_cast<int>(std::floor(N)) + 1;
      emxEnsureCapacity_real_T(y, k);
      ibtile = static_cast<int>(std::floor(N));
      for (k = 0; k <= ibtile; k++) {
        y->data[k] = k;
      }
    }

    emxInit_real_T(&b_y, 2);
    k = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = y->size[1];
    emxEnsureCapacity_real_T(b_y, k);
    nx = y->size[1];
    for (k = 0; k < nx; k++) {
      b_y->data[k] = rt_powd_snf(-1.0, y->data[k]);
    }

    emxInit_real_T(&r0, 1);
    k = r0->size[0];
    ibtile = static_cast<int>((N - 1.0));
    r0->size[0] = 2 + ibtile;
    emxEnsureCapacity_real_T(r0, k);
    r0->data[0] = 2.0;
    for (k = 0; k < ibtile; k++) {
      r0->data[k + 1] = 1.0;
    }

    emxInit_real_T(&c, 1);
    r0->data[1 + ibtile] = 2.0;
    k = c->size[0];
    c->size[0] = r0->size[0];
    emxEnsureCapacity_real_T(c, k);
    ibtile = r0->size[0];
    for (k = 0; k < ibtile; k++) {
      c->data[k] = r0->data[k] * b_y->data[k];
    }

    emxFree_real_T(&r0);
    emxFree_real_T(&b_y);
    emxInit_real_T(&X, 2);
    nx = x->size[0];
    k = X->size[0] * X->size[1];
    X->size[0] = nx;
    npages = static_cast<int>((N + 1.0));
    X->size[1] = npages;
    emxEnsureCapacity_real_T(X, k);
    nx = x->size[0];
    for (xpageoffset = 0; xpageoffset < npages; xpageoffset++) {
      ibtile = xpageoffset * nx;
      for (k = 0; k < nx; k++) {
        X->data[ibtile + k] = x->data[k];
      }
    }

    if (N + 1.0 < 0.0) {
      t = 0.0;
    } else {
      t = N + 1.0;
    }

    emxInit_real_T(&r1, 2);
    nx = static_cast<int>(t);
    k = r1->size[0] * r1->size[1];
    r1->size[0] = static_cast<int>(t);
    r1->size[1] = static_cast<int>(t);
    emxEnsureCapacity_real_T(r1, k);
    ibtile = static_cast<int>(t) * static_cast<int>(t);
    for (k = 0; k < ibtile; k++) {
      r1->data[k] = 0.0;
    }

    if (static_cast<int>(t) > 0) {
      for (k = 0; k < nx; k++) {
        r1->data[k + r1->size[0] * k] = 1.0;
      }
    }

    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = c->size[0];
    emxEnsureCapacity_real_T(y, k);
    ibtile = c->size[0];
    for (k = 0; k < ibtile; k++) {
      y->data[k] = 1.0 / c->data[k];
    }

    k = D->size[0] * D->size[1];
    D->size[0] = c->size[0];
    D->size[1] = y->size[1];
    emxEnsureCapacity_real_T(D, k);
    ibtile = c->size[0];
    for (k = 0; k < ibtile; k++) {
      nx = y->size[1];
      for (npages = 0; npages < nx; npages++) {
        t = c->data[k] * y->data[npages];
        D->data[k + D->size[0] * npages] = t / ((X->data[k + X->size[0] * npages]
          - X->data[npages + X->size[0] * k]) + r1->data[k + r1->size[0] *
          npages]);
      }
    }

    emxFree_real_T(&c);

    /* off-diagonal entries */
    k = X->size[0] * X->size[1];
    X->size[0] = D->size[1];
    X->size[1] = D->size[0];
    emxEnsureCapacity_real_T(X, k);
    ibtile = D->size[0];
    for (k = 0; k < ibtile; k++) {
      nx = D->size[1];
      for (npages = 0; npages < nx; npages++) {
        X->data[npages + X->size[0] * k] = D->data[k + D->size[0] * npages];
      }
    }

    ibtile = X->size[0];
    npages = X->size[1];
    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = X->size[1];
    emxEnsureCapacity_real_T(y, k);
    for (nx = 0; nx < npages; nx++) {
      xpageoffset = nx * X->size[0];
      y->data[nx] = X->data[xpageoffset];
      for (k = 2; k <= ibtile; k++) {
        y->data[nx] += X->data[(xpageoffset + k) - 1];
      }
    }

    emxFree_real_T(&X);
    npages = y->size[1];
    nx = y->size[1];
    xpageoffset = y->size[1];
    k = r1->size[0] * r1->size[1];
    r1->size[0] = nx;
    r1->size[1] = xpageoffset;
    emxEnsureCapacity_real_T(r1, k);
    ibtile = nx * xpageoffset;
    for (k = 0; k < ibtile; k++) {
      r1->data[k] = 0.0;
    }

    for (nx = 0; nx < npages; nx++) {
      r1->data[nx + r1->size[0] * nx] = y->data[nx];
    }

    emxFree_real_T(&y);
    k = D->size[0] * D->size[1];
    npages = D->size[0] * D->size[1];
    emxEnsureCapacity_real_T(D, npages);
    ibtile = k - 1;
    for (k = 0; k <= ibtile; k++) {
      D->data[k] -= r1->data[k];
    }

    emxFree_real_T(&r1);
  }

  /* diagonal entries */
}

/* End of code generation (cheb.cpp) */

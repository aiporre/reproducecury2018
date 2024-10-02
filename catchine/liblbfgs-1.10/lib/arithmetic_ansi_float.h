/*
 *      ANSI C implementation of vecfloattor operations.
 *
 * Copyright (c) 2007-2010 Naoaki Okazaki
 * All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* $Id$ */

#include <stdlib.h>
#include <memory.h>

#if     LBFGSFLOAT_IEEE_FLOAT
#define fsigndiff_float(x, y) (((*(uint32_t*)(x)) ^ (*(uint32_t*)(y))) & 0x80000000U)
#else
#define fsigndiff_float(x, y) (*(x) * (*(y) / fabs(*(y))) < 0.)
#endif/*LBFGS_IEEE_FLOAT*/

inline static void* vecfloatalloc(size_t size)
{
    void *memblock = malloc(size);
    if (memblock) {
        memset(memblock, 0, size);
    }
    return memblock;
}

inline static void vecfloatfree(void *memblock)
{
    free(memblock);
}

inline static void vecfloatset(lbfgsfloatfloatval_t *x, const lbfgsfloatfloatval_t c, const int n)
{
    int i;
    
    for (i = 0;i < n;++i) {
        x[i] = c;
    }
}

inline static void vecfloatcpy(lbfgsfloatfloatval_t *y, const lbfgsfloatfloatval_t *x, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] = x[i];
    }
}

inline static void vecfloatncpy(lbfgsfloatfloatval_t *y, const lbfgsfloatfloatval_t *x, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] = -x[i];
    }
}

inline static void vecfloatadd(lbfgsfloatfloatval_t *y, const lbfgsfloatfloatval_t *x, const lbfgsfloatfloatval_t c, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] += c * x[i];
    }
}

inline static void vecfloatdiff(lbfgsfloatfloatval_t *z, const lbfgsfloatfloatval_t *x, const lbfgsfloatfloatval_t *y, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        z[i] = x[i] - y[i];
    }
}

inline static void vecfloatscale(lbfgsfloatfloatval_t *y, const lbfgsfloatfloatval_t c, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] *= c;
    }
}

inline static void vecfloatmul(lbfgsfloatfloatval_t *y, const lbfgsfloatfloatval_t *x, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] *= x[i];
    }
}

inline static void vecfloatdot(float* s, const lbfgsfloatfloatval_t *x, const lbfgsfloatfloatval_t *y, const int n)
{
    int i;
    *s = 0.;
    for (i = 0;i < n;++i) {
        *s += x[i] * y[i];
    }
}

inline static void vecfloat2norm(float* s, const lbfgsfloatfloatval_t *x, const int n)
{
    vecfloatdot(s, x, x, n);
    *s = (float)sqrt(*s);
}

inline static void vecfloat2norminv(float* s, const lbfgsfloatfloatval_t *x, const int n)
{
    vecfloat2norm(s, x, n);
    *s = (float)(1.0 / *s);
}

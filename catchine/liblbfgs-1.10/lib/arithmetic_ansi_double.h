/*
 *      ANSI C implementation of vecdoubletor operations.
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

#define fsigndiff_double(x, y) (*(x) * (*(y) / fabs(*(y))) < 0.)

inline static void* vecdoublealloc(size_t size)
{
    void *memblock = malloc(size);
    if (memblock) {
        memset(memblock, 0, size);
    }
    return memblock;
}

inline static void vecdoublefree(void *memblock)
{
    free(memblock);
}

inline static void vecdoubleset(lbfgsdoublefloatval_t *x, const lbfgsdoublefloatval_t c, const int n)
{
    int i;
    
    for (i = 0;i < n;++i) {
        x[i] = c;
    }
}

inline static void vecdoublecpy(lbfgsdoublefloatval_t *y, const lbfgsdoublefloatval_t *x, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] = x[i];
    }
}

inline static void vecdoublencpy(lbfgsdoublefloatval_t *y, const lbfgsdoublefloatval_t *x, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] = -x[i];
    }
}

inline static void vecdoubleadd(lbfgsdoublefloatval_t *y, const lbfgsdoublefloatval_t *x, const lbfgsdoublefloatval_t c, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] += c * x[i];
    }
}

inline static void vecdoublediff(lbfgsdoublefloatval_t *z, const lbfgsdoublefloatval_t *x, const lbfgsdoublefloatval_t *y, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        z[i] = x[i] - y[i];
    }
}

inline static void vecdoublescale(lbfgsdoublefloatval_t *y, const lbfgsdoublefloatval_t c, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] *= c;
    }
}

inline static void vecdoublemul(lbfgsdoublefloatval_t *y, const lbfgsdoublefloatval_t *x, const int n)
{
    int i;

    for (i = 0;i < n;++i) {
        y[i] *= x[i];
    }
}

inline static void vecdoubledot(lbfgsdoublefloatval_t* s, const lbfgsdoublefloatval_t *x, const lbfgsdoublefloatval_t *y, const int n)
{
    int i;
    *s = 0.;
    for (i = 0;i < n;++i) {
        *s += x[i] * y[i];
    }
}

inline static void vecdouble2norm(lbfgsdoublefloatval_t* s, const lbfgsdoublefloatval_t *x, const int n)
{
    vecdoubledot(s, x, x, n);
    *s = (lbfgsdoublefloatval_t)sqrt(*s);
}

inline static void vecdouble2norminv(lbfgsdoublefloatval_t* s, const lbfgsdoublefloatval_t *x, const int n)
{
    vecdouble2norm(s, x, n);
    *s = (lbfgsdoublefloatval_t)(1.0 / *s);
}

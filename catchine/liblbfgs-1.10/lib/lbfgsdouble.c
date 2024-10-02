/*
 *      Limited memory BFGS (L-BFGS).
 *
 * Copyright (c) 1990, Jorge Nocedal
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

/*
This library is a C port of the FORTRAN implementation of Limited-memory
Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) method written by Jorge Nocedal.
The original FORTRAN source code is available at:
http://www.ece.northwestern.edu/~nocedal/lbfgsdouble.html

The L-BFGS algorithm is described in:
    - Jorge Nocedal.
      Updating Quasi-Newton Matrices with Limited Storage.
      <i>Mathematics of Computation</i>, Vol. 35, No. 151, pp. 773--782, 1980.
    - Dong C. Liu and Jorge Nocedal.
      On the limited memory BFGS method for large scale optimization.
      <i>Mathematical Programming</i> B, Vol. 45, No. 3, pp. 503-528, 1989.

The line search algorithms used in this implementation are described in:
    - John E. Dennis and Robert B. Schnabel.
      <i>Numerical Methods for Unconstrained Optimization and Nonlinear
      Equations</i>, Englewood Cliffs, 1983.
    - Jorge J. More and David J. Thuente.
      Line search algorithm with guaranteed sufficient decrease.
      <i>ACM Transactions on Mathematical Software (TOMS)</i>, Vol. 20, No. 3,
      pp. 286-307, 1994.

This library also implements Orthant-Wise Limited-memory Quasi-Newton (OWL-QN)
method presented in:
    - Galen Andrew and Jianfeng Gao.
      Scalable training of L1-regularized log-linear models.
      In <i>Proceedings of the 24th International Conference on Machine
      Learning (ICML 2007)</i>, pp. 33-40, 2007.

I would like to thank the original author, Jorge Nocedal, who has been
distributing the effieicnt and explanatory implementation in an open source
licence.
*/

#ifdef  HAVE_CONFIG_H
#include <config.h>
#endif/*HAVE_CONFIG_H*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <lbfgsdouble.h>

#ifdef  _MSC_VER
#define inline  __inline
#endif/*_MSC_VER*/

#if     defined(USE_SSE) && defined(__SSE2__) && LBFGSDOUBLE_FLOAT == 64
/* Use SSE2 optimization for 64bit double precision. */
#include "arithmetic_sse_double.h"

#elif   defined(USE_SSE) && defined(__SSE__) && LBFGSDOUBLE_FLOAT == 32
/* Use SSE optimization for 32bit float precision. */
#include "arithmetic_sse_float.h"

#else
/* No CPU specific optimization. */
#include "arithmetic_ansi_double.h"

#endif

#define min2(a, b)      ((a) <= (b) ? (a) : (b))
#define max2(a, b)      ((a) >= (b) ? (a) : (b))
#define max3(a, b, c)   max2(max2((a), (b)), (c));

struct tag_callback_double_data {
    int n;
    void *instance;
    lbfgsdouble_evaluate_t proc_evaluate;
    lbfgsdouble_progress_t proc_progress;
};
typedef struct tag_callback_double_data callback_double_data_t;

struct tag_iteration_double_data {
    lbfgsdoublefloatval_t alpha;
    lbfgsdoublefloatval_t *s;     /* [n] */
    lbfgsdoublefloatval_t *y;     /* [n] */
    lbfgsdoublefloatval_t ys;     /* vecdoubledot(y, s) */
};
typedef struct tag_iteration_double_data iteration_double_data_t;

static const lbfgs_parameter_t _defparam_double = {
    6, 1e-5, 0, 1e-5,
    0, LBFGS_LINESEARCH_DEFAULT, 40,
    1e-20, 1e20, 1e-4, 0.9, 0.9, 1.0e-16,
    0.0, 0, -1,
};

/* Forward function declarations. */

typedef int (*line_search_double_proc)(
    int n,
    lbfgsdoublefloatval_t *x,
    lbfgsdoublefloatval_t *f,
    lbfgsdoublefloatval_t *g,
    lbfgsdoublefloatval_t *s,
    lbfgsdoublefloatval_t *stp,
    const lbfgsdoublefloatval_t* xp,
    const lbfgsdoublefloatval_t* gp,
    lbfgsdoublefloatval_t *wa,
    callback_double_data_t *cd,
    const lbfgs_parameter_t *param
    );
    
static int line_search_double_backtracking(
    int n,
    lbfgsdoublefloatval_t *x,
    lbfgsdoublefloatval_t *f,
    lbfgsdoublefloatval_t *g,
    lbfgsdoublefloatval_t *s,
    lbfgsdoublefloatval_t *stp,
    const lbfgsdoublefloatval_t* xp,
    const lbfgsdoublefloatval_t* gp,
    lbfgsdoublefloatval_t *wa,
    callback_double_data_t *cd,
    const lbfgs_parameter_t *param
    );

static int line_search_double_backtracking_owlqn_double(
    int n,
    lbfgsdoublefloatval_t *x,
    lbfgsdoublefloatval_t *f,
    lbfgsdoublefloatval_t *g,
    lbfgsdoublefloatval_t *s,
    lbfgsdoublefloatval_t *stp,
    const lbfgsdoublefloatval_t* xp,
    const lbfgsdoublefloatval_t* gp,
    lbfgsdoublefloatval_t *wp,
    callback_double_data_t *cd,
    const lbfgs_parameter_t *param
    );

static int line_search_double_morethuente(
    int n,
    lbfgsdoublefloatval_t *x,
    lbfgsdoublefloatval_t *f,
    lbfgsdoublefloatval_t *g,
    lbfgsdoublefloatval_t *s,
    lbfgsdoublefloatval_t *stp,
    const lbfgsdoublefloatval_t* xp,
    const lbfgsdoublefloatval_t* gp,
    lbfgsdoublefloatval_t *wa,
    callback_double_data_t *cd,
    const lbfgs_parameter_t *param
    );

static int update_trial_interval_double(
    lbfgsdoublefloatval_t *x,
    lbfgsdoublefloatval_t *fx,
    lbfgsdoublefloatval_t *dx,
    lbfgsdoublefloatval_t *y,
    lbfgsdoublefloatval_t *fy,
    lbfgsdoublefloatval_t *dy,
    lbfgsdoublefloatval_t *t,
    lbfgsdoublefloatval_t *ft,
    lbfgsdoublefloatval_t *dt,
    const lbfgsdoublefloatval_t tmin,
    const lbfgsdoublefloatval_t tmax,
    int *brackt
    );

static lbfgsdoublefloatval_t owlqn_double_x1norm(
    const lbfgsdoublefloatval_t* x,
    const int start,
    const int n
    );

static void owlqn_double_pseudo_gradient(
    lbfgsdoublefloatval_t* pg,
    const lbfgsdoublefloatval_t* x,
    const lbfgsdoublefloatval_t* g,
    const int n,
    const lbfgsdoublefloatval_t c,
    const int start,
    const int end
    );

static void owlqn_double_project(
    lbfgsdoublefloatval_t* d,
    const lbfgsdoublefloatval_t* sign,
    const int start,
    const int end
    );


#if     defined(USE_SSE) && (defined(__SSE__) || defined(__SSE2__))
static int round_out_variables_double(int n)
{
    n += 7;
    n /= 8;
    n *= 8;
    return n;
}
#endif/*defined(USE_SSE)*/

lbfgsdoublefloatval_t* lbfgsdouble_malloc(int n)
{
#if     defined(USE_SSE) && (defined(__SSE__) || defined(__SSE2__))
    n = round_out_variables_double(n);
#endif/*defined(USE_SSE)*/
    return (lbfgsdoublefloatval_t*)vecdoublealloc(sizeof(lbfgsdoublefloatval_t) * n);
}

void lbfgsdouble_free(lbfgsdoublefloatval_t *x)
{
    vecdoublefree(x);
}

 #ifndef LIBFGSPARAMETERINIT
 #define LIBFGSPARAMETERINIT
void lbfgs_parameter_init(lbfgsfloat_parameter_t *param)
{
    memcpy(param, &_defparam_double, sizeof(*param));
}
#endif // LIBFGSPARAMETERINIT

int lbfgs(
    int n,
    lbfgsdoublefloatval_t *x,
    lbfgsdoublefloatval_t *ptr_fx,
    lbfgsdouble_evaluate_t proc_evaluate,
    lbfgsdouble_progress_t proc_progress,
    void *instance,
    lbfgs_parameter_t *_param
    )
{
    int ret;
    int i, j, k, ls, end, bound;
    lbfgsdoublefloatval_t step;
    /* Constant parameters and their default values. */
    lbfgs_parameter_t param = (_param != NULL) ? (*_param) : _defparam_double;
    const int m = param.m;
    lbfgsdoublefloatval_t *xp = NULL;
    lbfgsdoublefloatval_t *g = NULL, *gp = NULL, *pg = NULL;
    lbfgsdoublefloatval_t *d = NULL, *w = NULL, *pf = NULL;
    iteration_double_data_t *lm = NULL, *it = NULL;
    lbfgsdoublefloatval_t ys, yy;
    lbfgsdoublefloatval_t xnorm, gnorm, beta;
    lbfgsdoublefloatval_t fx = 0.;
    lbfgsdoublefloatval_t rate = 0.;
    line_search_double_proc linesearch = line_search_double_morethuente;
    /* Construct a callback data. */
    callback_double_data_t cd;
    cd.n = n;
    cd.instance = instance;
    cd.proc_evaluate = proc_evaluate;
    cd.proc_progress = proc_progress;
#if     defined(USE_SSE) && (defined(__SSE__) || defined(__SSE2__))
    /* Round out the number of variables. */
    n = round_out_variables_double(n);
#endif/*defined(USE_SSE)*/
    /* Check the input parameters for errors. */
    if (n <= 0) {
        return LBFGSDOUBLEERR_INVALID_N;
    }
#if     defined(USE_SSE) && (defined(__SSE__) || defined(__SSE2__))
    if (n % 8 != 0) {
        return LBFGSDOUBLEERR_INVALID_N_SSE;
    }
    if ((uintptr_t)(const void*)x % 16 != 0) {
        return LBFGSDOUBLEERR_INVALID_X_SSE;
    }
#endif/*defined(USE_SSE)*/
    if (param.epsilon < 0.) {
        return LBFGSDOUBLEERR_INVALID_EPSILON;
    }
    if (param.past < 0) {
        return LBFGSDOUBLEERR_INVALID_TESTPERIOD;
    }
    if (param.delta < 0.) {
        return LBFGSDOUBLEERR_INVALID_DELTA;
    }
    if (param.min_step < 0.) {
        return LBFGSDOUBLEERR_INVALID_MINSTEP;
    }
    if (param.max_step < param.min_step) {
        return LBFGSDOUBLEERR_INVALID_MAXSTEP;
    }
    if (param.ftol < 0.) {
        return LBFGSDOUBLEERR_INVALID_FTOL;
    }
    if (param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE ||
        param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE) {
        if (param.wolfe <= param.ftol || 1. <= param.wolfe) {
            return LBFGSDOUBLEERR_INVALID_WOLFE;
        }
    }
    if (param.gtol < 0.) {
        return LBFGSDOUBLEERR_INVALID_GTOL;
    }
    if (param.xtol < 0.) {
        return LBFGSDOUBLEERR_INVALID_XTOL;
    }
    if (param.max_linesearch <= 0) {
        return LBFGSDOUBLEERR_INVALID_MAXLINESEARCH;
    }
    if (param.orthantwise_c < 0.) {
        return LBFGSDOUBLEERR_INVALID_ORTHANTWISE;
    }
    if (param.orthantwise_start < 0 || n < param.orthantwise_start) {
        return LBFGSDOUBLEERR_INVALID_ORTHANTWISE_START;
    }
    if (param.orthantwise_end < 0) {
        param.orthantwise_end = n;
    }
    if (n < param.orthantwise_end) {
        return LBFGSDOUBLEERR_INVALID_ORTHANTWISE_END;
    }
    if (param.orthantwise_c != 0.) {
        switch (param.linesearch) {
        case LBFGS_LINESEARCH_BACKTRACKING:
            linesearch = line_search_double_backtracking_owlqn_double;
            break;
        default:
            /* Only the backtracking method is available. */
            return LBFGSDOUBLEERR_INVALID_LINESEARCH;
        }
    } else {
        switch (param.linesearch) {
        case LBFGS_LINESEARCH_MORETHUENTE:
            linesearch = line_search_double_morethuente;
            break;
        case LBFGS_LINESEARCH_BACKTRACKING_ARMIJO:
        case LBFGS_LINESEARCH_BACKTRACKING_WOLFE:
        case LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE:
            linesearch = line_search_double_backtracking;
            break;
        default:
            return LBFGSDOUBLEERR_INVALID_LINESEARCH;
        }
    }

    /* Allocate working space. */
    xp = (lbfgsdoublefloatval_t*)vecdoublealloc(n * sizeof(lbfgsdoublefloatval_t));
    g = (lbfgsdoublefloatval_t*)vecdoublealloc(n * sizeof(lbfgsdoublefloatval_t));
    gp = (lbfgsdoublefloatval_t*)vecdoublealloc(n * sizeof(lbfgsdoublefloatval_t));
    d = (lbfgsdoublefloatval_t*)vecdoublealloc(n * sizeof(lbfgsdoublefloatval_t));
    w = (lbfgsdoublefloatval_t*)vecdoublealloc(n * sizeof(lbfgsdoublefloatval_t));
    if (xp == NULL || g == NULL || gp == NULL || d == NULL || w == NULL) {
        ret = LBFGSDOUBLEERR_OUTOFMEMORY;
        goto lbfgsdouble_exit;
    }

    if (param.orthantwise_c != 0.) {
        /* Allocate working space for OW-LQN. */
        pg = (lbfgsdoublefloatval_t*)vecdoublealloc(n * sizeof(lbfgsdoublefloatval_t));
        if (pg == NULL) {
            ret = LBFGSDOUBLEERR_OUTOFMEMORY;
            goto lbfgsdouble_exit;
        }
    }

    /* Allocate limited memory storage. */
    lm = (iteration_double_data_t*)vecdoublealloc(m * sizeof(iteration_double_data_t));
    if (lm == NULL) {
        ret = LBFGSDOUBLEERR_OUTOFMEMORY;
        goto lbfgsdouble_exit;
    }

    /* Initialize the limited memory. */
    for (i = 0;i < m;++i) {
        it = &lm[i];
        it->alpha = 0;
        it->ys = 0;
        it->s = (lbfgsdoublefloatval_t*)vecdoublealloc(n * sizeof(lbfgsdoublefloatval_t));
        it->y = (lbfgsdoublefloatval_t*)vecdoublealloc(n * sizeof(lbfgsdoublefloatval_t));
        if (it->s == NULL || it->y == NULL) {
            ret = LBFGSDOUBLEERR_OUTOFMEMORY;
            goto lbfgsdouble_exit;
        }
    }

    /* Allocate an array for storing previous values of the objective function. */
    if (0 < param.past) {
        pf = (lbfgsdoublefloatval_t*)vecdoublealloc(param.past * sizeof(lbfgsdoublefloatval_t));
    }
	
    /* Evaluate the function value and its gradient. */
    fx = cd.proc_evaluate(cd.instance, x, g, cd.n, 0);
    if (0. != param.orthantwise_c) {
        /* Compute the L1 norm of the variable and add it to the object value. */
        xnorm = owlqn_double_x1norm(x, param.orthantwise_start, param.orthantwise_end);
        fx += xnorm * param.orthantwise_c;
        owlqn_double_pseudo_gradient(
            pg, x, g, n,
            param.orthantwise_c, param.orthantwise_start, param.orthantwise_end
            );
    }

    /* Store the initial value of the objective function. */
    if (pf != NULL) {
        pf[0] = fx;
    }

    /*
        Compute the direction;
        we assume the initial hessian matrix H_0 as the identity matrix.
     */
    if (param.orthantwise_c == 0.) {
        vecdoublencpy(d, g, n);
    } else {
        vecdoublencpy(d, pg, n);
    }

    /*
       Make sure that the initial variables are not a MINIMIZER_DOUBLE.
     */
    vecdouble2norm(&xnorm, x, n);
    if (param.orthantwise_c == 0.) {
        vecdouble2norm(&gnorm, g, n);
    } else {
        vecdouble2norm(&gnorm, pg, n);
    }
    if (xnorm < 1.0) xnorm = 1.0;
    if (gnorm / xnorm <= param.epsilon) {
        ret = LBFGSDOUBLE_ALREADY_MINIMIZED;
        goto lbfgsdouble_exit;
    }

    /* Compute the initial step:
        step = 1.0 / sqrt(vecdoubledot(d, d, n))
     */
    vecdouble2norminv(&step, d, n);

    k = 1;
    end = 0;
    for (;;) 
    {
        /* Store the current position and gradient vecdoubletors. */
        vecdoublecpy(xp, x, n);
        vecdoublecpy(gp, g, n);

        /* Search for an optimal step. */
        if (param.orthantwise_c == 0.) 
            ls = linesearch(n, x, &fx, g, d, &step, xp, gp, w, &cd, &param);
        else 
        {
            ls = linesearch(n, x, &fx, g, d, &step, xp, pg, w, &cd, &param);
            owlqn_double_pseudo_gradient(pg, x, g, n,param.orthantwise_c, param.orthantwise_start, param.orthantwise_end);
        }
        if (ls < 0) 
        {
            /* Revert to the previous point. */
            vecdoublecpy(x, xp, n);
            vecdoublecpy(g, gp, n);
            ret = ls;
            goto lbfgsdouble_exit;
        }

        /* Compute x and g norms. */
        vecdouble2norm(&xnorm, x, n);
        if (param.orthantwise_c == 0.) 
            vecdouble2norm(&gnorm, g, n);
        else 
            vecdouble2norm(&gnorm, pg, n);

        /* Report the progress. */
        if (cd.proc_progress)
            if ((ret = cd.proc_progress(cd.instance, x, g, fx, xnorm, gnorm, step, cd.n, k, ls)))
                goto lbfgsdouble_exit;

        /*
            Convergence test.
            The criterion is given by the following formula:
                |g(x)| / \max(1, |x|) < \epsilon
         */
        if (xnorm < 1.0) 
        		xnorm = 1.0;
        if (gnorm / xnorm <= param.epsilon) 
        {
            /* Convergence. */
            ret = LBFGSDOUBLE_SUCCESS;
            break;
        }

        /*
            Test for stopping criterion.
            The criterion is given by the following formula:
                (f(past_x) - f(x)) / f(x) < \delta
         */
        if (pf != NULL) 
        {
            /* We don't test the stopping criterion while k < past. */
            if (param.past <= k) 
            {
                /* Compute the relative improvement from the past. */
                rate = (pf[k % param.past] - fx) / fx;

                /* The stopping criterion. */
                if (rate < param.delta) 
                {
                    ret = LBFGSDOUBLE_STOP;
                    break;
                }
            }

            /* Store the current value of the objective function. */
            pf[k % param.past] = fx;
        }

        if (param.max_iterations != 0 && param.max_iterations < k+1) 
        {
            /* Maximum number of iterations. */
            ret = LBFGSDOUBLEERR_MAXIMUMITERATION;
            break;
        }

        /*
            Update vecdoubletors s and y:
                s_{k+1} = x_{k+1} - x_{k} = \step * d_{k}.
                y_{k+1} = g_{k+1} - g_{k}.
         */
        it = &lm[end];
        vecdoublediff(it->s, x, xp, n);
        vecdoublediff(it->y, g, gp, n);

        /*
            Compute scalars ys and yy:
                ys = y^t \cdot s = 1 / \rho.
                yy = y^t \cdot y.
            Notice that yy is used for scaling the hessian matrix H_0 (Cholesky factor).
         */
        vecdoubledot(&ys, it->y, it->s, n);
        vecdoubledot(&yy, it->y, it->y, n);
        it->ys = ys;

        /*
            Recursive formula to compute dir = -(H \cdot g).
                This is described in page 779 of:
                Jorge Nocedal.
                Updating Quasi-Newton Matrices with Limited Storage.
                Mathematics of Computation, Vol. 35, No. 151,
                pp. 773--782, 1980.
         */
        bound = (m <= k) ? m : k;
        ++k;
        end = (end + 1) % m;

        /* Compute the steepest direction. */
        if (param.orthantwise_c == 0.) 
            vecdoublencpy(d, g, n); /* Compute the negative of gradients. */
        else 
            vecdoublencpy(d, pg, n);

        j = end;
        for (i = 0;i < bound;++i) 
        {
            j = (j + m - 1) % m;    /* if (--j == -1) j = m-1; */
            it = &lm[j];
            /* \alpha_{j} = \rho_{j} s^{t}_{j} \cdot q_{k+1}. */
            vecdoubledot(&it->alpha, it->s, d, n);
            it->alpha /= it->ys;
            /* q_{i} = q_{i+1} - \alpha_{i} y_{i}. */
            vecdoubleadd(d, it->y, -it->alpha, n);
        }

        vecdoublescale(d, ys / yy, n);

        for (i = 0;i < bound;++i) 
        {
            it = &lm[j];
            /* \beta_{j} = \rho_{j} y^t_{j} \cdot \gamma_{i}. */
            vecdoubledot(&beta, it->y, d, n);
            beta /= it->ys;
            /* \gamma_{i+1} = \gamma_{i} + (\alpha_{j} - \beta_{j}) s_{j}. */
            vecdoubleadd(d, it->s, it->alpha - beta, n);
            j = (j + 1) % m;        /* if (++j == m) j = 0; */
        }

        /*
            Constrain the search direction for orthant-wise updates.
         */
        if (param.orthantwise_c != 0.) 
            for (i = param.orthantwise_start;i < param.orthantwise_end;++i) 
                if (d[i] * pg[i] >= 0) 
                    d[i] = 0;

        /*
            Now the search direction d is ready. We try step = 1 first.
         */
        step = 1.0;
    }

lbfgsdouble_exit:
    /* Return the final value of the objective function. */
    if (ptr_fx != NULL)
        *ptr_fx = fx;

    vecdoublefree(pf);

    /* Free memory blocks used by this function. */
    if (lm != NULL) 
    {
        for (i = 0;i < m;++i) 
        {
            vecdoublefree(lm[i].s);
            vecdoublefree(lm[i].y);
        }
        vecdoublefree(lm);
    }
    vecdoublefree(pg);
    vecdoublefree(w);
    vecdoublefree(d);
    vecdoublefree(gp);
    vecdoublefree(g);
    vecdoublefree(xp);

    return ret;
}



static int line_search_double_backtracking(
    int n,
    lbfgsdoublefloatval_t *x,
    lbfgsdoublefloatval_t *f,
    lbfgsdoublefloatval_t *g,
    lbfgsdoublefloatval_t *s,
    lbfgsdoublefloatval_t *stp,
    const lbfgsdoublefloatval_t* xp,
    const lbfgsdoublefloatval_t* gp,
    lbfgsdoublefloatval_t *wp,
    callback_double_data_t *cd,
    const lbfgs_parameter_t *param
    )
{
    int count = 0;
    lbfgsdoublefloatval_t width, dg;
    lbfgsdoublefloatval_t finit, dginit = 0., dgtest;
    const lbfgsdoublefloatval_t dec = 0.5, inc = 2.1;

    /* Check the input parameters for errors. */
    if (*stp <= 0.)
        return LBFGSDOUBLEERR_INVALIDPARAMETERS;

    /* Compute the initial gradient in the search direction. */
    vecdoubledot(&dginit, g, s, n);

    /* Make sure that s points to a descent direction. */
    if (0 < dginit)
        return LBFGSDOUBLEERR_INCREASEGRADIENT;

    /* The initial value of the objective function. */
    finit = *f;
    dgtest = param->ftol * dginit;

    for (;;) 
    {
        vecdoublecpy(x, xp, n);
        vecdoubleadd(x, s, *stp, n);

        /* Evaluate the function and gradient values. */
        *f = cd->proc_evaluate(cd->instance, x, g, cd->n, *stp);
        ++count;

        if (*f > finit + *stp * dgtest)
            width = dec;
        else 
        {
            /* The sufficient decrease condition (Armijo condition). */
            if (param->linesearch == LBFGS_LINESEARCH_BACKTRACKING_ARMIJO) 
            {
            		/* Exit with the Armijo condition. */
            		return count;
	       	  }

	          /* Check the Wolfe condition. */
	          vecdoubledot(&dg, g, s, n);
	        	if (dg < param->wolfe * dginit)
    		    		width = inc;
	        	else 
	        	{
		        		if(param->linesearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE) 
		        		{
		            		/* Exit with the regular Wolfe condition. */
		            		return count;
		        		}

		        		/* Check the strong Wolfe condition. */
		        		if(dg > -param->wolfe * dginit)
		            		width = dec;
		        		else
		        		{
		            		/* Exit with the strong Wolfe condition. */
		            		return count;
		        		}
            }
        }

        if (*stp < param->min_step) 
        {
            /* The step is the minimum value. */
            return LBFGSDOUBLEERR_MINIMUMSTEP;
        }
        if (*stp > param->max_step) 
        {
            /* The step is the maximum value. */
            return LBFGSDOUBLEERR_MAXIMUMSTEP;
        }
        if (param->max_linesearch <= count) 
        {
            /* Maximum number of iteration. */
            return LBFGSDOUBLEERR_MAXIMUMLINESEARCH;
        }

        (*stp) *= width;
    }
}



static int line_search_double_backtracking_owlqn_double(
    int n,
    lbfgsdoublefloatval_t *x,
    lbfgsdoublefloatval_t *f,
    lbfgsdoublefloatval_t *g,
    lbfgsdoublefloatval_t *s,
    lbfgsdoublefloatval_t *stp,
    const lbfgsdoublefloatval_t* xp,
    const lbfgsdoublefloatval_t* gp,
    lbfgsdoublefloatval_t *wp,
    callback_double_data_t *cd,
    const lbfgs_parameter_t *param
    )
{
    int i, count = 0;
    lbfgsdoublefloatval_t width = 0.5, norm = 0.;
    lbfgsdoublefloatval_t finit = *f, dgtest;

    /* Check the input parameters for errors. */
    if (*stp <= 0.) {
        return LBFGSDOUBLEERR_INVALIDPARAMETERS;
    }

    /* Choose the orthant for the new point. */
    for (i = 0;i < n;++i) {
        wp[i] = (xp[i] == 0.) ? -gp[i] : xp[i];
    }

    for (;;) {
        /* Update the current point. */
        vecdoublecpy(x, xp, n);
        vecdoubleadd(x, s, *stp, n);

        /* The current point is projected onto the orthant. */
        owlqn_double_project(x, wp, param->orthantwise_start, param->orthantwise_end);

        /* Evaluate the function and gradient values. */
        *f = cd->proc_evaluate(cd->instance, x, g, cd->n, *stp);

        /* Compute the L1 norm of the variables and add it to the object value. */
        norm = owlqn_double_x1norm(x, param->orthantwise_start, param->orthantwise_end);
        *f += norm * param->orthantwise_c;

        ++count;

        dgtest = 0.;
        for (i = 0;i < n;++i) {
            dgtest += (x[i] - xp[i]) * gp[i];
        }

        if (*f <= finit + param->ftol * dgtest) {
            /* The sufficient decrease condition. */
            return count;
        }

        if (*stp < param->min_step) {
            /* The step is the minimum value. */
            return LBFGSDOUBLEERR_MINIMUMSTEP;
        }
        if (*stp > param->max_step) {
            /* The step is the maximum value. */
            return LBFGSDOUBLEERR_MAXIMUMSTEP;
        }
        if (param->max_linesearch <= count) {
            /* Maximum number of iteration. */
            return LBFGSDOUBLEERR_MAXIMUMLINESEARCH;
        }

        (*stp) *= width;
    }
}



static int line_search_double_morethuente(
    int n,
    lbfgsdoublefloatval_t *x,
    lbfgsdoublefloatval_t *f,
    lbfgsdoublefloatval_t *g,
    lbfgsdoublefloatval_t *s,
    lbfgsdoublefloatval_t *stp,
    const lbfgsdoublefloatval_t* xp,
    const lbfgsdoublefloatval_t* gp,
    lbfgsdoublefloatval_t *wa,
    callback_double_data_t *cd,
    const lbfgs_parameter_t *param
    )
{
    int count = 0;
    int brackt, stage1, uinfo = 0;
    lbfgsdoublefloatval_t dg;
    lbfgsdoublefloatval_t stx, fx, dgx;
    lbfgsdoublefloatval_t sty, fy, dgy;
    lbfgsdoublefloatval_t fxm, dgxm, fym, dgym, fm, dgm;
    lbfgsdoublefloatval_t finit, ftest1, dginit, dgtest;
    lbfgsdoublefloatval_t width, prev_width;
    lbfgsdoublefloatval_t stmin, stmax;

    /* Check the input parameters for errors. */
    if (*stp <= 0.) {
        return LBFGSDOUBLEERR_INVALIDPARAMETERS;
    }

    /* Compute the initial gradient in the search direction. */
    vecdoubledot(&dginit, g, s, n);

    /* Make sure that s points to a descent direction. */
    if (0 < dginit)
        return LBFGSDOUBLEERR_INCREASEGRADIENT;

    /* Initialize local variables. */
    brackt = 0;
    stage1 = 1;
    finit = *f;
    dgtest = param->ftol * dginit;
    width = param->max_step - param->min_step;
    prev_width = 2.0 * width;

    /*
        The variables stx, fx, dgx contain the values of the step,
        function, and directional derivative at the best step.
        The variables sty, fy, dgy contain the value of the step,
        function, and derivative at the other endpoint of
        the interval of uncertainty.
        The variables stp, f, dg contain the values of the step,
        function, and derivative at the current step.
    */
    stx = sty = 0.;
    fx = fy = finit;
    dgx = dgy = dginit;

    for (;;) {
        /*
            Set the minimum and maximum steps to correspond to the
            present interval of uncertainty.
         */
        if (brackt) 
        {
            stmin = min2(stx, sty);
            stmax = max2(stx, sty);
        } 
        else 
        {
            stmin = stx;
            stmax = *stp + 4.0 * (*stp - stx);
        }

        /* Clip the step in the range of [stpmin, stpmax]. */
        if (*stp < param->min_step) *stp = param->min_step;
        if (param->max_step < *stp) *stp = param->max_step;

        /*
            If an unusual termination is to occur then let
            stp be the lowest point obtained so far.
         */
        if ((brackt && ((*stp <= stmin || stmax <= *stp) || param->max_linesearch <= count + 1 || uinfo != 0)) || (brackt && (stmax - stmin <= param->xtol * stmax))) {
            *stp = stx;
        }

        /*
            Compute the current value of x:
                x <- x + (*stp) * s.
         */
        vecdoublecpy(x, xp, n);
        vecdoubleadd(x, s, *stp, n);

        /* Evaluate the function and gradient values. */
        *f = cd->proc_evaluate(cd->instance, x, g, cd->n, *stp);
        vecdoubledot(&dg, g, s, n);

        ftest1 = finit + *stp * dgtest;
        ++count;

        /* Test for errors and convergence. */
        if (brackt && ((*stp <= stmin || stmax <= *stp) || uinfo != 0)) {
            /* Rounding errors prevent further progress. */
            return LBFGSDOUBLEERR_ROUNDING_ERROR;
        }
        if (*stp == param->max_step && *f <= ftest1 && dg <= dgtest) {
            /* The step is the maximum value. */
            return LBFGSDOUBLEERR_MAXIMUMSTEP;
        }
        if (*stp == param->min_step && (ftest1 < *f || dgtest <= dg)) {
            /* The step is the minimum value. */
            return LBFGSDOUBLEERR_MINIMUMSTEP;
        }
        if (brackt && (stmax - stmin) <= param->xtol * stmax) {
            /* Relative width of the interval of uncertainty is at most xtol. */
            return LBFGSDOUBLEERR_WIDTHTOOSMALL;
        }
        if (param->max_linesearch <= count) {
            /* Maximum number of iteration. */
            return LBFGSDOUBLEERR_MAXIMUMLINESEARCH;
        }
        if (*f <= ftest1 && fabs(dg) <= param->gtol * (-dginit)) {
            /* The sufficient decrease condition and the directional derivative condition hold. */
            return count;
        }

        /*
            In the first stage we seek a step for which the modified
            function has a nonpositive value and nonnegative derivative.
         */
        if (stage1 && *f <= ftest1 && min2(param->ftol, param->gtol) * dginit <= dg) {
            stage1 = 0;
        }

        /*
            A modified function is used to predict the step only if
            we have not obtained a step for which the modified
            function has a nonpositive function value and nonnegative
            derivative, and if a lower function value has been
            obtained but the decrease is not sufficient.
         */
        if (stage1 && ftest1 < *f && *f <= fx) {
            /* Define the modified function and derivative values. */
            fm = *f - *stp * dgtest;
            fxm = fx - stx * dgtest;
            fym = fy - sty * dgtest;
            dgm = dg - dgtest;
            dgxm = dgx - dgtest;
            dgym = dgy - dgtest;

            /*
                Call update_trial_interval_double() to update the interval of
                uncertainty and to compute the new step.
             */
            uinfo = update_trial_interval_double(
                &stx, &fxm, &dgxm,
                &sty, &fym, &dgym,
                stp, &fm, &dgm,
                stmin, stmax, &brackt
                );

            /* Reset the function and gradient values for f. */
            fx = fxm + stx * dgtest;
            fy = fym + sty * dgtest;
            dgx = dgxm + dgtest;
            dgy = dgym + dgtest;
        } else {
            /*
                Call update_trial_interval_double() to update the interval of
                uncertainty and to compute the new step.
             */
            uinfo = update_trial_interval_double(
                &stx, &fx, &dgx,
                &sty, &fy, &dgy,
                stp, f, &dg,
                stmin, stmax, &brackt
                );
        }

        /*
            Force a sufficient decrease in the interval of uncertainty.
         */
        if (brackt) {
            if (0.66 * prev_width <= fabs(sty - stx)) {
                *stp = stx + 0.5 * (sty - stx);
            }
            prev_width = width;
            width = fabs(sty - stx);
        }
    }

    return LBFGSDOUBLEERR_LOGICERROR;
}



/**
 * Define the local variables for computing MINIMIZER_DOUBLEs.
 */
#define USES_MINIMIZER_DOUBLE \
    lbfgsdoublefloatval_t a, d, gamma, theta, p, q, r, s;

/**
 * Find a MINIMIZER_DOUBLE of an interpolated cubic function.
 *  @param  cm      The MINIMIZER_DOUBLE of the interpolated cubic.
 *  @param  u       The value of one point, u.
 *  @param  fu      The value of f(u).
 *  @param  du      The value of f'(u).
 *  @param  v       The value of another point, v.
 *  @param  fv      The value of f(v).
 *  @param  du      The value of f'(v).
 */
#define CUBIC_MINIMIZER_DOUBLE(cm, u, fu, du, v, fv, dv) \
    d = (v) - (u); \
    theta = ((fu) - (fv)) * 3 / d + (du) + (dv); \
    p = fabs(theta); \
    q = fabs(du); \
    r = fabs(dv); \
    s = max3(p, q, r); \
    /* gamma = s*sqrt((theta/s)**2 - (du/s) * (dv/s)) */ \
    a = theta / s; \
    gamma = s * sqrt(a * a - ((du) / s) * ((dv) / s)); \
    if ((v) < (u)) gamma = -gamma; \
    p = gamma - (du) + theta; \
    q = gamma - (du) + gamma + (dv); \
    r = p / q; \
    (cm) = (u) + r * d;

/**
 * Find a MINIMIZER_DOUBLE of an interpolated cubic function.
 *  @param  cm      The MINIMIZER_DOUBLE of the interpolated cubic.
 *  @param  u       The value of one point, u.
 *  @param  fu      The value of f(u).
 *  @param  du      The value of f'(u).
 *  @param  v       The value of another point, v.
 *  @param  fv      The value of f(v).
 *  @param  du      The value of f'(v).
 *  @param  xmin    The maximum value.
 *  @param  xmin    The minimum value.
 */
#define CUBIC_MINIMIZER_DOUBLE2(cm, u, fu, du, v, fv, dv, xmin, xmax) \
    d = (v) - (u); \
    theta = ((fu) - (fv)) * 3 / d + (du) + (dv); \
    p = fabs(theta); \
    q = fabs(du); \
    r = fabs(dv); \
    s = max3(p, q, r); \
    /* gamma = s*sqrt((theta/s)**2 - (du/s) * (dv/s)) */ \
    a = theta / s; \
    gamma = s * sqrt(max2(0, a * a - ((du) / s) * ((dv) / s))); \
    if ((u) < (v)) gamma = -gamma; \
    p = gamma - (dv) + theta; \
    q = gamma - (dv) + gamma + (du); \
    r = p / q; \
    if (r < 0. && gamma != 0.) { \
        (cm) = (v) - r * d; \
    } else if (a < 0) { \
        (cm) = (xmax); \
    } else { \
        (cm) = (xmin); \
    }

/**
 * Find a MINIMIZER_DOUBLE of an interpolated quadratic function.
 *  @param  qm      The MINIMIZER_DOUBLE of the interpolated quadratic.
 *  @param  u       The value of one point, u.
 *  @param  fu      The value of f(u).
 *  @param  du      The value of f'(u).
 *  @param  v       The value of another point, v.
 *  @param  fv      The value of f(v).
 */
#define QUARD_MINIMIZER_DOUBLE(qm, u, fu, du, v, fv) \
    a = (v) - (u); \
    (qm) = (u) + (du) / (((fu) - (fv)) / a + (du)) / 2 * a;

/**
 * Find a MINIMIZER_DOUBLE of an interpolated quadratic function.
 *  @param  qm      The MINIMIZER_DOUBLE of the interpolated quadratic.
 *  @param  u       The value of one point, u.
 *  @param  du      The value of f'(u).
 *  @param  v       The value of another point, v.
 *  @param  dv      The value of f'(v).
 */
#define QUARD_MINIMIZER_DOUBLE2(qm, u, du, v, dv) \
    a = (u) - (v); \
    (qm) = (v) + (dv) / ((dv) - (du)) * a;

/**
 * Update a safeguarded trial value and interval for line search.
 *
 *  The parameter x represents the step with the least function value.
 *  The parameter t represents the current step. This function assumes
 *  that the derivative at the point of x in the direction of the step.
 *  If the bracket is set to true, the MINIMIZER_DOUBLE has been bracketed in
 *  an interval of uncertainty with endpoints between x and y.
 *
 *  @param  x       The pointer to the value of one endpoint.
 *  @param  fx      The pointer to the value of f(x).
 *  @param  dx      The pointer to the value of f'(x).
 *  @param  y       The pointer to the value of another endpoint.
 *  @param  fy      The pointer to the value of f(y).
 *  @param  dy      The pointer to the value of f'(y).
 *  @param  t       The pointer to the value of the trial value, t.
 *  @param  ft      The pointer to the value of f(t).
 *  @param  dt      The pointer to the value of f'(t).
 *  @param  tmin    The minimum value for the trial value, t.
 *  @param  tmax    The maximum value for the trial value, t.
 *  @param  brackt  The pointer to the predicate if the trial value is
 *                  bracketed.
 *  @retval int     Status value. Zero indicates a normal termination.
 *  
 *  @see
 *      Jorge J. More and David J. Thuente. Line search algorithm with
 *      guaranteed sufficient decrease. ACM Transactions on Mathematical
 *      Software (TOMS), Vol 20, No 3, pp. 286-307, 1994.
 */
static int update_trial_interval_double(
    lbfgsdoublefloatval_t *x,
    lbfgsdoublefloatval_t *fx,
    lbfgsdoublefloatval_t *dx,
    lbfgsdoublefloatval_t *y,
    lbfgsdoublefloatval_t *fy,
    lbfgsdoublefloatval_t *dy,
    lbfgsdoublefloatval_t *t,
    lbfgsdoublefloatval_t *ft,
    lbfgsdoublefloatval_t *dt,
    const lbfgsdoublefloatval_t tmin,
    const lbfgsdoublefloatval_t tmax,
    int *brackt
    )
{
    int bound;
    int dsign = fsigndiff_double(dt, dx);
    lbfgsdoublefloatval_t mc; /* MINIMIZER_DOUBLE of an interpolated cubic. */
    lbfgsdoublefloatval_t mq; /* MINIMIZER_DOUBLE of an interpolated quadratic. */
    lbfgsdoublefloatval_t newt;   /* new trial value. */
    USES_MINIMIZER_DOUBLE;     /* for CUBIC_MINIMIZER_DOUBLE and QUARD_MINIMIZER_DOUBLE. */

    /* Check the input parameters for errors. */
    if (*brackt) {
        if (*t <= min2(*x, *y) || max2(*x, *y) <= *t) {
            /* The trival value t is out of the interval. */
            return LBFGSDOUBLEERR_OUTOFINTERVAL;
        }
        if (0. <= *dx * (*t - *x)) {
            /* The function must decrease from x. */
            return LBFGSDOUBLEERR_INCREASEGRADIENT;
        }
        if (tmax < tmin) {
            /* Incorrect tmin and tmax specified. */
            return LBFGSDOUBLEERR_INCORRECT_TMINMAX;
        }
    }

    /*
        Trial value selection.
     */
    if (*fx < *ft) {
        /*
            Case 1: a higher function value.
            The minimum is brackt. If the cubic MINIMIZER_DOUBLE is closer
            to x than the quadratic one, the cubic one is taken, else
            the average of the MINIMIZER_DOUBLEs is taken.
         */
        *brackt = 1;
        bound = 1;
        CUBIC_MINIMIZER_DOUBLE(mc, *x, *fx, *dx, *t, *ft, *dt);
        QUARD_MINIMIZER_DOUBLE(mq, *x, *fx, *dx, *t, *ft);
        if (fabs(mc - *x) < fabs(mq - *x)) {
            newt = mc;
        } else {
            newt = mc + 0.5 * (mq - mc);
        }
    } else if (dsign) {
        /*
            Case 2: a lower function value and derivatives of
            opposite sign. The minimum is brackt. If the cubic
            MINIMIZER_DOUBLE is closer to x than the quadratic (secant) one,
            the cubic one is taken, else the quadratic one is taken.
         */
        *brackt = 1;
        bound = 0;
        CUBIC_MINIMIZER_DOUBLE(mc, *x, *fx, *dx, *t, *ft, *dt);
        QUARD_MINIMIZER_DOUBLE2(mq, *x, *dx, *t, *dt);
        if (fabs(mc - *t) > fabs(mq - *t)) {
            newt = mc;
        } else {
            newt = mq;
        }
    } else if (fabs(*dt) < fabs(*dx)) {
        /*
            Case 3: a lower function value, derivatives of the
            same sign, and the magnitude of the derivative decreases.
            The cubic MINIMIZER_DOUBLE is only used if the cubic tends to
            infinity in the direction of the MINIMIZER_DOUBLE or if the minimum
            of the cubic is beyond t. Otherwise the cubic MINIMIZER_DOUBLE is
            defined to be either tmin or tmax. The quadratic (secant)
            MINIMIZER_DOUBLE is also computed and if the minimum is brackt
            then the the MINIMIZER_DOUBLE closest to x is taken, else the one
            farthest away is taken.
         */
        bound = 1;
        CUBIC_MINIMIZER_DOUBLE2(mc, *x, *fx, *dx, *t, *ft, *dt, tmin, tmax);
        QUARD_MINIMIZER_DOUBLE2(mq, *x, *dx, *t, *dt);
        if (*brackt) {
            if (fabs(*t - mc) < fabs(*t - mq)) {
                newt = mc;
            } else {
                newt = mq;
            }
        } else {
            if (fabs(*t - mc) > fabs(*t - mq)) {
                newt = mc;
            } else {
                newt = mq;
            }
        }
    } else {
        /*
            Case 4: a lower function value, derivatives of the
            same sign, and the magnitude of the derivative does
            not decrease. If the minimum is not brackt, the step
            is either tmin or tmax, else the cubic MINIMIZER_DOUBLE is taken.
         */
        bound = 0;
        if (*brackt) {
            CUBIC_MINIMIZER_DOUBLE(newt, *t, *ft, *dt, *y, *fy, *dy);
        } else if (*x < *t) {
            newt = tmax;
        } else {
            newt = tmin;
        }
    }

    /*
        Update the interval of uncertainty. This update does not
        depend on the new step or the case analysis above.

        - Case a: if f(x) < f(t),
            x <- x, y <- t.
        - Case b: if f(t) <= f(x) && f'(t)*f'(x) > 0,
            x <- t, y <- y.
        - Case c: if f(t) <= f(x) && f'(t)*f'(x) < 0, 
            x <- t, y <- x.
     */
    if (*fx < *ft) {
        /* Case a */
        *y = *t;
        *fy = *ft;
        *dy = *dt;
    } else {
        /* Case c */
        if (dsign) {
            *y = *x;
            *fy = *fx;
            *dy = *dx;
        }
        /* Cases b and c */
        *x = *t;
        *fx = *ft;
        *dx = *dt;
    }

    /* Clip the new trial value in [tmin, tmax]. */
    if (tmax < newt) newt = tmax;
    if (newt < tmin) newt = tmin;

    /*
        Redefine the new trial value if it is close to the upper bound
        of the interval.
     */
    if (*brackt && bound) {
        mq = *x + 0.66 * (*y - *x);
        if (*x < *y) {
            if (mq < newt) newt = mq;
        } else {
            if (newt < mq) newt = mq;
        }
    }

    /* Return the new trial value. */
    *t = newt;
    return 0;
}





static lbfgsdoublefloatval_t owlqn_double_x1norm(
    const lbfgsdoublefloatval_t* x,
    const int start,
    const int n
    )
{
    int i;
    lbfgsdoublefloatval_t norm = 0.;

    for (i = start;i < n;++i) {
        norm += fabs(x[i]);
    }

    return norm;
}

static void owlqn_double_pseudo_gradient(
    lbfgsdoublefloatval_t* pg,
    const lbfgsdoublefloatval_t* x,
    const lbfgsdoublefloatval_t* g,
    const int n,
    const lbfgsdoublefloatval_t c,
    const int start,
    const int end
    )
{
    int i;

    /* Compute the negative of gradients. */
    for (i = 0;i < start;++i) {
        pg[i] = g[i];
    }

    /* Compute the psuedo-gradients. */
    for (i = start;i < end;++i) {
        if (x[i] < 0.) {
            /* Differentiable. */
            pg[i] = g[i] - c;
        } else if (0. < x[i]) {
            /* Differentiable. */
            pg[i] = g[i] + c;
        } else {
            if (g[i] < -c) {
                /* Take the right partial derivative. */
                pg[i] = g[i] + c;
            } else if (c < g[i]) {
                /* Take the left partial derivative. */
                pg[i] = g[i] - c;
            } else {
                pg[i] = 0.;
            }
        }
    }

    for (i = end;i < n;++i) {
        pg[i] = g[i];
    }
}

static void owlqn_double_project(
    lbfgsdoublefloatval_t* d,
    const lbfgsdoublefloatval_t* sign,
    const int start,
    const int end
    )
{
    int i;

    for (i = start;i < end;++i) {
        if (d[i] * sign[i] <= 0) {
            d[i] = 0;
        }
    }
}

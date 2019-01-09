// Copyright (C) 2016 Yixuan Qiu <yixuan.qiu@cos.name>
// Under MIT license

#ifndef LINE_SEARCHMT_H
#define LINE_SEARCHMT_H

#define min2(a, b)      ((a) <= (b) ? (a) : (b))
#define max2(a, b)      ((a) >= (b) ? (a) : (b))
/**
 * Define the local variables for computing minimizers.
 */
#define USES_MINIMIZER \
    Scalar a, d, gamma, theta, p, q, r, s;

/**
 * Find a minimizer of an interpolated cubic function.
 *  cm      The minimizer of the interpolated cubic.
 *  u       The value of one point, u.
 *  fu      The value of f(u).
 *  du      The value of f'(u).
 *  v       The value of another point, v.
 *  fv      The value of f(v).
 *  du      The value of f'(v).
 */
#define CUBIC_MINIMIZER(cm, u, fu, du, v, fv, dv) \
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
 * Find a minimizer of an interpolated cubic function.
 *  cm      The minimizer of the interpolated cubic.
 *  u       The value of one point, u.
 *  fu      The value of f(u).
 *  du      The value of f'(u).
 *  v       The value of another point, v.
 *  fv      The value of f(v).
 *  du      The value of f'(v).
 *  xmin    The maximum value.
 *  xmin    The minimum value.
 */
#define CUBIC_MINIMIZER2(cm, u, fu, du, v, fv, dv, xmin, xmax) \
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
 * Find a minimizer of an interpolated quadratic function.
 *  qm      The minimizer of the interpolated quadratic.
 *  u       The value of one point, u.
 *  fu      The value of f(u).
 *  du      The value of f'(u).
 *  v       The value of another point, v.
 *  fv      The value of f(v).
 */
#define QUADR_MINIMIZER(qm, u, fu, du, v, fv) \
    a = (v) - (u); \
    (qm) = (u) + (du) / (((fu) - (fv)) / a + (du)) / 2 * a;

/**
 * Find a minimizer of an interpolated quadratic function.
 *  qm      The minimizer of the interpolated quadratic.
 *  u       The value of one point, u.
 *  du      The value of f'(u).
 *  v       The value of another point, v.
 *  dv      The value of f'(v).
 */
#define QUADR_MINIMIZER2(qm, u, du, v, dv) \
    a = (u) - (v); \
    (qm) = (v) + (dv) / ((dv) - (du)) * a;

#include <Eigen/Core>
#include <stdexcept>  // std::runtime_error 

namespace LBFGSpp {


    //
    // Line search algorithms for LBFGS. Mainly for internal use.
    //
    template <typename Scalar>
        class LineSearchMT
        {
            private:
                typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

            public:
                //
                // Line search by backtracking.
                //
                // \param f      A function object such that `f(x, grad)` returns the
                //               objective function value at `x`, and overwrites `grad` with
                //               the gradient.
                // \param fout   In: The objective function value at the current point.
                //               Out: The function value at the new point.
                // \param x      Out: The new point moved to.
                // \param grad   In: The current gradient vector. 
                //               Out: The gradient at the new point.
                // \param step   In: The initial step length.
                //               Out: The calculated step length.
                // \param drt    The current moving direction.
                // \param xp     The current point.
                // \param param  Parameters for the LBFGS algorithm
                //
                template <typename Foo>
                    static void MoreThuente(
                            Foo& f,
                            Scalar& fout,
                            Vector& x,
                            Vector& grad,
                            Scalar& step,
                            const Vector& drt,
                            const Vector& xp,
                            const LBFGSParam<Scalar>& param)
                    {
                        // Check the value of step
                        if(step <= Scalar(0))
                            std::invalid_argument("'step' must be positive");

                        // Save the function value at the current x
                        const Scalar f_init = fout;
                        // Projection of gradient on the search direction
                        const Scalar dg_init = grad.dot(drt);
                        // Make sure d points to a descent direction
                        if(dg_init > 0)
                            std::logic_error("the moving direction increases the objective function value");
                        //
                        //

                        const Scalar dg_test = param.wolfec1 * dg_init;
                        /*
                           The variables stx, fx, dgx contain the values of the step,
                           function, and directional derivative at the best step.
                           The variables sty, fy, dgy contain the value of the step,
                           function, and derivative at the other endpoint of
                           the interval of uncertainty.
                           The variables stp, f, dg contain the values of the step,
                           function, and derivative at the current step.
                           */

                        Scalar stx {0}, sty {0}, stmin, stmax;
                        Scalar fx, fy, dgx, dgy;
                        Scalar fxm, dgxm, fym, dgym, fm, dgm;
                        Scalar ftest1;
                        Scalar width = param.max_step - param.min_step;
                        Scalar prev_width = 2.0 * width;
                        int bracketed {0}, stage1 {1}, uinfo {0};
                        //
                        dgx = dgy = dg_init;
                        fx = fy = f_init;
                        // finit is f_init;

                        for(int iter = 0; iter < param.max_linesearch; iter++)
                        {

                            if (bracketed) {
                                stmin = min2(stx, sty);
                                stmax = max2(stx, sty);
                            } else {
                                stmin = stx;
                                stmax = step + 4.0 * (step - stx);
                            }
                            /* Clip the step in the range of [stpmin, stpmax]. */
                            if (step < param.min_step) step = param.min_step;
                            if (param.max_step < step) step = param.max_step;
                            /*
                               If an unusual termination is to occur then let
                               step be the lowest point obtained so far.
                               */
                            if ((bracketed && ((step <= stmin || step >= stmax) || uinfo != 0 )) || (bracketed && (stmax - stmin <= param.delta * stmax))) {
                                step = stx;
                            }
                            /*
                               Compute the current value of x:
                               x <- x + (*stp) * s.
                               */
                            x.noalias() = xp + step * drt;

                            /* Evaluate the function and gradient values. */
                            fout = f(x, grad);
                            const Scalar dg = grad.dot(drt);

                            ftest1 = f_init + step * dg_test;

                            if (bracketed && ((step <= stmin || stmax <= step) ))
                                throw std::runtime_error("Rounding errors prevent further progress.");
                            if (step == param.max_step && fout <= ftest1 && dg <= dg_test)
                                throw std::runtime_error("The step is the maximum value.");
                            if (step == param.min_step && (ftest1 < fout || dg_test <= dg))
                                throw std::runtime_error("The step is the minimum value.");
                            if (bracketed && (stmax - stmin) <= param.delta * stmax)
                                throw std::runtime_error("Relative width of the interval of uncertainty is at most param.delta.");
                            if (fout <= ftest1 && abs(dg) <= param.wolfec2 * (-dg_init)) 
                                /* The strong Wolfe conditions hold. */
                                return;
                            /*
                               In the first stage we seek a step for which the modified
                               function has a nonpositive value and nonnegative derivative.
                               */
                            if (stage1 && fout <= ftest1 && dg >= dg_test) {
                                stage1 = 0;
                            }

                            /*
                               A modified function is used to predict the step only if
                               we have not obtained a step for which the modified
                               function has a nonpositive function value and nonnegative
                               derivative, and if a lower function value has been
                               obtained but the decrease is not sufficient.
                               */
                            if (stage1 && fout > ftest1 && fout <= fx) {
                                /* Define the modified function and derivative values. */
                                fm = fout - step * dg_test; 
                                fxm = fx - stx * dg_test;
                                fym = fy - sty * dg_test;
                                dgm = dg - dg_test;
                                dgxm = dgx - dg_test;
                                dgym = dgy - dg_test;

                                /*
                                   Call update_trial_interval() to update the interval of
                                   uncertainty and to compute the new step.
                                   */
                                uinfo = update_trial_interval(
                                        stx, fxm, dgxm,
                                        sty, fym, dgym,
                                        step, fm, dgm,
                                        stmin, stmax, bracketed
                                        );

                                /* Reset the function and gradient values for f. */
                                fx = fxm + stx * dg_test;
                                fy = fym + sty * dg_test;
                                dgx = dgxm + dg_test;
                                dgy = dgym + dg_test;
                            } else {
                                /*
                                   Call update_trial_interval() to update the interval of
                                   uncertainty and to compute the new step.
                                   */
                                uinfo = update_trial_interval(
                                        stx, fx, dgx,
                                        sty, fy, dgy,
                                        step, f, dg,
                                        stmin, stmax, bracketed
                                        );
                            }

                            /*
                               Force a sufficient decrease in the interval of uncertainty.
                               */
                            if (bracketed) {
                                // TODO: replace fabs
                                if (0.66 * prev_width <= fabs(sty - stx)) {
                                    step = stx + 0.5 * (sty - stx);
                                }
                                prev_width = width;
                                width = fabs(sty - stx);
                            }



                            // // TODO: remove this
                            // // x_{k+1} = x_k + step * d_k
                            // // Evaluate this candidate
                            // fx = f(x, grad);

                            // if(fx > fx_init + step * dg_test)
                            // {
                            //     width = dec;
                            // } else {
                            //     // Armijo condition is met
                            //     if(param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)
                            //         break;

                            //     const Scalar dg = grad.dot(drt);
                            //     if(dg < param.wolfec2 * dg_init)
                            //     {
                            //         width = inc;
                            //     } else {
                            //         // Regular Wolfe condition is met
                            //         if(param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE)
                            //             break;

                            //         if(dg > -param.wolfec2 * dg_init)
                            //         {
                            //             width = dec;
                            //         } else {
                            //             // Strong Wolfe condition is met
                            //             break;
                            //         }
                            //     }
                            // }

                            // if(step < param.min_step)
                            //     throw std::runtime_error("the line search step became smaller than the minimum value allowed");

                            // if(step > param.max_step)
                            //     throw std::runtime_error("the line search step became larger than the maximum value allowed");

                            // step *= width;
                        }
                    } 

                /**
                 * update_trial_interval
                 * Update a safeguarded trial value and interval for line search.
                 *
                 *  The parameter x represents the step with the least function value.
                 *  The parameter t represents the current step. This function assumes
                 *  that the derivative at the point of x in the direction of the step.
                 *  If the bracket is set to true, the minimizer has been bracketed in
                 *  an interval of uncertainty with endpoints between x and y.
                 *
                 *  x       The value of one endpoint.
                 *  fx      The value of f(x).
                 *  dx      The value of f'(x).
                 *  y       The value of another endpoint.
                 *  fy      The value of f(y).
                 *  dy      The value of f'(y).
                 *  t       The value of the trial value, t.
                 *  ft      The value of f(t).
                 *  dt      The value of f'(t).
                 *  tmin    The minimum value for the trial value, t.
                 *  tmax    The maximum value for the trial value, t.
                 *  bracketed  The predicate if the trial value is bracketed.
                 *  @retval int     Status value. Zero indicates a normal termination.
                 *  
                 *  @see
                 *      Jorge J. More and David J. Thuente. Line search algorithm with
                 *      guaranteed sufficient decrease. ACM Transactions on Mathematical
                 *      Software (TOMS), Vol 20, No 3, pp. 286-307, 1994.
                 */
                static int update_trial_interval(
                        Scalar &x,
                        Scalar &fx,
                        Scalar &dx,
                        Scalar &y,
                        Scalar &fy,
                        Scalar &dy,
                        Scalar &t,
                        Scalar &ft,
                        Scalar &dt,
                        const Scalar tmin,
                        const Scalar tmax,
                        int &bracketed
                        )
                {
                    int bound;
                    int dsign = fsigndiff(dt, dx);
                    Scalar mc; /* minimizer of an interpolated cubic. */
                    Scalar mq; /* minimizer of an interpolated quadratic. */
                    Scalar newt;   /* new trial value. */
                    USES_MINIMIZER;     /* for CUBIC_MINIMIZER and QUADR_MINIMIZER. */

                    /* Check the input parameters for errors. */
                    if (bracketed) {
                        if (t <= min2(x, y) || max2(x, y) <= t) {
                            /* The trial value t is out of the interval. */
                            return 1;
                        }
                        if (0. <= dx * (t - x)) {
                            /* The function must decrease from x. */
                            return 2;
                        }
                        if (tmax < tmin) {
                            /* Incorrect tmin and tmax specified. */
                            return 3;
                        }
                    }

                    /*
                       Trial value selection.
                       */
                    if (fx < ft) {
                        /*
                           Case 1: a higher function value.
                           The minimum is bracketed. If the cubic minimizer is closer
                           to x than the quadratic one, the cubic one is taken, else
                           the average of the minimizers is taken.
                           */
                        bracketed = 1;
                        bound = 1;
                        CUBIC_MINIMIZER(mc, x, fx, dx, t, ft, dt);
                        QUADR_MINIMIZER(mq, x, fx, dx, t, ft);
                        // TODO: replace fabs
                        if (fabs(mc - x) < fabs(mq - x)) {
                            newt = mc;
                        } else {
                            newt = mc + 0.5 * (mq - mc);
                        }
                    } else if (dsign) {
                        /*
                           Case 2: a lower function value and derivatives of
                           opposite sign. The minimum is bracketed. If the cubic
                           minimizer is closer to x than the quadratic (secant) one,
                           the cubic one is taken, else the quadratic one is taken.
                           */
                        bracketed = 1;
                        bound = 0;
                        CUBIC_MINIMIZER(mc, x, fx, dx, t, ft, dt);
                        QUADR_MINIMIZER2(mq, x, dx, t, dt);
                        if (fabs(mc - t) > fabs(mq - t)) {
                            newt = mc;
                        } else {
                            newt = mq;
                        }
                    } else if (fabs(dt) < fabs(dx)) {
                        /*
                           Case 3: a lower function value, derivatives of the
                           same sign, and the magnitude of the derivative decreases.
                           The cubic minimizer is only used if the cubic tends to
                           infinity in the direction of the minimizer or if the minimum
                           of the cubic is beyond t. Otherwise the cubic minimizer is
                           defined to be either tmin or tmax. The quadratic (secant)
                           minimizer is also computed and if the minimum is bracketed
                           then the the minimizer closest to x is taken, else the one
                           farthest away is taken.
                           */
                        bound = 1;
                        CUBIC_MINIMIZER2(mc, x, fx, dx, t, ft, dt, tmin, tmax);
                        QUADR_MINIMIZER2(mq, x, dx, t, dt);
                        if (bracketed) {
                            if (fabs(t - mc) < fabs(t - mq)) {
                                newt = mc;
                            } else {
                                newt = mq;
                            }
                        } else {
                            if (fabs(t - mc) > fabs(t - mq)) {
                                newt = mc;
                            } else {
                                newt = mq;
                            }
                        }
                    } else {
                        /*
                           Case 4: a lower function value, derivatives of the
                           same sign, and the magnitude of the derivative does
                           not decrease. If the minimum is not bracketed, the step
                           is either tmin or tmax, else the cubic minimizer is taken.
                           */
                        bound = 0;
                        if (bracketed) {
                            CUBIC_MINIMIZER(newt, t, ft, dt, y, fy, dy);
                        } else if (x < t) {
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
                       - Case b: if f(t) <= f(x) && f'(t)f'(x) > 0,
                       x <- t, y <- y.
                       - Case c: if f(t) <= f(x) && f'(t)f'(x) < 0, 
                       x <- t, y <- x.
                       */
                    if (fx < ft) {
                        /* Case a */
                        y = t;
                        fy = ft;
                        dy = dt;
                    } else {
                        /* Case c */
                        if (dsign) {
                            y = x;
                            fy = fx;
                            dy = dx;
                        }
                        /* Cases b and c */
                        x = t;
                        fx = ft;
                        dx = dt;
                    }

                    /* Clip the new trial value in [tmin, tmax]. */
                    if (tmax < newt) newt = tmax;
                    if (newt < tmin) newt = tmin;

                    /*
                       Redefine the new trial value if it is close to the upper bound
                       of the interval.
                       */
                    if (bracketed && bound) {
                        mq = x + 0.66 * (y - x);
                        if (x < y) {
                            if (mq < newt) newt = mq;
                        } else {
                            if (newt < mq) newt = mq;
                        }
                    }

                    /* Return the new trial value. */
                    t = newt;
                    return 0;
                }

        };


} // namespace LBFGSpp

#endif // LINE_SEARCHMT_H

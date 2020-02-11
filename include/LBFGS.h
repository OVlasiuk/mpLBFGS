// Copyright (c) 2016 Yixuan Qiu
// Copyright (c) 2019 Alex Vlasiuk <oleksandr.vlasiuk@gmail.com>

#ifndef LBFGS_H
#define LBFGS_H

#include <eigen3/Eigen/Core>
#include "LBFGS/LineSearchMT.h"


namespace LBFGSpp {


    //
    // LBFGS solver for unconstrained numerical optimization
    //
    template <typename Scalar>
        class LBFGSSolver
        {
            private:
                typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
                typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
                typedef Eigen::Map<Vector> MapVec;

                const LBFGSParam<Scalar>& m_param;  // Parameters to control the LBFGS algorithm
                Matrix                    m_s;      // History of the s vectors
                Matrix                    m_y;      // History of the y vectors
                Vector                    m_ys;     // History of the s'y values
                Vector                    m_alpha;  // History of the step lengths
                Vector                    m_fx;     // History of the objective function values
                Vector                    x_old;     // Old x
                Vector                    grad_new;   // New gradient
                Vector                    grad_old;  // Old gradient
                Vector                    drt_mov;    // Moving direction

                inline void reset(int n)
                {
                    const int m = m_param.m;
                    m_s.resize(n, m);
                    m_y.resize(n, m);
                    m_ys.resize(m);
                    m_alpha.resize(m);
                    x_old.resize(n);
                    grad_new.resize(n);
                    grad_old.resize(n);
                    drt_mov.resize(n);
                    if(m_param.past > 0)
                        m_fx.resize(m_param.past);
                }

            public:
                //
                // Constructor for LBFGS solver.
                //
                // param An object of LBFGSParam to store parameters for the algorithm
                //
                LBFGSSolver(const LBFGSParam<Scalar>& param) :
                    m_param(param)
            {
                m_param.check_param();
            }

                //
                // Minimizing a multivariate function using LBFGS algorithm.
                // Exceptions will be thrown if an error occurs.
                //
                // f  A function object such that `f(x, grad)` returns the
                //    objective function value at `x`, and overwrites `grad` with
                //    the gradient.
                // x  In: An initial guess of the optimal point. Out: The best point
                //    found.
                // fx Out: The objective function value at `x`.
                //
                // @retval Number of iterations used.
                //
                template <typename Foo>
                    int minimize(Foo f, Vector& x, Scalar& fx, Scalar& gnorm_return)
                    {
                        const int n = x.size();
                        const int past_param = m_param.past;
                        reset(n);

                        // Evaluate function and compute gradient
                        fx = f(x, grad_new);
                        Scalar xnorm = x.norm();
                        Scalar gnorm = grad_new.norm();
                        if(past_param > 0)
                            m_fx[0] = fx;

                        // Early exit if the initial x is already a minimizer
                        if(gnorm <= m_param.epsilon * std::max(xnorm, Scalar(1.0)))
                        {
                            gnorm_return = gnorm;
                            return 1;
                        }

                        // Initial direction
                        drt_mov.noalias() = -grad_new;
                        // Initial step
                        Scalar step = Scalar(1.0); // /drt_mov.norm();

                        int k = 1;
                        int end = 0;
                        for( ; ; )
                        {
                            // Save the curent x and gradient
                            x_old.noalias() = x;
                            grad_old.noalias() = grad_new;

                            // Line search to update x, fx and gradient
                            //
                            LineSearchMT<Scalar>::MoreThuente(f, fx, x, grad_new, step, drt_mov, x_old, m_param);

                            // New x norm and gradient norm
                            xnorm = x.norm();
                            gnorm = grad_new.norm();

                            // Convergence test -- gradient
                            if(gnorm <= m_param.epsilon * std::max(xnorm, Scalar(1.0)))
                            {
                                gnorm_return = gnorm;
                                return k;
                            }
                            // Convergence test -- objective function value
                            if(past_param > 0)
                            {
                                if(k >= past_param && mpfr::abs((m_fx[k % past_param] - fx) / fx) < m_param.delta)
                                    return k;

                                m_fx[k % past_param] = fx;
                            }
                            // Maximum number of iterations
                            if(m_param.max_iterations != 0 && k >= m_param.max_iterations)
                            {
                                return k;
                            }

                            // Update s and y
                            // s_{k+1} = x_{k+1} - x_k
                            // y_{k+1} = g_{k+1} - g_k
                            MapVec svec(&m_s(0, end), n);
                            MapVec yvec(&m_y(0, end), n);
                            svec.noalias() = x - x_old;
                            yvec.noalias() = grad_new - grad_old;

                            // ys = y's = 1/rho
                            // yy = y'y
                            Scalar ys = yvec.dot(svec);
                            Scalar yy = yvec.squaredNorm();
                            m_ys[end] = ys;

                            // Recursive formula to compute d = -H * g
                            drt_mov.noalias() = -grad_new;
                            int bound = std::min(m_param.m, k);
                            end = (end + 1) % m_param.m;
                            int j = end;
                            for(int i = 0; i < bound; i++)
                            {
                                j = (j + m_param.m - 1) % m_param.m;
                                MapVec sj(&m_s(0, j), n);
                                MapVec yj(&m_y(0, j), n);
                                m_alpha[j] = sj.dot(drt_mov) / m_ys[j];
                                drt_mov.noalias() -= m_alpha[j] * yj;
                            }

                            drt_mov *= (ys / yy);

                            for(int i = 0; i < bound; i++)
                            {
                                MapVec sj(&m_s(0, j), n);
                                MapVec yj(&m_y(0, j), n);
                                Scalar beta = yj.dot(drt_mov) / m_ys[j];
                                drt_mov.noalias() += (m_alpha[j] - beta) * sj;
                                j = (j + 1) % m_param.m;
                            }

                            // initial guess for step
                            step = Scalar(1.0); //  /drt_mov.norm();
                            k++;
                        }

                        return k;
                    }
        };


} // namespace LBFGSpp

#endif // LBFGS_H

// Copyright (c) 2019 Alex Vlasiuk <oleksandr.vlasiuk@gmail.com>

#ifndef DESCENT_H
#define DESCENT_H

#include <eigen3/Eigen/Core>
#include "../Param.h"
#include "../linesearch/LineSearchMT.h"
#include "../linesearch/LineSearch.h"


namespace LBFGSpp {


    //
    // LBFGS solver for unconstrained numerical optimization
    //
    template <typename Scalar>
        class GradSolver
        {
            private:
                typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
                typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
                typedef Eigen::Map<Vector> MapVec;

                const LBFGSParam<Scalar>& m_param;  // Parameters to control the gradient descent algorithm
                Vector                    x_old;     // Old x
                Vector                    grad_new;   // New gradient
                Vector                    grad_old;  // Old gradient
                Vector                    drt_mov;    // Moving direction

                inline void reset(int n)
                {
                    x_old.resize(n);
                    grad_new.resize(n);
                    grad_old.resize(n);
                    drt_mov.resize(n);
                }

            public:
                GradSolver(const LBFGSParam<Scalar>& param) :
                    m_param(param)
            {
                m_param.check_param();
            } 
                template <typename Foo>
                    int minimize(Foo f, Vector& x, Scalar& fout, Scalar& gnorm_return);
        };

    template <typename Scalar> 
        template <typename Foo>
        int GradSolver<Scalar>::minimize(Foo f, Vector& x, Scalar& fout, Scalar& gnorm_return)
        //
        // Minimizing a multivariate function using gradient descent.
        // Exceptions will be thrown if an error occurs.
        //
        // f  A function object such that `f(x, grad)` returns the
        //    objective function value at `x`, and overwrites `grad` with
        //    the gradient.
        // x  In: An initial guess of the optimal point. Out: The best point
        //    found.
        // fout: The objective function value at `x`.
        //
        // return: Number of iterations used.
        //
        { 
            const int N = x.size(); 
            reset(N);
            Scalar xnorm;
            Scalar gnorm;
            fout = f(x, grad_new);
            for(int i=1; i<=m_param.max_iterations; i++)
            {
                // Search direction
                drt_mov.noalias() = -grad_new;
                // Initial step
                Scalar step = fout / grad_new.norm();
                x_old.noalias() = x;
                LineSearch<Scalar>::Backtracking(f, fout, x, grad_new, step, drt_mov, x_old, m_param);
                gnorm = grad_new.norm();
                xnorm = x.norm(); 
                if(gnorm <= m_param.epsilon * std::max(xnorm, Scalar(1.0)))
                {
                    gnorm_return = gnorm;
                    return i;
                }
            }
            return -1;

        }

} // namespace LBFGSpp

#endif // DESCENT_H

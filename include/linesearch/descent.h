// Copyright (c) 2019 Alex Vlasiuk <oleksandr.vlasiuk@gmail.com>

#ifndef DESCENT_H
#define DESCENT_H

#include <eigen3/Eigen/Core>
#include "LBFGS/Param.h"
#include "LBFGS/LineSearch.h"
#include "LBFGS/LineSearchMT.h"


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
                //
                // Constructor for LBFGS solver.
                //
                // param An object of LBFGSParam to store parameters for the algorithm
                //
                GradSolver(const LBFGSParam<Scalar>& param) :
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
                // fout: The objective function value at `x`.
                //
                // @retval Number of iterations used.
                //
                template <typename Foo>
                    int minimize(Foo f, Vector& x, Scalar& fout, Scalar& gnorm_return, int dim)
                    { 
                        const int N = x.size();

                        const int n = N / (dim+1);
                        reset(N);
                        Scalar xnorm;
                        Scalar gnorm;
                        std::cout.precision(10);
                        std::cout << "N = " << n*dim << "\t n = " << n << std::endl;
                        for(int i=1; i<=m_param.max_iterations; i++)
                        {
                            fout = f(x, grad_new);
                            // Search direction
                            drt_mov.noalias() = -grad_new;
                            // Initial step
                            Scalar step = Scalar(1);
                            x_old.noalias() = x;
                            LineSearchMT<Scalar>::MoreThuente(f, fout, x, grad_new, step, drt_mov, x_old, m_param);
                            gnorm = grad_new.norm();
                            std::cout << step << '\t' << gnorm  << std::endl;
                            if(gnorm <= m_param.epsilon * std::max(xnorm, Scalar(1.0)))
                            {
                                gnorm_return = gnorm;
                                return i;
                            }
                            xnorm = x.norm();

                            mpreal temp = x.segment(n*dim, n).norm(); 
                            x.segment(n*dim, n) /= temp;
                            
                        }
                        return -1;

                    }
        };


} // namespace LBFGSpp

#endif // DESCENT_H

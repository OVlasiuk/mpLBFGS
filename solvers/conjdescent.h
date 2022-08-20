// Copyright (c) 2019 Alex Vlasiuk <oleksandr.vlasiuk@gmail.com>

#ifndef CONJDESCENT_H
#define CONJDESCENT_H

#include <eigen3/Eigen/Core>
#include "../Param.h"
#include "../riesz.h"
#include "../linesearch/LineSearchMT.h"
#include "../linesearch/LineSearch.h"


namespace mpopt {

void normalize(VectorXmp &x, const int n, const int dim)
{
    mpreal temp2;
    for(int i = 0; i < n; i++)  // TODO: this->normalize does not work with mpreals?
    {
        temp2 = mpfr::sqrt(x.segment(i*dim, dim).dot(x.segment(i*dim, dim))); 
        x.segment(i*dim, dim) /= temp2;
    }
}


    //
    // LBFGS solver for unconstrained numerical optimization
    //
    template <typename Scalar>
        class ConjGradSolver
        {
            private:
                typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
                typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
                typedef Eigen::Map<Vector> MapVec;

                const LBFGSParam<Scalar>& m_param;  // Parameters to control the gradient descent algorithm
                Vector        x_old;       // Old x
                Vector        grad_new;    // New gradient
                Vector        grad_old;    // Old gradient
                Vector        drt_mov;     // Moving direction
                Vector        drt_mov_old; // Old moving direction

                inline void reset(int n)
                {
                    x_old.resize(n);
                    grad_new.resize(n);
                    grad_old.resize(n);
                    drt_mov.resize(n);
                    drt_mov_old.resize(n);
                }

            public:
                ConjGradSolver(const LBFGSParam<Scalar>& param) :
                    m_param(param)
            {
                m_param.check_param();
            } 
                template <typename Foo>
                    int minimize(Foo f, Vector& x, Scalar& fout, Scalar& gnorm_return);
        };

    template <typename Scalar> 
        template <typename Foo>
        int ConjGradSolver<Scalar>::minimize(Foo f, Vector& x, Scalar& fout, Scalar& gnorm_return)
        //
        // Polak-Ribiere conjugate gradient.
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
            unsigned n = x.size();
            Scalar xnorm, gnorm, beta;
            fout = f(x, grad_new);
            grad_old.noalias() = grad_new;
            drt_mov_old = Vector::Zero(n);
            for(int i=1; i<=m_param.max_iterations; i++)
            {
                // beta_n in Polak-Ribiere:
                beta = std::max(Scalar(0.0), grad_new.dot(grad_new - grad_old) / grad_old.dot(grad_old));
                // Search direction
                drt_mov.noalias() = -grad_new + beta*drt_mov_old;
                grad_old.noalias() = grad_new;
                // Initial search step
                Scalar step = fout / grad_new.norm();
                x_old.noalias() = x;
                LineSearchMT<Scalar>::MoreThuente(f, fout, x, grad_new, step, drt_mov, x_old, m_param);
                normalize(x, n/3, 3);
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

} 

#endif 

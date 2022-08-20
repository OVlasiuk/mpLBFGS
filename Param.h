// Copyright (c) 2016 Yixuan Qiu
// Copyright (c) 2019 Alex Vlasiuk <oleksandr.vlasiuk@gmail.com>
// ## The MIT License
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:

// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

#ifndef PARAM_H
#define PARAM_H

#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/MPRealSupport>
#include "./mpreal.h"
#include <stdexcept>  // std::invalid_argument

using mpfr::mpreal;


namespace mpopt {


//
//
// Enumeration types for line search.
//

//
//
// The enumeration of line search algorithms.
//
enum LINE_SEARCH_ALGORITHM
{
    //
    // Backtracking method with the Armijo condition.
    // The backtracking method finds the step length such that it satisfies
    // the sufficient decrease (Armijo) condition,
    // f(x + step \cdot drt) \le f(x) + \beta' \cdot step \cdot g(x)^T drt,
    // where x is the current point, drt is the current search direction,
    // step is the step length, and \beta' is the value specified by
    //  LBFGSParam::wolfec1. f and g are the function and gradient values respectively.
    //
    LBFGS_LINESEARCH_BACKTRACKING_ARMIJO = 1,

    //
    // The backtracking method with the defualt (regular Wolfe) condition.
    // An alias of `LBFGS_LINESEARCH_BACKTRACKING_WOLFE`.
    //
    LBFGS_LINESEARCH_BACKTRACKING = 2,

    //
    // Backtracking method with regular Wolfe condition.
    // The backtracking method finds the step length such that it satisfies
    // both the Armijo condition (`LBFGS_LINESEARCH_BACKTRACKING_ARMIJO`)
    // and the curvature condition,
    // g(x + a \cdot d)^T d \ge \beta \cdot g(x)^T d, where \beta
    // is the value specified by  LBFGSParam::wolfec2.
    //
    LBFGS_LINESEARCH_BACKTRACKING_WOLFE = 2,

    //
    // Backtracking method with strong Wolfe condition.
    // The backtracking method finds the step length such that it satisfies
    // both the Armijo condition (`LBFGS_LINESEARCH_BACKTRACKING_ARMIJO`)
    // and the following condition,
    // \vert g(x + a \cdot d)^T d\vert \le \beta \cdot \vert g(x)^T d\vert,
    // where \beta is the value specified by  LBFGSParam::wolfec2.
    //
    LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE = 3
};


//
// Parameters to control the LBFGS algorithm.
//
template <typename Scalar = mpfr::mpreal>
class LBFGSParam
{
public:
    //
    // The number of corrections to approximate the inverse hessian matrix.
    // The L-BFGS routine stores the computation results of previous  m
    // iterations to approximate the inverse hessian matrix of the current
    // iteration. This parameter controls the size of the limited memories
    // (corrections). The default value is  6. Values less than  3 are
    // not recommended. Large values will result in excessive computing time.
    //
    int    m;
    //
    // Tolerance for convergence test.
    // This parameter determines the accuracy with which the solution is to
    // be found. A minimization terminates when
    // ||g|| < \epsilon * \max(1, ||x||),
    // where ||.|| denotes the Euclidean (L2) norm. The default value is
    //  1e-5.
    //
    Scalar epsilon;
    //
    // Distance for delta-based conergence test.
    // This parameter determines the distance d to compute the
    // rate of decrease of the objective function,
    // (f_{k-d}(x)-f_k(x))/f_k(x), where k is the current iteration
    // step. If the value of this parameter is zero, the delta-based convergence
    // test will not be performed. The default value is  0.
    //
    int    past;
    //
    // Delta for convergence test.
    // The algorithm stops when the following condition is met,
    // (f_{k-d}(x)-f_k(x))/f_k(x)<\delta, where f_k(x) is
    // the current function value, f_{k-d}(x) is the function value
    // d iterations ago (specified by the  past parameter).
    // The default value is  0.
    //
    Scalar delta;
    //
    // The maximum number of iterations.
    // The optimization process is terminated when the iteration count
    // exceedes this parameter. Setting this parameter to zero continues an
    // optimization process until a convergence or error. The default value
    // is  0.
    //
    int    max_iterations;
    //
    // The line search algorithm.
    // This parameter specifies the line search algorithm that will be used
    // by the LBFGS routine. The default value is `LBFGS_LINESEARCH_BACKTRACKING_ARMIJO`.
    //
    int    linesearch;
    //
    // The maximum number of trials for the line search.
    // This parameter controls the number of function and gradients evaluations
    // per iteration for the line search routine. The default value is  20.
    //
    int    max_linesearch;
    //
    // The minimum step length allowed in the line search.
    // The default value is  1e-20. Usually this value does not need to be
    // modified.
    //
    Scalar min_step;
    //
    // The maximum step length allowed in the line search.
    // The default value is  1e+20. Usually this value does not need to be
    // modified.
    //
    Scalar max_step;
    //
    // A parameter to control the accuracy of the line search routine.
    // The default value is  1e-4. This parameter should be greater
    // than zero and smaller than  0.5.
    //
    Scalar wolfec1;
    //
    // A coefficient for the Wolfe condition.
    // This parameter is valid only when the backtracking line-search
    // algorithm is used with the Wolfe condition.
    // The default value is  0.9. This parameter should be greater
    // the  wolfec1 parameter and smaller than  1.0.
    //
    Scalar wolfec2;

public:
    //
    // Constructor for LBFGS parameters.
    // Default values for parameters will be set when the object is created.
    //
    LBFGSParam()
    {
        m              = 6;
        epsilon        = Scalar(1e-5);
        past           = 0;
        delta          = Scalar(0);
        max_iterations = 0;
        linesearch     = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
        max_linesearch = 20;
        min_step       = Scalar(1e-8);
        max_step       = Scalar(1e+8);
        wolfec1         = Scalar(1e-4);
        wolfec2         = Scalar(0.9);
    }

    //
    // Checking the validity of LBFGS parameters.
    // An `std::invalid_argument` exception will be thrown if some parameter
    // is invalid.
    //
    inline void check_param() const
    {
        if(m <= 0)
            throw std::invalid_argument("'m' must be positive");
        if(epsilon <= 0)
            throw std::invalid_argument("'epsilon' must be positive");
        if(past < 0)
            throw std::invalid_argument("'past' must be non-negative");
        if(delta < 0)
            throw std::invalid_argument("'delta' must be non-negative");
        if(max_iterations < 0)
            throw std::invalid_argument("'max_iterations' must be non-negative");
        if(linesearch < LBFGS_LINESEARCH_BACKTRACKING_ARMIJO ||
           linesearch > LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE)
           throw std::invalid_argument("unsupported line search algorithm");
        if(max_linesearch <= 0)
            throw std::invalid_argument("'max_linesearch' must be positive");
        if(min_step < 0)
            throw std::invalid_argument("'min_step' must be positive");
        if(max_step < min_step )
            throw std::invalid_argument("'max_step' must be greater than 'min_step'");
        if(wolfec1 <= 0 || wolfec1 >= 0.5)
            throw std::invalid_argument("'wolfec1' must satisfy 0 < wolfec1 < 0.5");
        if(wolfec2 <= wolfec1 || wolfec2 >= 1)
            throw std::invalid_argument("'wolfec2' must satisfy wolfec1 < wolfec2 < 1");
    }
};


}

#endif 

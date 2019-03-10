# mpLBFGS

**mpLBFGS** is a header-only C++ library that implements the limited-memory
BFGS algorithm (L-BFGS) with multiple precision. The code is
derived from [libLBFGS](https://github.com/chokkan/liblbfgs)
and [LBFGS++](http://yixuan.cos.name/LBFGSpp/doc/) libraries; the original
algorithm and implementation goes back to Jorge Nocedal. 
The multiple precision class mpreal uses [MPFR C++](http://www.holoborodko.com/pavel/mpfr/).

**mpLBFGS** is a header-only C++ library, whose only dependency, dependent on [Eigen](http://eigen.tuxfamily.org/), which is also header-only.

## A short example

To use **mpLBFGS**, one needs to first define a functor to represent the
multivariate function to be minimized. It should return the objective function
value on a vector `x` and overwrite the vector `grad` with the gradient
evaluated on `x`. For example we could define the
[Rosenbrock function](https://en.wikipedia.org/wiki/Rosenbrock_function) in the
following way:

```cpp
#include <Eigen/Core>
#include <iostream>
#include <LBFGS.h>

using Eigen::VectorXd;
using namespace LBFGSpp;

class Rosenbrock
{
private:
    int n;
public:
    Rosenbrock(int n_) : n(n_) {}
    double operator()(const VectorXd& x, VectorXd& grad)
    {
        double fx = 0.0;
        for(int i = 0; i < n; i += 2)
        {
            double t1 = 1.0 - x[i];
            double t2 = 10 * (x[i + 1] - x[i] * x[i]);
            grad[i + 1] = 20 * t2;
            grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
            fx += t1 * t1 + t2 * t2;
        }
        return fx;
    }
};
```

Then we just need to set up parameters, create solver object,
provide initial guess, and then run the minimization function.

```cpp
int main()
{
    const int n = 10;
    // Set up parameters
    LBFGSParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 100;

    // Create solver and function object
    LBFGSSolver<double> solver(param);
    Rosenbrock fun(n);

    // Initial guess
    VectorXd x = VectorXd::Zero(n);
    // x will be overwritten to be the best point found
    double fx;
    int niter = solver.minimize(fun, x, fx);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return 0;
}
```

The example can then be compiled and run.

```bash
$ g++ -I/path/to/eigen -I/path/to/lbfgspp/include -O2 example.cpp
$ ./a.out
23 iterations
x =
1 1 1 1 1 1 1 1 1 1
f(x) = 1.87948e-19
```

## License

**mpLBFGS** is licensed under the GPL license.

Copyright (c) 1990, Jorge Nocedal
Copyright (c) 2007–2010, Naoaki Okazaki
Copyright (c) 2016, Yixuan Qiu
Copyright (c) 2018–2019, Alex Vlasiuk 

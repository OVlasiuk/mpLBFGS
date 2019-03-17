#include <iostream>
#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/MPRealSupport>
#include "include/mpreal.h"
#include "include/LBFGS.h"

using Eigen::Matrix;
using Eigen::Dynamic;

using namespace mpfr;
using namespace LBFGSpp;

typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;
typedef Matrix<mpreal,Dynamic,1>        VectorXmp;

mpreal foo(const VectorXmp& x, VectorXmp& grad)
{
    const int n = x.size();
    VectorXmp d(n);
    for(int i = 0; i < n; i++)
        d[i] = i;

    mpreal f = (x - d).squaredNorm();
    grad.noalias() = 2.0 * (x - d);
    return f;
}

int main()
{
    mpreal::set_default_prec(1024);
    const int n = 10;
    LBFGSParam<mpreal> param;
    param.epsilon = mpreal(1e-20);
    LBFGSSolver<mpreal> solver(param);

    VectorXmp x = VectorXmp::Zero(n);
    mpreal fx, gn;
    int niter = solver.minimize(foo, x, fx, gn);

    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return 0;
}

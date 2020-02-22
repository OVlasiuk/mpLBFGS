#include <eigen3/Eigen/Core>
// #include <eigen3/unsupported/Eigen/mprealSupport>
// #include "include/mpreal.h"
#include "include/solvers/conjdescent.h"

using Eigen::Matrix;
using Eigen::Dynamic;

using namespace mpfr;
using namespace mpopt;

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
    ConjGradSolver<mpreal> solver(param);
    param.max_iterations = 1000;

    VectorXmp x = VectorXmp::Zero(n);
    mpreal fx, gn;
    int niter = solver.minimize(foo, x, fx, gn);

    if (niter > 0)
        std::cout << niter << " iterations" << std::endl;
    else
        std::cout << "Maximum number of iterations reached: "<< param.max_iterations << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return 0;
}

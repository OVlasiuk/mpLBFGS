// pframe energy optimization (unconstrained version)
#include <iostream>
#include <eigen3/Eigen/Core>
#include "include/mpreal.h"
#include "include/LBFGS.h"
#include <eigen3/unsupported/Eigen/MPRealSupport>

using namespace mpfr;
using namespace LBFGSpp;

using Eigen::Array;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::RowMajor;

typedef Matrix<mpreal,Dynamic,Dynamic, RowMajor>  MatrixXmp;
typedef Matrix<mpreal,Dynamic,1>        VectorXmp;
typedef Matrix<double,Dynamic,1>        VectorX;
typedef Array<mpreal, Dynamic, Dynamic> ArrayXmp;

class Pframe
{
private:
    const unsigned n;
    const unsigned dim;
    VectorXmp tempp;  //stores pth powers of norms
    VectorXmp temp2;  //stores squares of norms
    mpreal p;
    mpreal f(mpreal t) { return mpfr::pow(mpfr::abs(t), p); };
    mpreal fprime(mpreal t) { return p * t * mpfr::pow(mpfr::abs(t), p-2); };
public:
    Pframe(unsigned n_, unsigned dim_, mpreal p_) : n{n_}, dim{dim_}, p{p_} {
        tempp = VectorXmp::Zero(n);
        temp2 = VectorXmp::Zero(n);
    }
    mpreal operator()(const VectorXmp& x, VectorXmp& grad)
    {
        // store norms in the data member
        for(int i = 0; i < n; i++)
        {
            temp2[i] = x.segment(i*dim, dim).dot(x.segment(i*dim, dim)); // compare with a loop
            tempp[i] = mpfr::pow(temp2[i], p/2.0); 
        }
        // energy
        mpreal fx {0.0}; 
        for(int i = 0; i < n; i++)
            for (int j = 0; j < i; j++)
            {
                fx += f(x.segment(j*dim, dim).dot(x.segment(i*dim, dim))) /tempp[i] /tempp[j];
            }
        fx *= 2.0;
        // gradient
        for (int l = 0; l < dim; l++)
        {
            for(int i = 0; i < n; i++)
            {
                mpreal gr {0};
                for (int j = 0; j<n; j++)
                {
                    mpreal prod = x.segment(j*dim, dim).dot(x.segment(i*dim, dim));
                    gr += (
                            fprime(prod) * x(j*dim+l) - p * x(i*dim+l) * f(prod) / temp2[i]
                          ) / tempp[i] / tempp[j]; 
                }
                grad(i*dim+l) = 2 * gr;
            }
        }
        return fx;
    }
};

int main()
{
    const unsigned n  {56};
    const unsigned dim  {7};
    mpreal p {3};
    mpreal::set_default_prec(mpfr::digits2bits(100));
    VectorXmp grad = VectorXmp::Zero(n*dim);
    std::cout.precision(8);
    std::srand((unsigned int) time(0));

    LBFGSParam<mpreal> param;
    param.epsilon = mpreal(1e-10);
    LBFGSSolver<mpreal> solver(param);

    Pframe fun(n, dim, p);

    // std::cout << "mpfr random\n" << mpfr::random((unsigned int) time(0)) << std::endl;


    // VectorXmp dz = VectorXmp(n*dim);
    // y << 1, 3, 0, 1 ; 
    // dz << 1e-12, 0, 0, 0; 

    // VectorXmp x = VectorXmp::Random(n*dim);
    VectorX y = VectorX::Random(n*dim);
    VectorXmp x = y.cast<mpreal>(); 
    mpreal fx;
    int niter = solver.minimize(fun, x, fx);

    Eigen::Map<MatrixXmp> M(x.data(), n,dim);

    // .normalize does not work with mpreals?
    VectorXmp temp2 = VectorXmp::Zero(n);
    for(int i = 0; i < n; i++)
    {
        temp2[i] = mpfr::sqrt(x.segment(i*dim, dim).dot(x.segment(i*dim, dim))); 
    }
    for(int i = 0; i < n; i++)
        M.row(i) /= temp2[i];
    //
    
    // std::cout << "M:" << std::endl << M << std::endl; 
    //
    // VectorXmp G = M.rowwise().norm(); 
    // std::cout << "x:" << std::endl << x << std::endl;
    //
    std::cout << "Gram matrix = \n"  <<  M*M.transpose() << std::endl;
    std::cout << niter << " iterations" << std::endl;
    std::cout << "f(x) = " << fx << std::endl;
    std::cout <<  "Gradient norm \n" << grad.norm()  << std::endl; 

    return 0;
}





    // // Tests
    // VectorX y = VectorX(n*dim);
    // VectorXmp dz = VectorXmp(n*dim);
    //     // ::Random(n*dim);
    // y << 1, 3, 0, 1 ; 
    // dz << 1e-12, 0, 0, 0; 

    // VectorXmp z = y.cast<mpreal>(); 
    // dz = dz.cast<mpreal>(); 
    // for (int i=0; i<2; i++)
    //     z(i) = mpfr::sqrt(z(i)); 
    // std::cout << "Vector value \n" << z <<  std::endl;
    // // std::cout << "evaluation at y \n" << fun(y, grad) <<  std::endl;
    // std::cout << "Evaluation at z:\t" << fun(z, grad) <<  std::endl;
    // std::cout << "Gradient at z \n" << grad <<  std::endl;
    // std::cout << "Numerical gradient at z \n" << (fun(z+dz, grad)-fun(z, grad))/dz.norm() <<  std::endl;
    // // End tests

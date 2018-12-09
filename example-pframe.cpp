// pframe energy optimization (unconstrained version)
#include <iostream>
#include "include/mpreal.h"
#include "include/LBFGS.h"
#include <eigen3/Eigen/Core>
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
    mpreal f(mpreal t) { return mpfr::pow(mpfr::abs(t), p); };  // p-th power of the absolute value
    mpreal fprime(mpreal t) { return p * t * mpfr::pow(mpfr::abs(t), p-2); };   // derivative of f
public:
    Pframe(unsigned n_, unsigned dim_, mpreal p_) : n{n_}, dim{dim_}, p{p_} {
        tempp = VectorXmp::Zero(n);
        temp2 = VectorXmp::Zero(n);
    }
    mpreal operator()(const VectorXmp& x, VectorXmp& grad)
    {
        // Store norms in the data members
        for(int i = 0; i < n; i++)
        { 
            temp2[i] = x.segment(i*dim, dim).dot(x.segment(i*dim, dim));    // squared norms of the configuration 
            tempp[i] = mpfr::pow(temp2[i], p/2.0);  // p-th powers of norms of the configuration
        }
        // Energy
        mpreal fx {n/2.0}; 
        for(int i = 0; i < n; i++)
            for (int j = 0; j < i; j++)
            {
                fx += f(x.segment(j*dim, dim).dot(x.segment(i*dim, dim))) /tempp[i] /tempp[j];
            }
        fx *= 2.0 / n / n;
        // Gradient
        grad = VectorXmp::Zero(n*dim);
        for (int l = 0; l < dim; l++)   // iterate over coordinates
        {
            for(int i = 0; i < n; i++)  // iterate over vectors
            {
                mpreal gr {0};
                for (int j = 0; j<n; j++)
                {
                    mpreal prod = x.segment(j*dim, dim).dot(x.segment(i*dim, dim));
                    gr += (
                            fprime(prod) * x(j*dim+l) - p * x(i*dim+l) * f(prod) / temp2[i]
                          ) / tempp[i] / tempp[j]; 
                }
                grad(i*dim+l) = 2 * gr / n / n; // "2" to account for i-th column and i-th row
            }
        }
        return fx;
    }
};

int main(int argc, char *argv[])
{
    // Initialize scalar parameters
    unsigned n  {12};
    unsigned dim  {3};
    mpreal p {3};
    std::string epsilon {"20"};
    mpreal fx, gn;
    switch (argc) {
        case 5 ... INT_MAX: // only gcc and clang-compatible
            epsilon = argv[4];
        case 4: // only gcc and clang-compatible
            p = atof(argv[3]);
        case 3:
            dim = atoi(argv[2]);
        case 2:
            n = atoi(argv[1]);
    }
    std::cout << "n = " << n <<  std::endl;
    std::cout << "dim = " << dim <<  std::endl;
    std::cout << "p = " << p <<  std::endl;
    std::cout << "-log10(epsilon) = " << epsilon <<  std::endl;

    // Initialize vector parameters
    VectorX y(n*dim);
    VectorXmp x(n*dim);

    // Initialize system parameters
    mpreal::set_default_prec(mpfr::digits2bits(100));
    std::srand((unsigned int) time(0));

    // Initialize the solver class
    LBFGSParam<mpreal> param;
    param.epsilon = mpreal("1e-"+epsilon);
    LBFGSSolver<mpreal> solver(param); 
    Pframe fun(n, dim, p);

    // std::cout << "mpfr random\n" << mpfr::random((unsigned int) time(0)) << std::endl; 
    y = VectorX::Random(n*dim); // random starting config
    x += y.cast<mpreal>();  // TODO: VectorXmp::Random does not work?  

    int niter = solver.minimize(fun, x, fx, gn);
    //

    Eigen::Map<MatrixXmp> M(x.data(), n,dim); 
    VectorXmp temp2 = VectorXmp::Zero(n);
    for(int i = 0; i < n; i++)  // TODO: this->normalize does not work with mpreals?
    {
        temp2[i] = mpfr::sqrt(x.segment(i*dim, dim).dot(x.segment(i*dim, dim))); 
    }
    for(int i = 0; i < n; i++)
        M.row(i) /= temp2[i];
    //
    
    // std::cout << "Gram matrix = \n"  <<  M*M.transpose() << std::endl;
    std::cout << niter << " iterations" << std::endl;
    std::cout.precision(24);
    std::cout << "f(x) = " << fx << std::endl;
    std::cout <<  "Gradient norm \n" << gn  << std::endl; 

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

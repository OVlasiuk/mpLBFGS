#include "include/mpreal.h"
#include "include/LBFGS.h"
#include "include/descent.h"
#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/MPRealSupport>

using Eigen::Array;
using Eigen::Matrix;
using Eigen::Dynamic;
using Eigen::RowMajor;

typedef Matrix<mpreal,Dynamic,Dynamic, RowMajor>  MatrixXmp;
typedef Matrix<mpreal,Dynamic,1>        VectorXmp;
typedef Matrix<double,Dynamic,1>        VectorX;
typedef Array<mpreal, Dynamic, Dynamic> ArrayXmp;

class wtPframe
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
    wtPframe(unsigned n_, unsigned dim_, mpreal p_) : n{n_}, dim{dim_}, p{p_} {
        tempp = VectorXmp::Zero(n);
        temp2 = VectorXmp::Zero(n);
    }
    mpreal operator()(const VectorXmp& x, VectorXmp& grad)
    {
        // Store norms in the data members
        for(int i = 0; i < n; i++)
        { 
            temp2[i] = x.segment(i*dim, dim).dot(x.segment(i*dim, dim));    // squared norms of the configuration 
            tempp[i] = mpfr::pow(temp2[i], p/mpreal("2.0"));  // p-th powers of norms of the configuration
        }
        // Energy
        mpreal fx {mpreal("0.0")}; 
        for(int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
            {
                fx += x[n*dim+i]* x[n*dim+i] * x[n*dim+j] * x[n*dim+j] *
                    f(x.segment(j*dim, dim).dot(x.segment(i*dim, dim))) /tempp[i] /tempp[j];
            }
        // fx *= mpreal("2.0");
        // Gradient: coordinate terms
        grad = VectorXmp::Zero(n*(dim+1));
        for (int l = 0; l < dim; l++)   // iterate over coordinates
        {
            for(int i = 0; i < n; i++)  // iterate over vectors
            {
                mpreal gr  {"0"};
                for (int j = 0; j<n; j++)
                {
                    mpreal prod = x.segment(j*dim, dim).dot(x.segment(i*dim, dim));
                    gr += (
                            fprime(prod) * x(j*dim+l) - p * x(i*dim+l) * f(prod) / temp2[i]
                          ) / tempp[i] / tempp[j] * x[n*dim+j] * x[n*dim+j]; 
                }
                grad(i*dim+l) = mpreal("2") * gr * x[n*dim+i]* x[n*dim+i]; // "2" to account for i-th column and i-th row
            }
        }
        // Gradient: weight terms
        for(int i = 0; i < n; i++)  // iterate over vectors
        {
            mpreal gr  {"0"};
            for (int j = 0; j<n; j++)
            {
                gr +=  f(x.segment(j*dim, dim).dot(x.segment(i*dim, dim))) /tempp[i] /tempp[j] * x[n*dim+j] * x[n*dim+j]; 
            }
            grad(n*dim+i) = mpreal("4") * gr * x[n*dim+i]; // "4" to account for i-th column and i-th row
        }
        return fx;
    }

    void testgrad(VectorXmp x0)
    {
        // Test gradient
        VectorXmp   h(x0.size()), z(x0.size()), y0(x0.size()), y1(x0.size()); 
        mpreal hterm = 1e-12;
        for (int j = 0; j < x0.size(); j++)
        {
            h[j] = hterm;
            z[j] = (this->operator()(x0+h, y1) - this->operator()(x0, y0))/hterm;
            h[j] = 0.0;
        }
        std::cout << (z - y0).norm() << std::endl;
    }
};

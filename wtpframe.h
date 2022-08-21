#include "include/mpreal.h"
#include "solvers/LBFGS.h"
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

void parseinput(const int __argc, char* const * __argv, unsigned& n, unsigned& dim, double& p, std::string& epsilon, int& paramm, long& inits, bool& gram, bool& coord, bool& stage2, bool& help);

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
        mpreal fx {"0.0"}, wtsum {"0.0"}, jterms; 
        for(int i = 0; i < n; i++)
        {
            jterms = 0.0;
            for (int j = 0; j < n; j++)
            {
                jterms += x[n*dim+j] * x[n*dim+j] *
                    f(x.segment(j*dim, dim).dot(x.segment(i*dim, dim))) / tempp[j];
            }
            jterms *= x[n*dim+i] * x[n*dim+i] / tempp[i];  
            fx += jterms;
            wtsum += x[n*dim+i]* x[n*dim+i];
        }
        fx /= wtsum*wtsum;

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
                          ) / tempp[j] * x[n*dim+j] * x[n*dim+j]; 
                }
                gr /= tempp[i] * wtsum * wtsum;
                grad(i*dim+l) = mpreal("2") * gr * x[n*dim+i]* x[n*dim+i]; // "2" to account for i-th column and i-th row
            }
        }
        // Gradient: weight terms
        for(int i = 0; i < n; i++)  // iterate over vectors
        {
            mpreal gr  {"0"};
            for (int j = 0; j<n; j++)
            {
                gr += f(x.segment(j*dim, dim).dot(x.segment(i*dim, dim)))  /tempp[j] *  x[n*dim+j] * x[n*dim+j]; 
            }
            gr /= tempp[i] * wtsum * wtsum;
            gr -= fx /wtsum; // save flops by reusing the energy value 
            grad(n*dim+i) = 4 * x[n*dim+i] * gr; 
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
        std::cout << "Discrepancy between the numerical and given gradients:" << std::endl;
        std::cout << (z - y0).norm() << std::endl;
    }
};

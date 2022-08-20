// weighted pframe energy optimization 
#include <iostream>
#include <fstream>
#include <random>
#include <getopt.h>
#include "wtpframe.h"

using namespace mpfr;
using namespace mpopt;


void normalize(VectorXmp &x, const int n, const int dim)
{
    mpreal temp2;
    for(int i = 0; i < n; i++)  // TODO: this->normalize does not work with mpreals?
    {
        temp2 = mpfr::sqrt(x.segment(i*dim, dim).dot(x.segment(i*dim, dim))); 
        x.segment(i*dim, dim) /= temp2;
        x[n*dim+i] = mpfr::abs(x[n*dim+i]);
    }
    // x.segment(n*dim, n) /= x.segment(n*dim, n).dot(x.segment(n*dim, n));
}


int main(int argc, char *argv[])
{
    // Optimization parameters
    unsigned n  {12};
    unsigned dim  {3};
    double p {3};
    std::string epsilon {"20"}, epsilon_step {"40"};
    long inits {1}; 
    const unsigned history {5};
    // Output parameters
    mpreal fx {mpfr::const_infinity()}, gx {mpfr::const_infinity()}, gn;
    int niter; 
    // Solver parameter class
    LBFGSParam<mpreal> param;
    param.m = 6;
    param.max_iterations = 10000;
    // Output options
    bool gram {false}, coord {false}, help {false}, stage2 {true};

    // // Parse command line input
    parseinput(argc, argv, n, dim, p, epsilon, param.m, inits, gram, coord, stage2, help);
    if (help)
        return 0;
    // Done parsing command line input


    // Initialize vector parameters
    VectorXmp x(n*(dim+1));
    MatrixXmp X(n*(dim+1), history);
    VectorX y(n*(dim+1));

    // Initialize library parameters
    mpreal::set_default_prec(mpfr::digits2bits(100));

    // Initialize random device
    std::random_device rd{};
    std::mt19937 gen{rd()}; // Mersenne twister
    std::normal_distribution<double> d{0,1}; 


    // std::cout << "mpfr random\n" << mpfr::random((unsigned int) time(0)) << std::endl; 
    if (stage2)
    {
        // Initialize the solver class
        param.epsilon = mpfr::sqrt(mpreal("1e-"+epsilon));
        param.min_step = mpreal("1e-"+epsilon_step);
        LBFGSSolver<mpreal> solver(param); 
        wtPframe fun(n, dim, p);
        unsigned hist = 0;
        for (int i = 0; i < inits; i++)
        {
            mpreal f0;
            VectorXmp  x0(n*(dim+1));
            for (int j = 0; j < n*(dim+1); j++)
                x0[j] = d(rd);
            normalize(x0, n, dim);

            niter = solver.minimize(fun, x0, f0, gn);
            if (f0 < fx)
            {
                fx = f0;
                X.col(hist++ % history) = x0;
            }
            std::cout << niter << " iterations" << std::endl;
        }

        std::cout << std::endl << "Stage 2 optimization using " << std::min(hist, history) << " lowest values from stage 1:" << std::endl;

        // New solver class
        param.epsilon = mpreal("1e-"+epsilon);
        param.min_step = mpreal("1e-"+epsilon_step);
        LBFGSSolver<mpreal> solver1(param); 
        fx = mpfr::const_infinity();
        for (int i = 0; i<std::min(hist, history); i++)
        {
            mpreal f0; 
            VectorXmp  x0(n*(dim+1));
            x0 = X.col(i); 
            normalize(x0, n, dim);

            niter = solver1.minimize(fun, x0, f0, gn);
            if (f0 < fx)
            {
                fx = f0;
                x = x0;
                gx = gn;
            }
            std::cout << niter << " iterations" << std::endl;
        }
    }
    else
    {
        // Initialize the solver class
        param.epsilon = mpreal("1e-"+epsilon);
        param.min_step = mpreal("1e-"+epsilon_step);
        LBFGSSolver<mpreal> solver(param); 
        wtPframe fun(n, dim, p);
        for (int i = 0; i < inits; i++)
        {
            mpreal f0;
            VectorXmp  x0(n*(dim+1));
            for (int j = 0; j < n*(dim+1); j++)
                x0[j] = d(rd);
            normalize(x0, n, dim);

            niter = solver.minimize(fun, x0, f0, gn);
            if (f0 < fx)
            {
                fx = f0;
                x = x0;
                gx = gn;
            }
            std::cout << niter << " iterations" << std::endl;
        }

    }

    normalize(x, n, dim);
    for(int i=0; i < n; i++)
        x.tail(n)[i] *= x.tail(n)[i];
    x.segment(n*dim, n) /= x.segment(n*dim, n).sum();
    Eigen::Map<MatrixXmp> M(x.head(n*dim).data(), n,dim), W(x.tail(n).data(), n,1); 
    
    std::cout.precision(24);
    std::cout << "f(x) = " << fx << std::endl;
    std::cout <<  "Gradient norm \n" << gx  << std::endl; 

    char pstring[10];
    sprintf(pstring, "%.3f", p);

    if ( gram )
    {
        std::string filename;
        filename =  "grammatrix_" 
            + std::to_string(dim) + "_"
            + std::to_string(n) + "_"
            + std::string(pstring) + ".txt";
        std::ofstream ostrm;
        ostrm.open(filename);
        ostrm.precision(24);
        ostrm << M*M.transpose() << '\n'; 
        std::cout << "\nThe gram matrix is written to file\n  "  << filename  << std::endl;
    }
    if ( coord )
    {
        std::string filename, wtname;
        filename =  "coords_" 
            + std::to_string(dim) + "_"
            + std::to_string(n) + "_"
            + std::string(pstring) + ".txt";
        wtname =  "weights_" 
            + std::to_string(dim) + "_"
            + std::to_string(n) + "_"
            + std::string(pstring) + ".txt";
        std::ofstream ostrm;
        ostrm.open(filename);
        ostrm.precision(24);
        ostrm << M << '\n'; 
        ostrm.close();

        ostrm.open(wtname);
        ostrm.precision(24);
        ostrm << W << '\n'; 
        std::cout << "Coordinates written to file\n  "  << filename  << std::endl;
        std::cout << "Weights written to file\n  "  << wtname  << std::endl;
        ostrm.close();
    }

    return 0;
}

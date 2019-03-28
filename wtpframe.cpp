// pframe energy optimization (unconstrained version)
#include <iostream>
#include <fstream>
#include "wtpframe.h"

using namespace mpfr;
using namespace LBFGSpp;

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
    // Initialize scalar parameters
    unsigned n  {12};
    unsigned dim  {3};
    mpreal p {"3"};
    std::string epsilon {"20"}, epsilon_step {"20"};
    long inits {1};
    mpreal fx {mpfr::const_infinity()}, gn;
    int niter;
    LBFGSParam<mpreal> param;
    param.m = 6;
    param.max_iterations = 5000;
    switch (argc) {
        case 7 ... INT_MAX: // only gcc and clang-compatible
            inits = atol(argv[6]);
        case 6: 
            param.m = atoi(argv[5]);
        case 5: 
            epsilon = argv[4];
        case 4: 
            p = atof(argv[3]);
        case 3:
            dim = atoi(argv[2]);
        case 2:
            n = atoi(argv[1]);
    }
    std::cout << "Parameter values:" << std::endl;
    std::cout << "n = " << n <<  "\t\t" << "dim = " << dim <<  "\t\t" << "p = " << p <<  "\t\t" << "-log10(epsilon) = " << epsilon <<  std::endl;
    std::cout << "Performing initializations: " << inits <<  std::endl << std::endl;
    if (argc<2)
    {
        std::cout << "This set of parameters corresponds to the following command line input: " << std::endl;
        std::cout << "./pframe " << n << " " << dim << " " << p << " " << epsilon << " " << param.m << " " << inits  << std::endl ;
        std::cout <<  "n" << " " << "dim" << " " << "p" << " " << "epsilon" << " " << "param.m" << " " << "inits"  << std::endl << std::endl;
    }
    param.epsilon = mpreal("1e-"+epsilon);
    param.min_step = mpreal("1e-"+epsilon_step);


    // Initialize vector parameters
    VectorXmp x(n*(dim+1));
    VectorX y(n*(dim+1));

    // Initialize system parameters
    mpreal::set_default_prec(mpfr::digits2bits(100));
    std::srand((unsigned int) time(0));

    // Initialize the solver class
    LBFGSSolver<mpreal> solver(param); 
    wtPframe fun(n, dim, p);

    // std::cout << "mpfr random\n" << mpfr::random((unsigned int) time(0)) << std::endl; 
    for (int i = 0; i < inits; i++)
    {
        mpreal f0;
        VectorXmp  x0(n*(dim+1));
        y = VectorX::Random(n*(dim+1)); // random starting config
        x0 = y.cast<mpreal>();  // TODO: VectorXmp::Random does not work?  
        normalize(x0, n, dim);

        niter = solver.minimize(fun, x0, f0, gn);
        if (f0 < fx)
        {
            fx = f0;
            x = x0;
        }
        std::cout << niter << " iterations" << std::endl;
    }

    normalize(x, n, dim);
    for(int i=0; i < n; i++)
        x.tail(n)[i] *= x.tail(n)[i];
    x.segment(n*dim, n) /= x.segment(n*dim, n).sum();
    Eigen::Map<MatrixXmp> M(x.head(n*dim).data(), n,dim), W(x.tail(n).data(), n,1); 
    
    std::cout.precision(24);
    std::cout << "f(x) = " << fx << std::endl;
    std::cout <<  "Gradient norm \n" << gn  << std::endl; 
    std::cout <<  "Output the Gram matrix y/N? \n" ; 
    char yn;
    std::cin >> yn;
    if (yn == 'y')
    {
        char filename[] {"grammatrix.txt"};
        std::ofstream ostrm;
        ostrm.open(filename);
        ostrm.precision(24);
        ostrm << M*M.transpose() << '\n'; 
        std::cout << "Written to file "  << filename  << std::endl;
    }
    std::cout <<  "Output coordinates and weights y/N? \n" ; 
    std::cin >> yn;
    if (yn == 'y')
    {
        char filename[] {"coords.txt"}, wtname[] {"weights.txt"};
        std::ofstream ostrm;
        ostrm.open(filename);
        ostrm.precision(24);
        ostrm << M << '\n'; 
        ostrm.close();

        ostrm.open(wtname);
        ostrm.precision(24);
        ostrm << W << '\n'; 
        std::cout << "Coordinates written to file "  << filename  << std::endl;
        std::cout << "Weights written to file "  << wtname  << std::endl;
        ostrm.close();
    }

    return 0;
}


        ////
        //// for (int j = 0; j < z.size(); j++)
        ////     std::cout << z[j]  << '\t' << y0[j] << std::endl;
        //// Test gradient
        ////
        //return 0;



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

// weighted pframe energy optimization 
#include <iostream>
#include <fstream>
#include <getopt.h>
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
    // Optimization parameters
    unsigned n  {12};
    unsigned dim  {3};
    mpreal p {"3"};
    std::string epsilon {"20"}, epsilon_step {"20"};
    long inits {1}; 
    // Output parameters
    mpreal fx {mpfr::const_infinity()}, gn;
    int niter; 
    // Solver parameter class
    LBFGSParam<mpreal> param;
    param.m = 6;
    param.max_iterations = 5000;
    // Output options
    bool gram {false}, coord {false}, help {false};


    // Parse command line input
    int c; 
    while (1) {
        int option_index = 0;
        static struct option long_options[] = {
            {"n"       , required_argument , 0 , 'n' } ,
            {"dim"     , required_argument , 0 , 'd' } ,
            {"p"       , required_argument , 0 , 'p' } ,
            {"epsilon" , required_argument , 0 , 'e' } ,
            {"m"       , required_argument , 0 , 'm' } ,
            {"inits"   , required_argument , 0 , 'i' } ,
            {"gram"    , no_argument       , 0 , 'g' } ,
            {"coord"   , no_argument       , 0 , 'c' } ,
            {"help"    , no_argument       , 0 , 'h' } ,
            {0         , 0                 , 0 , 0   }
        };

        c = getopt_long(argc, argv, "n:d:p:e:m:i:gch",
                long_options, &option_index);
        if (c == -1)
            break;

        switch (c) {
            case 0:
                printf("option %s", long_options[option_index].name);
                if (optarg)
                    printf(" with arg %s", optarg);
                printf("\n");
                break;

            case 'n':
                n = atoi(optarg); 
                break; 
            case 'd':
                dim = atoi(optarg);
                break; 
            case 'p':
                p = atof(optarg);
                break; 
            case 'e':
                epsilon = optarg;
                break; 
            case 'm':
                param.m = atoi(optarg);
                break; 
            case 'i':
                inits = atol(optarg);
                break; 
            case 'g':
                gram = true;
                break; 
            case 'c':
                coord = true;
                break; 
            case 'h':
                help = true;
                break; 
            case '?':
                break; 
            default:
                printf("?? getopt returned character code 0%o ??\n", c);
        }
    }
    if (optind < argc) {
        printf("Non-recognized arguments:");
        while (optind < argc)
            printf("%s ", argv[optind++]);
        printf("\n");
    }
    char pstring[10];
    sprintf(pstring, "%.3f",double(p));
    // Done parsing command line input


    // Talk to the user
    std::cout << "Parameter values:" << std::endl;
    std::printf("n = %d\t\tdim = %d \tp = %s\
            \n-log10(epsilon) = %s\t\tparam.m = %d\n",
            n, dim, pstring, 
            epsilon.c_str(), param.m);
    // std::cout << "n = " << n <<  "\t\t" << "dim = " << dim <<  "\t\t" << "p = " << p <<  "\t\t" << "-log10(epsilon) = " << epsilon <<  std::endl;
    std::cout << "Initializations: " << inits  << std::endl;
    std::cout << "Write gram matrix: " << (gram ?  "Yes" : "No") <<   std::endl;
    std::cout << "Write coordinate matrix: " << (coord ?  "Yes" : "No") << std::endl <<  std::endl;

    if (argc<2 || help)
    {
        std::cout << "This set of parameters corresponds to the following command line input: " << std::endl;
        std::printf("./pframe -n %d -d %d -p %s -e %s -m %d -i %ld %s%s\n",
                n, dim, pstring, epsilon.c_str(), param.m, inits,
                (gram ?  "-g " : ""), (coord ?  "-c" : ""));
        std::printf("Long options are also supported:\n\
                --n       --dim     --p      \n\
                --epsilon --m       --inits  \n\
                --gram    --coord   --help\n");
    }
    if (help)
        return 0;
    // Done talking to the user

    param.epsilon = mpreal("1e-"+epsilon);
    param.min_step = mpreal("1e-"+epsilon_step);


    // Initialize vector parameters
    VectorXmp x(n*(dim+1));
    VectorX y(n*(dim+1));

    // Initialize library parameters
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

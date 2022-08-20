#include <iostream>
#include <getopt.h>
#include "wtpframe.h"

// (int ___argc, char *const *___argv, const char *__shortopts, const struct
// option *__longopts, int *__longind)

static const char short_options[]="1n:d:p:e:m:i:gch";
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

void parseinput(const int __argc, char* const * __argv, unsigned& n, unsigned& dim, double& p, std::string& epsilon, int& paramm, long& inits, bool& gram, bool& coord, bool& stage2, bool& help)
{ 
    // Parse command line input
    int c; 
    int option_index = 0;
    while ((c = getopt_long(__argc, __argv, short_options, long_options, &option_index)) != -1) { 
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
                paramm = atoi(optarg);
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
            case '1':
                stage2 = false;
                break; 
            case 'h':
                help = true;
                break; 
            case '?':
                break; 

            default:
                printf("?? getopt returned character code 0%o ??\n", c);
        }
    }; 
    
    if (optind < __argc) {
        printf("Non-recognized arguments:");
        while (optind < __argc)
            printf("%s ", __argv[optind++]);
        printf("\n");
    }

    // Talk to the user
    std::cout << "Parameter values:" << std::endl;
    std::printf("n = %d\t\tdim = %d \tp = %.3f\
            \n-log10(epsilon) = %s\t\tparam.m = %d\n",
            n, dim, p, 
            epsilon.c_str(), paramm);
    // std::cout << "n = " << n <<  "\t\t" << "dim = " << dim <<  "\t\t" << "p = " << p <<  "\t\t" << "-log10(epsilon) = " << epsilon <<  std::endl;
    std::cout << "Initializations: " << inits  << std::endl;
    std::cout << "2-stage optimization: " << (stage2 ?  "Yes" : "No") <<   std::endl;
    std::cout << "Write gram matrix: " << (gram ?  "Yes" : "No") <<   std::endl;
    std::cout << "Write coordinate matrix: " << (coord ?  "Yes" : "No") << std::endl <<  std::endl;

    if (__argc<2 || help)
    {
        std::cout << "This set of parameters corresponds to the following command line input: " << std::endl;
        std::printf("./pframe -n %d -d %d -p %.3f -e %s -m %d -i %ld %s%s\n",
                n, dim, p, epsilon.c_str(), paramm, inits,
                (gram ?  "-g " : ""), (coord ?  "-c" : ""));
        std::printf("Long options are also supported:\n\
                --n       --dim     --p      \n\
                --epsilon --m       --inits  \n\
                --gram    --coord   --help\n");
    }
    // Done talking to the user
}

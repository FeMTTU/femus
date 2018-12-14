#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "petsc.h"
#include "petscmat.h"
#include "PetscMatrix.hpp"

#include "slepceps.h"
#include "sparseGrid.hpp"

//THIS EXTENDS WHAT IS IN EX10 TO SPARSE GRIDS

using namespace femus;



//BEGIN stochastic data

unsigned alpha = 1;
unsigned M = pow (10, alpha); //number of samples
unsigned N = 2; //dimension of the parameter space (each of the M samples has N entries)

//FOR NORMAL DISTRIBUTION
boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)
boost::normal_distribution<> nd (0., 0.3);
boost::variate_generator < boost::mt19937&,
      boost::normal_distribution<> > var_nor (rng, nd);

//FOR UNIFORM DISTRIBUTION
boost::mt19937 rng1; // I don't seed it on purpouse (it's not relevant)
boost::random::uniform_real_distribution<> un (- 1., 1.);
boost::variate_generator < boost::mt19937&,
      boost::random::uniform_real_distribution<> > var_unif (rng1, un);

//FOR LAPLACE DISTRIBUTION
boost::mt19937 rng2; // I don't seed it on purpouse (it's not relevant)
boost::random::uniform_real_distribution<> un1 (- 0.5, 0.49999999999);
boost::variate_generator < boost::mt19937&,
      boost::random::uniform_real_distribution<> > var_unif1 (rng2, un1);
double b = 2.;
//END

int main (int argc, char** argv)
{

    std::vector < std::vector < double > >  samples;
    samples.resize (M);
    for (unsigned m = 0; m < M; m++) {
        samples[m].resize (N);
    }


    for (unsigned m = 0; m < M; m++) {
        for (unsigned n = 0; n < N; n++) {

            double var = var_nor();
            double varunif = var_unif();
            double U = var_unif1();
//     samples[m][n] = var * var * var;
//     samples[m][n] = exp(var);
//         samples[m][n] = exp (varunif);
            samples[m][n] = varunif;

            //exp of truncated gaussian
//     if(fabs(var) <= 1.) {
//       samples[m][n] = var / (0.5 * ((1. + erf((1. / 0.3) / sqrt(2))) - (1. + erf((- 1. / 0.3) / sqrt(2))))) ;    //truncated Gaussian
//     }
//     else samples[m][n] = 0.;

            //laplace distribution
//     double signU = 0.;
//     if(U < 0) signU = - 1.;
//     else if(U > 0) signU = 1.;
//     samples[m][n] = 0. - b * signU * log(1. - 2. * fabs(U)) ;

//     std::cout << "samples[" << m << "][" << n << "]=" << samples[m][n] << std::endl;

        }
    }


    sparseGrid spg (N, M, samples);
    
    double phi;
    double x = 0.;
    unsigned nn = 0;
    unsigned ll = 1;
    unsigned ii = 0;
    bool scale = false;
    spg.EvaluateOneDimensionalPhi(phi, x, nn, ll, ii, scale);
    
    std::cout << "phi = " << phi << std::endl;

    return 0;

} //end main




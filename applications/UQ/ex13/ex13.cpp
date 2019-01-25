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


#include <vector>


using namespace femus;



//BEGIN stochastic data

//bool sparse = true;
bool sparse = false;

unsigned alpha = 7;
unsigned M = pow(10, alpha);    //number of samples
unsigned N = 2; //dimension of the parameter space (each of the M samples has N entries)
unsigned L = 3; //max refinement level
bool output = false; //for debugging
bool matlabView = true;

double xmin = - 5.5;   //-1.5 for uniform // -5.5 for Gaussian
double xmax = 5.5;     //1.5 for uniform // 5.5 for Gaussian

//FOR NORMAL DISTRIBUTION
boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)
boost::normal_distribution<> nd(0., 1.);
boost::variate_generator < boost::mt19937&,
      boost::normal_distribution<> > var_nor(rng, nd);

//FOR UNIFORM DISTRIBUTION
boost::mt19937 rng1; // I don't seed it on purpouse (it's not relevant)
boost::random::uniform_real_distribution<> un(- 1., 1.);
boost::variate_generator < boost::mt19937&,
      boost::random::uniform_real_distribution<> > var_unif(rng1, un);

//FOR LAPLACE DISTRIBUTION
boost::mt19937 rng2; // I don't seed it on purpouse (it's not relevant)
boost::random::uniform_real_distribution<> un1(- 0.5, 0.49999999999);
boost::variate_generator < boost::mt19937&,
      boost::random::uniform_real_distribution<> > var_unif1(rng2, un1);
double b = 2.;
//END

void PrintBaseCoefficients(const std::vector < double > & c, const std::vector <unsigned> &I, const std::vector < unsigned > &n);

int main(int argc, char** argv) {

  //BEGIN construction of the sample set
  std::vector < std::vector < double > >  samples;
  samples.resize(M);

  for (unsigned m = 0; m < M; m++) {
    samples[m].resize(N);
  }

  for (unsigned m = 0; m < M; m++) {
    for (unsigned n = 0; n < N; n++) {
      double var = var_nor();
      double varunif = var_unif();
      double U = var_unif1();
      samples[m][n] = var;
      //samples[m][n] = varunif;
    }
  }
  //END


  std::vector< unsigned > n1D(L + 1);
  n1D[0] = 1u;
  for (unsigned i = 1; i < n1D.size(); i++) {
    n1D[i] = n1D[i - 1] * 2u;
  }

  //Build Histogram
  std::vector <double> cI(static_cast< unsigned >(pow(n1D[L], N)));

  double h = (xmax - xmin) / n1D[L];
  double h2 = h * h;

  for (unsigned m = 0; m < M; m++) {

    std::vector < unsigned > I(N);
    for (unsigned k = 0; k < N; k++) {
      double x = (samples[m][k] - xmin) / h;
      if (x < 0.) I[k] = 0;
      else if (x >= n1D[L]) I[k] = n1D[L] - 1;
      else I[k] = static_cast < unsigned >(floor(x));
    }

    unsigned INDEX = I[0];
    for (unsigned k = 1; k < N; k++) {
      INDEX = INDEX * n1D[L] + I[k];
    }
    cI[INDEX] += 1. / (M * h2);
  }

  std::vector <unsigned> IL(N, L);
  std::cout << "Histogram " << std::endl;
  PrintBaseCoefficients(cI, IL, n1D);

  ///////////////////////////////////////////////////////////////////////////////

  //Build Hierarchical Bases

  std::vector < std::vector <std::vector <double> > > CH(L + 1);
  for (unsigned iL = 0; iL < CH.size(); iL++) {
    CH[iL].resize(L + 1);
    for (unsigned jL = 0; iL * sparse + jL < CH[iL].size(); jL++) {
      CH[iL][jL].resize(n1D[iL] * n1D[jL]);
    }
  }

  double Hx = (xmax - xmin);
  double Hy = (xmax - xmin);

  for (unsigned m = 0; m < M; m++) {
    double X = (samples[m][0] - xmin) / Hx;
    double Y = (samples[m][1] - xmin) / Hy;

    if (X < 0.)  X = 0.;
    if (Y >= 1.) X = 0.9999999999999999999999999999999999;

    if (Y < 0.)  Y = 0.;
    if (Y >= 1.) Y = 0.9999999999999999999999999999999999;

    for (unsigned iL = 0; iL < CH.size(); iL++) {
      double hx = Hx / n1D[iL];
      double x = X * n1D[iL];
      for (unsigned jL = 0; iL * sparse + jL < CH[iL].size(); jL++) {
        double hy = Hy / n1D[jL];
        double y = Y * n1D[jL];
        unsigned i = static_cast < unsigned >(floor(x));
        unsigned j = static_cast < unsigned >(floor(y));
        CH[iL][jL][i * n1D[jL] + j] += 1. / (M * hx * hy);
      }
    }
  }

  
  std::vector < unsigned > JL(N);
  for (unsigned iL = 0; iL < CH.size(); iL++) {
    for (unsigned jL = 0; iL * sparse + jL < CH[iL].size(); jL++) {
      int iL1 = iL;
      while (iL1 >= 0) {
        int jL1 = (iL1 == iL) ? jL : jL + 1;
        while (jL1 >= 1) {
          jL1--;
          IL[0] = iL;  IL[1] = jL;   
          JL[0] = iL1; JL[1] = jL1;
          for (unsigned  i = 0; i < n1D[iL]; i++) {
            unsigned i1 = static_cast < unsigned >(floor(i / n1D[IL[0] - JL[0]]));
            for (unsigned  j = 0; j < n1D[jL]; j++) {
              unsigned j1 = static_cast < unsigned >(floor(j / n1D[IL[1] - JL[1] ]));
              CH[IL[0]][IL[1]][i * n1D[IL[1]] + j] -= CH[ JL[0] ][ JL[1] ][i1 * n1D[ JL[1] ] + j1];
            }
          }
        }
        iL1--;
      }
    }
  }

  std::cout << std::endl;
  std::cout << "Hierarchical bases " << std::endl;
  for (unsigned iL = 0; iL < CH.size(); iL++) {
    for (unsigned jL = 0; iL * sparse + jL < CH[iL].size(); jL++) {
      IL[0] = iL;
      IL[1] = jL;
      PrintBaseCoefficients(CH[iL][jL], IL, n1D);
    }
  }

  //////////////////////////////////////////////////////////////////////

  //Reconstruct Histogram from Hierarchical Bases

  std::vector <double> cIr(static_cast< unsigned >(pow(n1D[L], N)));


  for (unsigned i = 0; i < n1D[L]; i++) {
    std::vector < unsigned > I(N);
    I[0] = i;
    for (unsigned j = 0; j < n1D[L]; j++) {

      I[1] = j;
      unsigned INDEX = I[0];
      for (unsigned k = 1; k < N; k++) {
        INDEX = INDEX * n1D[L] + I[k];
      }

      int iL = L;
      while (iL >= 0) {
        unsigned i1 = static_cast < unsigned >(floor(i / n1D[L - iL]));
        int jL = L;
        while (jL >= 0) {
          if (iL * sparse + jL < CH[iL].size()) {
            unsigned j1 = static_cast < unsigned >(floor(j / n1D[L - jL]));
            cIr[INDEX] += CH[iL][jL][i1 * n1D[jL] + j1];
          }
          jL--;
        }
        iL--;
      }
    }
  }

  std::cout << "Histogram reconstructed from Hierarchical Bases" << std::endl;
  IL.assign(N, L);
  PrintBaseCoefficients(cIr, IL, n1D);

  std::cout << "Difference between Histogram and Histogram reconstructed from Hierarchical Bases" << std::endl;
  for (unsigned i = 0; i < cIr.size(); i++) cIr[i] -= cI[i];
  PrintBaseCoefficients(cIr, IL, n1D);

  return 0;

} //end main


void PrintBaseCoefficients(const std::vector < double > & c, const std::vector <unsigned> &I, const std::vector < unsigned > &n) {

  for (unsigned k = 0; k < N; k++) {
    std::cout << "I[" << k << "] = " << I[k] << ", ";
  }
  std::cout << std::endl;

  unsigned N = I.size();
  unsigned Nm1 = N - 1u;
  unsigned Nm2 = N - 2u;
  unsigned nNm1 = n[ I[Nm1] ];
  for (unsigned i = 0; i < c.size(); i++) {
    std::cout << c[i] << " ";
    unsigned size = nNm1;
    unsigned ip1 = i + 1u;
    for (unsigned k = 0; k < N; k++) {
      if (ip1 % size == 0) std::cout << std::endl;
      if (k < Nm1) size *= n[ I[Nm2 - k]];
    }
  }
}

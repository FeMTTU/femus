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

bool sparse = false;

unsigned alpha = 4;
unsigned M = pow (10, alpha);   //number of samples
unsigned N = 2; //dimension of the parameter space (each of the M samples has N entries)
unsigned L = 4; //max refinement level
bool output = false; //for debugging
bool matlabView = true;

double xmin = - 5.5;   //-1.5 for uniform // -5.5 for Gaussian
double xmax = 5.5;     //1.5 for uniform // 5.5 for Gaussian

//FOR NORMAL DISTRIBUTION
boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)
boost::normal_distribution<> nd (0., 1.);
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

int main (int argc, char** argv) {

  //BEGIN construction of the sample set
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
      samples[m][n] = var;
      //samples[m][n] = varunif;
    }
  }
  //END

  
  std::vector< double > n1(L);
  n1[0] = 1u;
  for(unsigned i = 1; i < L; i++){
    n1[i] = n1[i-1] * 2u;    
  }

  //Build Histogram
   
  std::vector < std::vector <double> > cI(n1[L-1]);
  for (unsigned i = 0; i < n1[L-1]; i++) cI[i].assign (n1[L-1], 0.);

  double h = (xmax - xmin) / n1[L-1];
  double h2 = h * h;

  for (unsigned m = 0; m < M; m++) {
      
    std::vector < unsigned > i(N);
    for(unsigned k = 0; k < N; k++){
      double x = (samples[m][k] - xmin) / h;
      if (x < 0.) i[k] = 0;
      else if (x >= n1[L-1]) i[k] = n1[L-1] - 1;
      else i[k] = static_cast < unsigned > (floor (x));
    }
    cI[i[0]][i[1]] += 1. / (M * h2);
  }


//   std::cout << "Histogram "<<std::endl;
//   for (unsigned i = 0; i < n1[L-1]; i++) {
//     for (unsigned j = 0; j < n1[L-1]; j++) {
//       std::cout << cI[i][j] << " ";
//     }
//     std::cout << std::endl;
//   }


  //Build Hierarchical Bases

  std::vector < std::vector < std::vector <std::vector <double> > > > cH (L);
  for (unsigned iL = 0; iL < cH.size(); iL++) {
    cH[iL].resize (L);
    unsigned iDim = n1[iL];
    for (unsigned jL = 0; iL * sparse + jL < cH[iL].size(); jL++) {
      cH[iL][jL].resize (iDim);
      unsigned jDim = n1[jL];
      for (unsigned i = 0; i < cH[iL][jL].size(); i++) {
        cH[iL][jL][i].assign (jDim, 0.);
      }
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

    for (unsigned iL = 0; iL < cH.size(); iL++) {
      double hx = Hx / n1[iL];
      double x = X * n1[iL];
      for (unsigned jL = 0; iL * sparse + jL < cH[iL].size(); jL++) {
        double hy = Hy / n1[jL];
        double y = Y * n1[jL];
        unsigned i = static_cast < unsigned > (floor (x));
        unsigned j = static_cast < unsigned > (floor (y));

        cH[iL][jL][i][j] += 1. / (M * hx * hy);
      }
    }
  }



  for (unsigned iL = 0; iL < cH.size(); iL++) {
    for (unsigned jL = 0; iL * sparse + jL < cH[iL].size(); jL++) {
      //std::cout << iL << " " << jL << std::endl;
      int iL1 = iL;
      while (iL1 >= 0) {
        int jL1 = (iL1 == iL) ? jL : jL + 1;
        while (jL1 >= 1) {
          jL1--;
          //std::cout << "\t\t" << iL1 << " " << jL1 << std::endl;
          for (unsigned  i = 0; i < cH[iL][jL].size(); i++) {
            unsigned i1 = static_cast < unsigned > (floor (i / n1[iL-iL1]));
            for (unsigned  j = 0; j < cH[iL][jL][i].size(); j++) {
              unsigned j1 = static_cast < unsigned > (floor (j / n1[jL - jL1]));
              cH[iL][jL][i][j] -= cH[iL1][jL1][i1][j1];
            }
          }
        }
        iL1--;
        //std::cout<<std::endl;
      }
    }
  }



//   std::cout << std::endl;
//   std::cout << "Hierarchical bases "<<std::endl;
//   for (unsigned iL = 0; iL < cH.size(); iL++) {
//     for (unsigned jL = 0; jL < cH[iL].size(); jL++) {
//       std::cout << "iL = " << iL << " jL = "<< jL <<std::endl;
//       for (unsigned i = 0; i < cH[iL][iL].size(); i++) {
//         for (unsigned j = 0; j < cH[iL][jL][i].size(); j++) {
//           std::cout << cH[iL][jL][i][j] << " ";
//         }
//         std::cout << std::endl;
//       }
//       std::cout << std::endl;
//     }
//     std::cout << std::endl;
//   }



  //Reconstruct Histogram from Hierarchical Bases
  std::vector < std::vector <double> > cIr (n1[L-1]);
  for (unsigned i = 0; i < n1[L-1]; i++) cIr[i].assign (n1[L-1], 0.);

  for (unsigned i = 0; i < cIr.size(); i++) {
    for (unsigned j = 0; j < cIr[i].size(); j++) {
      int iL = L - 1;
      while (iL >= 0) {
        unsigned i1 = static_cast < unsigned > (floor (i / n1[L - 1 - iL]));
        int jL = L - 1;
        while (jL >= 0) {
          if (iL * sparse + jL < cH[iL].size()) {
            unsigned j1 = static_cast < unsigned > (floor (j / n1[L - 1 - jL]));
            cIr[i][j] += cH[iL][jL][i1][j1];
          }
          jL--;
        }
        iL--;
      }
    }
  }

  std::cout<<"Histogram reconstructed from Hierarchical Bases"<<std::endl;
  for (unsigned i = 0; i < n1[L-1]; i++) {
    for (unsigned j = 0; j < n1[L-1]; j++) {
      std::cout << cIr[i][j] << " ";
    }
    std::cout << std::endl;
  }
//   std::cout << std::endl;
//   std::cout << std::endl;


  std::cout<<"Difference between Histogram and Histogram reconstructed from Hierarchical Bases"<<std::endl;
  for (unsigned i = 0; i < n1[L-1]; i++) {
    for (unsigned j = 0; j < n1[L-1]; j++) {
      std::cout << (cI[i][j] - cIr[i][j]) << " ";
    }
    std::cout << std::endl;
  }
//   std::cout << std::endl;
//   std::cout << std::endl;

  return 0;

} //end main




#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "Marker.hpp"
#include "Line.hpp"

#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "../include/mpmFem.hpp"


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
//double NeumannFactor = 0.;

using namespace femus;

// OLD BEST RESULT WITH E = 4.2 * 1.e6, 5 levels, dt= 0.01, NR = 300, R0 = 1.5, factor = 1.3
// MOST BEST RESULT WITH E = 4.2 * 1.e6, 4 levels, dt= 0.01, NR = 300, R0 = 1.4, factor = 1.14,  beta = 0.3, Gamma = 0.5

void ComputeIndexSet (std::vector < std::vector <unsigned> > & Jp,
                      const unsigned & degree, const unsigned & dimension, const bool &output = false);

void GaussianElemination (std::vector<std::vector < double > > & A, std::vector < double> &x, const bool &output = false);

void GetChebyshev (std::vector<double> &T, const unsigned &n, const double &x, const bool &output = false);

int main (int argc, char** args) {


  bool output = true;

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, args, MPI_COMM_WORLD);

  std::vector < std::vector <unsigned> > alphaIndex;
  unsigned pOrder = 7;
  unsigned dim = 1;
  ComputeIndexSet (alphaIndex, pOrder, dim, output);

  unsigned nve = 5u;
  unsigned nel = nve - 1u;

  std::vector < double > Xv (nve);
  std::vector < double > hv (nve);


  Xv[0] = 0.;
  Xv[1] = 0.1;
  Xv[2] = 0.5;
  Xv[3] = 1.;
  Xv[4] = 1.3;

  hv[0] = 0.1;
  hv[1] = 0.4;
  hv[2] = 0.5;
  hv[3] = 0.5;
  hv[4] = 0.3;

  unsigned nDofs = 2;

  std::vector < std::vector < unsigned > > elemDofs (nel);
  for (unsigned iel = 0; iel < nel; iel++) {
    elemDofs[iel].resize (nDofs);
    elemDofs[iel][0] = iel;
    elemDofs[iel][1] = iel + 1;
  }

  unsigned Np = alphaIndex.size() + 5;

  std::vector<std::vector< double> >Xp (nel);

  for (unsigned iel = 0; iel < nel; iel++) {
    Xp[iel].resize (Np);
    for (unsigned p = 0; p < Np; p++) {
      Xp[iel][p] = Xv[iel] + 1.0 * rand() / RAND_MAX * (Xv[iel + 1] - Xv[iel]);
    }
  }

  std::vector < std::vector < std::vector< double> > > M (nve); // array of matrices

  for (unsigned i = 0; i < nve; i++) {
    M[i].resize (alphaIndex.size());
    for (unsigned k = 0; k < alphaIndex.size(); k++) {
      M[i][k].assign (alphaIndex.size() + 1, 0.);
    }
  }

  std::vector < double > T;
  GetChebyshev (T, pOrder, 0., output);

  for (unsigned i = 0; i < nve; i++) {
    for (unsigned j = 0; j < T.size(); j++) {
      M[i][j][alphaIndex.size()] = T[j];
    }
  }

  for (unsigned iel = 0; iel < nel; iel++) {
    for (unsigned p = 0; p < Np; p++) {

      for (unsigned idof = 0; idof < nDofs; idof++) {

        unsigned i = elemDofs[iel][idof];

        unsigned jdof = (idof == 0) ? 1 : 0;
        unsigned j = elemDofs[iel][jdof];

        GetChebyshev (T, pOrder, (Xv[i] - Xp[iel][p]) / hv[i]);

        double W = (1. - (Xv[i] - Xp[iel][p]) / (Xv[i] - Xv[j]));
        
        for (unsigned k = 0; k < alphaIndex.size(); k++) {
          for (unsigned l = 0; l < alphaIndex.size(); l++) {
            M[i][k][l] +=  W * T[k] * T[l];
          }
        }
        
      }
    }
  }

  std::vector < std::vector< double> > alpha (nve);

  for (unsigned i = 0; i < nve; i++) {
    alpha[i].resize (alphaIndex.size());
    GaussianElemination (M[i], alpha[i], output);
  }

  std::vector < double > Ur (nve, 0.);

  for (unsigned iel = 0; iel < nel; iel++) {
    for (unsigned p = 0; p < Np; p++) {

      for (unsigned idof = 0; idof < nDofs; idof++) {

        unsigned i = elemDofs[iel][idof];
        
        unsigned jdof = (idof == 0) ? 1 : 0;
        unsigned j = elemDofs[iel][jdof];
        
        GetChebyshev (T, pOrder, (Xv[i] - Xp[iel][p]) / hv[i]);

        double sumAlphaT = 0.;
        for (unsigned k = 0; k < alphaIndex.size(); k++) {
          sumAlphaT += alpha[i][k] * T[k];
        }
        
        double W = (1. - (Xv[i] - Xp[iel][p]) / (Xv[i] - Xv[j]));
        
        Ur[i] += W * sumAlphaT  * pow (Xp[iel][p], pOrder) ;
      }
    }
  }


  std::cout << std::endl;
  for (unsigned i = 0; i < nve; i++) {
    std::cout << pow (Xv[i], pOrder) << " " << Ur[i] << std::endl;
  }
  std::cout << std::endl;

}


void GaussianElemination (std::vector<std::vector < double > > & A, std::vector < double> &x, const bool &output) {

  unsigned n = A.size();

  if (output) {
    std::cout << "Before LU\n";
    for (unsigned i = 0; i < n; i++) {
      for (unsigned j = 0; j < n; j++) {
        std::cout << A[i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  for (unsigned i = 0; i < n - 1; i++) {
    unsigned p = i;
    while (A[p][i] == 0) {
      p++;
      if (p == n) {
        std::cout << "The Matrix A is singular\n";
        exit (0);
      }
    }
    if (p != i) {
      for (unsigned j = 0; j < n + 1; j++) {
        double tmp;
        tmp = A[i][j];
        A[i][j] = A[p][j];
        A[p][j] = tmp;
      }
    }
    for (unsigned j = i + 1; j < n; j++) {
      double mji = A[j][i] / A[i][i];
      for (unsigned k = i; k < n + 1; k++) {
        A[j][k] -= mji * A[i][k];
      }
    }
  }
  if (A[n - 1][n - 1] == 0) {
    std::cout << "The Matrix A is singular\n";
    exit (0);
  }
  else {
    x[n - 1] = A[n - 1][n] / A[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--) {
      x[i] = A[i][n];
      for (unsigned j = i + 1; j < n; j++) {
        x[i] -= A[i][j] * x[j];
      }
      x[i] /= A[i][i];
    }
  }

  if (output) {
    std::cout << "After LU\n";
    for (unsigned i = 0; i < n; i++) {
      for (unsigned j = 0; j < n; j++) {
        std::cout << A[i][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  return;
}

void ComputeIndexSet (std::vector < std::vector <unsigned> > & Jp,
                      const unsigned & degree, const unsigned & dimension, const bool &output) { //p is max poly degree


  unsigned dimJp = static_cast <unsigned> (boost::math::binomial_coefficient<double> (dimension + degree, degree));

  Jp.resize (dimJp);
  for (unsigned i = 0; i < dimJp; i++) {
    Jp[i].resize (dimension);
  }

  unsigned index = 0;
  unsigned counters[dimension + 1];
  memset (counters, 0, sizeof (counters));

  while (!counters[dimension]) {

    unsigned entrySum = 0;
    for (unsigned j = 0; j < dimension; j++) {
      entrySum += counters[j];
    }

    if (entrySum <= degree) {
      for (unsigned j = 0; j < dimension; j++) {
        Jp[index][j] = counters[dimension - 1 - j];
        if (output) {
          std::cout << "alpha[" << index << "][" << j << "]= " << Jp[index][j] ;
        }
      }
      if (output) {
        std::cout << std::endl;
      }
      index++;
    }
    unsigned i;
    for (i = 0; counters[i] == degree; i++) {   // inner loops that are at maxval restart at zero
      counters[i] = 0;
    }
    ++counters[i];  // the innermost loop that isn't yet at maxval, advances by 1
  }
  if (output) {
    std::cout << std::endl;
  }
}


void GetChebyshev (std::vector<double> &T, const unsigned &n, const double &x, const bool &output) {
  T.resize (n + 1);
  T[0] = 1;
  T[1] = x;
  for (unsigned i = 2; i < n + 1; i++) {
    T[i] = 2 * x * T[ i - 1] - T[ i - 2];
  }
  if (output) {
    std::cout << "Chebyshev Polynomilas at x = " << x << std::endl;
    for (unsigned i = 0; i < n + 1; i++) {
      std::cout << "T" << i << " [x] = " << T[i] << std::endl;
    }
    std::cout << std::endl;

  }

}


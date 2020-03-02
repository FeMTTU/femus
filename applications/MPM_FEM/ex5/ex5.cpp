
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
                      const unsigned & degree, const unsigned & dimension, bool output = false);

void GaussianElemination (std::vector<std::vector < double > > & A, std::vector < double> &x);

void GaussianEleminationWithPivoting (std::vector<std::vector<double>> &A, std::vector<double> &x);

int main (int argc, char** args) {


  // init Petsc-MPI communicator
  FemusInit mpinit (argc, args, MPI_COMM_WORLD);

  std::vector < std::vector <unsigned> > alphaIndex;
  unsigned pOrder = 6;
  unsigned dim = 1;
  ComputeIndexSet (alphaIndex, pOrder, dim, true);

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


  unsigned Np = alphaIndex.size();


  std::vector<std::vector< double> >Xp (nel);

  for (unsigned iel = 0; iel < nel; iel++) {
    Xp[iel].resize (Np);
    for (unsigned p = 0; p < Np; p++) {
      Xp[iel][p] = Xv[iel] + 1.0 * rand() / RAND_MAX * (Xv[iel + 1] - Xv[iel]);
      std::cout << Xp[iel][p] << " ";
    }
    std::cout << "\n" << std::endl;
  }


  std::vector < std::vector < std::vector< double> > > M (nve); // array of matrices


  for (unsigned i = 0; i < nve; i++) {
    M[i].resize (alphaIndex.size());
    for (unsigned k = 0; k < alphaIndex.size(); k++) {
      M[i][k].assign (alphaIndex.size() + 1, 0.);
    }
  }

  for (unsigned i = 0; i < nve; i++) {
    M[i][0][alphaIndex.size()] = 1.;
  }


  for (unsigned iel = 0; iel < nel; iel++) {
    for (unsigned p = 0; p < Np; p++) {
      for (unsigned k = 0; k < alphaIndex.size(); k++) {
        for (unsigned l = 0; l < alphaIndex.size(); l++) {
          M[iel][k][l] += (1. - (Xv[iel] - Xp[iel][p]) / (Xv[iel] - Xv[iel + 1])) * pow ( (Xv[iel] - Xp[iel][p]) / hv[iel], alphaIndex[k][0] + alphaIndex[l][0]);
          M[iel + 1][k][l] += (1. - (Xv[iel + 1] - Xp[iel][p]) / (Xv[iel + 1] - Xv[iel])) * pow ( (Xv[iel + 1] - Xp[iel][p]) / hv[iel + 1], alphaIndex[k][0] + alphaIndex[l][0]);
        }
      }
    }
  }


  std::vector < std::vector< double> > alpha (nve);

  for (unsigned i = 0; i < nve; i++) {
    alpha[i].resize (alphaIndex.size());
    GaussianEleminationWithPivoting (M[i], alpha[i]);
  }

  std::vector < double > Ur (nve, 0.);


  for (unsigned iel = 0; iel < nel; iel++) {

    for (unsigned p = 0; p < Np; p++) {

      std::vector < double > h (alphaIndex.size());

      double det = (Xv[iel] - Xp[iel][p]) / hv[iel];
      for (unsigned k = 0; k < alphaIndex.size(); k++) {
        h[k] = pow (det, alphaIndex[k][0]);
      }

      det = 0.;
      for (unsigned k = 0; k < alphaIndex.size(); k++) {
        det += alpha[iel][k] * h[k];
      }

      Ur[iel]     += (1. - (Xv[iel] - Xp[iel][p]) / (Xv[iel] - Xv[iel + 1])) * det  * pow (Xp[iel][p], pOrder) ;

      det = (Xv[iel + 1] - Xp[iel][p]) / hv[iel + 1];
      for (unsigned k = 0; k < alphaIndex.size(); k++) {
        h[k] = pow (det, alphaIndex[k][0]);
      }

      det = 0.;
      for (unsigned k = 0; k < alphaIndex.size(); k++) {
        det += alpha[iel + 1][k] * h[k];
      }

      Ur[iel + 1] += (1. - (Xv[iel + 1] - Xp[iel][p]) / (Xv[iel + 1] - Xv[iel])) * det * pow (Xp[iel][p], pOrder);
    }

  }


  std::cout << std::endl;

  for (unsigned i = 0; i < nve; i++) {
    std::cout << pow (Xv[i], pOrder) << " " << Ur[i] << std::endl;

  }
  std::cout << std::endl;

}


void GaussianElemination (std::vector<std::vector < double > > & A, std::vector < double> &x) {

  unsigned n = A.size();

//   for(unsigned i = 0; i < n; i++){
//     for(unsigned j =0; j< n; j++){
//       std::cout << A[i][j] <<" ";
//     }
//     std::cout<<std::endl;
//   }
//   std::cout<<std::endl;
//
  for (unsigned i = 0; i < n - 1; i++) {
    unsigned p = i;
    while (A[p][i] < pow (10, -8)) {
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

  for (unsigned i = 0; i < n; i++) {
    for (unsigned j = 0; j <= n; j++) {
      std::cout << A[i][j] << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  return;
}

///----------------------



void GaussianEleminationWithPivoting (std::vector<std::vector<double>> &A, std::vector<double> &x) {
// A is nx(n+1) augmented matrix
  int n  = A.size();

  std::cout << "Before pivoting\n A=\n" ;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= n; j++) {
      std::cout << A[i][j] << " " ;
    }
    std::cout << "\n" ;
  }
  std::cout << "\n";

  const double tol = 1e-6;
  double max_row = 0.0, abs;
  for (int i = 0; i < n; i++) {
    for (int k = i + 1; k < n; k++) {
      max_row = fabs (A[i][i]);
      if ( (abs = fabs (A[k][i])) > max_row) {
        max_row = abs;
        for (int j = 0; j < n + 1; j++) {
          double temp = A[i][j];
          A[i][j] = A[k][j];
          A[k][j] = temp;
        }
      }
    }
    if (max_row < tol) {
      std::cout << "Degenerate matrix.";
      exit (0);
    }

    for (int j = i + 1; j < n; j++) {
      double mji = A[j][i] / A[i][i];
      for (int k = i; k < n + 1; k++) {
        A[j][k] -= mji * A[i][k];
      }
    }

  }
  if (fabs (A[n - 1][n - 1]) < tol) {
    std::cout << "Zero row! Matrix is singular \n";
    exit (0);
  }

  x[n - 1] = A[n - 1][n] / A[n - 1][n - 1];
  for (int i = n - 2; i >= 0; i--) {
    x[i] = A[i][n];
    for (int j = i + 1; j < n; j++) {
      x[i] -= A[i][j] * x[j];
    }
    x[i] /= A[i][i];
  }

  std::cout << "After pivoting \nA=\n" ;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= n; j++) {
      std::cout << A[i][j] << " " ;
    }
    std::cout << "\n" ;
  }
  std::cout << "\n\n";

  std::cout << "Solution is: \n";
  for (int i = 0; i < n; i++) {
    std::cout << x[i] << "\n";
  }
  return;

}


///---------------------







void ComputeIndexSet (std::vector < std::vector <unsigned> > & Jp,
                      const unsigned & degree, const unsigned & dimension, bool output) { //p is max poly degree


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


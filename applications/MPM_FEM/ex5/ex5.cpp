
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

int main (int argc, char** args) {


  std::vector < std::vector <unsigned> > Jp;
  unsigned polynomial_order = 4;
  unsigned dim = 2;
  ComputeIndexSet (Jp, polynomial_order, dim, true);

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, args, MPI_COMM_WORLD);

  unsigned N = 5u;
  std::vector < double > Xv (N);
  Xv[0] = 0.;
  Xv[1] = 0.1;
  Xv[2] = 0.5;
  Xv[3] = 1.;
  Xv[4] = 1.3;

  unsigned Np = 10;
  std::vector<std::vector< double> >Xp (N - 1);

  for (unsigned iel = 0; iel < N - 1; iel++) {
    Xp[iel].resize (Np);
    for (unsigned p = 0; p < Np; p++) {
      Xp[iel][p] = Xv[iel] + 1.0 * rand() / RAND_MAX * (Xv[iel + 1] - Xv[iel]);
      std::cout << Xp[iel][p] << " ";
    }
    std::cout << std::endl;
  }

  std::vector < std::vector < std::vector< double> > >M (N);

  unsigned Nr = 7;

  for (unsigned i = 0; i < N; i++) {
    M[i].resize (Nr);
    for (unsigned k = 0; k < Nr; k++) {
      M[i][k].assign (Nr + 1, 0.);
    }
  }

  for (unsigned i = 0; i < N; i++) {
    M[i][0][Nr] = 1.;
  }


  for (unsigned iel = 0; iel < N - 1; iel++) {
    for (unsigned p = 0; p < Np; p++) {
      for (unsigned k = 0; k < Nr; k++) {
        for (unsigned l = 0; l < Nr; l++) {
          M[iel][k][l] += (1. - (Xv[iel] - Xp[iel][p]) / (Xv[iel] - Xv[iel + 1])) * pow ( (Xv[iel] - Xp[iel][p]), k + l);
          M[iel + 1][k][l] += (1. - (Xv[iel + 1] - Xp[iel][p]) / (Xv[iel + 1] - Xv[iel])) * pow ( (Xv[iel + 1] - Xp[iel][p]), k + l);
        }
      }
    }
  }


  std::vector < std::vector< double> > alpha (N);

  for (unsigned i = 0; i < N; i++) {
    alpha[i].resize (Nr);
    GaussianElemination (M[i], alpha[i]);
  }

  std::vector < double > Ur (N, 0.);

  double ptest = 6;
  for (unsigned iel = 0; iel < N - 1; iel++) {
    for (unsigned p = 0; p < Np; p++) {

      std::vector < double > h (Nr);
      h[0] = 1;
      double det = (Xv[iel] - Xp[iel][p]);
      for (unsigned k = 1; k < Nr; k++) {
        h[k] = pow (det, k);
      }
      det = 0.;

      for (unsigned k = 0; k < Nr; k++) {
        det += alpha[iel][k] * h[k];
      }

      Ur[iel]     += (1. - (Xv[iel] - Xp[iel][p]) / (Xv[iel] - Xv[iel + 1])) * det  * pow (Xp[iel][p], ptest) ;

      h[0] = 1;
      det = (Xv[iel + 1] - Xp[iel][p]);
      for (unsigned k = 1; k < Nr; k++) {
        h[k] = pow (det, k);
      }
      det = 0.;

      for (unsigned k = 0; k < Nr; k++) {
        det += alpha[iel + 1][k] * h[k];
      }


      Ur[iel + 1] += (1. - (Xv[iel + 1] - Xp[iel][p]) / (Xv[iel + 1] - Xv[iel])) * det * pow (Xp[iel][p], ptest);
    }
  }


  std::cout << std::endl;
  for (unsigned i = 0; i < N; i++) {
    std::cout << pow (Xv[i], ptest) << " " << Ur[i] << std::endl;
  }
  std::cout << std::endl;

}


void GaussianElemination (std::vector<std::vector < double > > & A, std::vector < double> &x) {

  unsigned n = A.size();
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
  return;
}

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
          std::cout << " Jp[" << index << "][" << j << "]= " << Jp[index][j] ;
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
}

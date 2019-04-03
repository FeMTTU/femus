
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



int main (int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit (argc, args, MPI_COMM_WORLD);

  unsigned N = 5u;
  std::vector < double > Xv (N);
  Xv[0] = 0.;
  Xv[1] = 0.1;
  Xv[2] = 0.5;
  Xv[3] = 1.;
  Xv[4] = 1.3;

  unsigned Np = 3;
  std::vector<std::vector< double> >Xp (N - 1);

  for (unsigned iel = 0; iel < N - 1; iel++) {
    Xp[iel].resize (Np);
    for (unsigned p = 0; p < Np; p++) {
      Xp[iel][p] = Xv[iel] + (p + 1) * (Xv[iel + 1] - Xv[iel]) / (Np + 1); // particle points
      std::cout << Xp[iel][p] << " ";
    }
    std::cout << std::endl;
  }

  std::vector < std::vector < std::vector< double> > >M (N); // array of matrices

  for (unsigned i = 0; i < N; i++) {
    M[i].resize (3);
    for (unsigned k = 0; k < 3; k++) {
      M[i][k].assign (3, 0.); //assigns new contents to the vector, replacing its current contents, and modifying its size accordingly.
    }
  }


  for (unsigned iel = 0; iel < N - 1; iel++) {
    for (unsigned p = 0; p < Np; p++) {
      for (unsigned k = 0; k < 3; k++) {
        for (unsigned l = 0; l < 3; l++) {
          M[iel][k][l] += (1. - (Xv[iel] - Xp[iel][p]) / (Xv[iel] - Xv[iel + 1])) * pow ( (Xv[iel] - Xp[iel][p]), k + l);
          M[iel + 1][k][l] += (1. - (Xv[iel + 1] - Xp[iel][p]) / (Xv[iel + 1] - Xv[iel])) * pow ( (Xv[iel + 1] - Xp[iel][p]), k + l);
        }
      }
    }
  }
  
  std::vector < std::vector < std::vector< double> > >MI (N);
  
  for (unsigned i = 0; i < N; i++) {
    MI[i].resize (3);
    for (unsigned k = 0; k < 3; k++) {
      MI[i][k].assign (3, 0.);
    }
  }
  
  for (unsigned i = 0; i < N; i++) {
    double det = (M[i][0][0] * (M[i][1][1] * M[i][2][2] - M[i][1][2] * M[i][2][1]) +
    M[i][0][1] * (M[i][1][2] * M[i][2][0] - M[i][1][0] * M[i][2][2]) +
    M[i][0][2] * (M[i][1][0] * M[i][2][1] - M[i][1][1] * M[i][2][0]));
  
    MI[i][0][0] = (-M[i][1][2] * M[i][2][1] + M[i][1][1] * M[i][2][2]) / det;
    MI[i][0][1] = (M[i][0][2] * M[i][2][1] - M[i][0][1] * M[i][2][2]) / det;
    MI[i][0][2] = (-M[i][0][2] * M[i][1][1] + M[i][0][1] * M[i][1][2]) / det;
    MI[i][1][0] = (M[i][1][2] * M[i][2][0] - M[i][1][0] * M[i][2][2]) / det;
    MI[i][1][1] = (-M[i][0][2] * M[i][2][0] + M[i][0][0] * M[i][2][2]) / det;
    MI[i][1][2] = (M[i][0][2] * M[i][1][0] - M[i][0][0] * M[i][1][2]) / det;
    MI[i][2][0] = (-M[i][1][1] * M[i][2][0] + M[i][1][0] * M[i][2][1]) / det;
    MI[i][2][1] = (M[i][0][1] * M[i][2][0] - M[i][0][0] * M[i][2][1]) / det;
    MI[i][2][2] = (-M[i][0][1] * M[i][1][0] + M[i][0][0] * M[i][1][1]) / det;
  }
  
  
  std::vector < std::vector < double > > I(3);
  for (unsigned k = 0; k < 3; k++) {
    I[k].resize (3);
  }
  
  for (unsigned i = 0; i < N; i++) {
    
    for (unsigned j = 0; j < 3; j++) {
      for (unsigned k = 0; k < 3; k++) {
        I[j][k] = 0.;
        for (unsigned l = 0; l < 3; l++) {
          I[j][k] += M[i][j][l] * MI[i][l][k];
        }
      }
    }
    for (unsigned j = 0; j < 3; j++) {
      for (unsigned k = 0; k < 3; k++) {
        std::cout << I[j][k]<<" ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  
  std::vector < double > h0(3,0.);
  h0[0] = 1.;
  
  
  std::vector < double > Ur(N,0.);
  
  double ptest = 3;
  for (unsigned iel = 0; iel < N - 1; iel++) {
    for (unsigned p = 0; p < Np; p++) {
      
      std::vector < double > h(3);
      h[0] = 1;
      h[1] = (Xv[iel] - Xp[iel][p]);
      h[2] = h[1] * h[1];
      double det = 0.;
      for (unsigned j = 0; j < 3; j++) {
        for (unsigned k = 0; k < 3; k++) {
          det += h0[j] * MI[iel][j][k] * h[k];
        }
      }
      
      Ur[iel]     += (1. - (Xv[iel] - Xp[iel][p]) / (Xv[iel] - Xv[iel + 1])) * det  * pow(Xp[iel][p], ptest) ;
      
      h[0] = 1;
      h[1] = (Xv[iel+1] - Xp[iel][p]);
      h[2] = h[1] * h[1];
      det = 0.;
      for (unsigned j = 0; j < 3; j++) {
        for (unsigned k = 0; k < 3; k++) {
          det += h0[j] * MI[iel+1][j][k] * h[k];
        }
      }
  }
  
  for (unsigned i = 0; i < N; i++) {
    std::cout << pow(Xv[i],ptest) << " " << Ur[i] << std::endl;
  }
  std::cout << std::endl;

}

          
      Ur[iel + 1] += (1. - (Xv[iel + 1] - Xp[iel][p]) / (Xv[iel + 1] - Xv[iel])) * det * pow(Xp[iel][p], ptest);
    }

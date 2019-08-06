#include "FemusInit.hpp"
#include "./include/gmpm.hpp"

using namespace femus;

int main (int argc, char** args) {

  FemusInit (argc, args, MPI_COMM_WORLD);

  //std::vector < double > x {0., 0.43, 0.81, 1.25, 1.6};
  unsigned nve = 5;
  double L = 2.;
  std::vector < double > x (nve);

  double H = L / (nve - 1.);

  x[0] = 0.;
  for (unsigned i = 1; i < nve; i++) {
    x[i] = x[i - 1] + H;
  }

  std::vector < std::vector < double > > P (nve);
  std::vector < std::vector < double > > M (nve);
  std::vector < std::vector < double > > K (nve);
  std::vector < std::vector < double > > MP (nve);
  std::vector < std::vector < double > > Kt (nve);


  for (unsigned i = 0; i < nve; i++) {
    P[i].assign (nve, 0.);
    M[i].assign (nve, 0.);
    K[i].assign (nve + 1, 0.);
    MP[i].assign (nve, 0.);
    Kt[i].assign (nve + 1, 0.);
  }

  for (unsigned iel = 0; iel < nve - 1; iel++) {

    unsigned i0 = iel;
    unsigned i1 = iel + 1;
    double x0 = x[i0];
    double x1 = x[i1];
    double h = x1 - x0;

    M[i0][i0] += 1. / 6. * h * 2;
    M[i0][i1] += 1. / 6. * h;
    M[i1][i0] += 1. / 6. * h;
    M[i1][i1] += 1. / 6. * h * 2;

    K[i0][i0] +=  1. / h;
    K[i0][i1] += -1. / h;
    K[i1][i0] += -1. / h;
    K[i1][i1] +=  1. / h;

  }
  P[0][0] = -1. / (x[1] - x[0]);
  P[0][1] =  1. / (x[1] - x[0]);
  for (unsigned i = 1; i < nve - 1; i++) {
    double h = x[i + 1] - x[i - 1];
    P[i][i - 1] = -1. / h;
    P[i][i + 1] =  1. / h;
  }
  P[nve - 1][nve - 2] = -1. / (x[nve - 1] - x[nve - 2]);
  P[nve - 1][nve - 1] =  1. / (x[nve - 1] - x[nve - 2]);

  for (unsigned i = 0; i < nve; i++) {
    for (unsigned j = 0; j < nve; j++) {
      for (unsigned k = 0; k < nve; k++) {
        MP[i][j] += M[i][k] * P[k][j];
      }
    }
  }

  for (unsigned i = 0; i < nve; i++) {
    for (unsigned j = 0; j < nve; j++) {
      for (unsigned k = 0; k < nve; k++) {
        Kt[i][j] += P[k][i] * MP[k][j];
      }
    }
  }

  /*//point load
  for(unsigned j = 1; j< nve + 1; j++){
    K[0][j] = 0.;
    Kt[0][j] = 0.;
  }
  K[nve - 1][nve] = 1.;
  Kt[nve - 1][nve] = 1.;
  */


  //parabola
  for (unsigned j = 0; j < nve + 1; j++) {
    K[0][j] = 0.;
    K[nve - 1][j] = 0.;

    Kt[0][j] = 0.;
    Kt[nve - 1][j] = 0.;
  }
  K[0][0] = 1.;
  K[nve - 1][nve - 1] = 1.;
  Kt[0][0] = 1.;
  Kt[nve - 1][nve - 1] = 1.;

  for (unsigned i = 1; i < nve - 1; i++) {
    K[i][nve] = 2. * H;
    Kt[i][nve] = 2. * H;
  }

  std::vector < double> u (nve);
  GaussianEleminationWithPivoting (K, u, true);
  std::cout << u[ (nve - 1) / 2] << std::endl;

  GaussianEleminationWithPivoting (Kt, u, true);
  std::cout << u[ (nve - 1) / 2] << std::endl;

//   for (unsigned i = 0; i < nve; i++) {
//     for (unsigned j = 0; j < nve + 1; j++) {
//       std::cout << Kt[i][j] <<" ";
//     }
//     std::cout << std::endl;
//   }
//   std::cout << std::endl;
//


}





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

#include<PolynomialBases.hpp>

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

void GaussianEleminationWithPivoting (std::vector<std::vector < double > > & A, std::vector < double> &x, const bool &output = false);
void LUwithPivoting (std::vector<std::vector < double > > & A, std::vector < unsigned > &idx, const bool &output = false);
void LUbackward (const std::vector<std::vector < double > > & A, const std::vector < unsigned > &idx,
                 std::vector < double > &b, std::vector < double > &x, const bool &output = false);
void LUsolve (std::vector<std::vector < double > > & A, std::vector < unsigned > &idx,
              std::vector < double > &b, std::vector < double > &x, const bool &output = false);

void GetChebyshev (std::vector<double> &T, std::vector<double> &dT, const unsigned &n, const double &x, const bool &output = false);
void GetPolynomial (std::vector<double> &P, std::vector<double> &dP, const unsigned &n, const double &x, const bool &output = false);
void GetMultiIndex (std::vector <unsigned> &idx, const unsigned &dim, const unsigned& n, const unsigned &i);
void SetElementDofs (std::vector <unsigned> &elementDofs, const std::vector < unsigned > & idx, const unsigned & nve1d);

void PrintSolution (const std::vector <double> &U, const std::vector <double> &Ur, const double *x, const unsigned &dim, const unsigned &n);
void PrintGnuplotScript (const double & xmin, const double & xmax, const double & ymin, const double & ymax, const unsigned &nve, const bool&printDerivative = false);

class GMPM {
  public:
    GMPM (const unsigned &dim, const unsigned &Nv) {
      _dim = dim;
      _xp.resize (dim);
      _node.reserve (Nv);
      _s.reserve (Nv);
    };

    unsigned _dim;
    std::vector < unsigned > _node;
    std::vector < double > _xp;
    std::vector < std::vector < double > > _s;
};

int main (int argc, char** args) {

  //FemusInit mpinit (argc, args, MPI_COMM_WORLD);

  bool nonLocal = true;
  
  bool output = true;

  std::vector < std::vector <unsigned> > aIdx;
  unsigned pOrder = 2;
  unsigned dim = 1;
  ComputeIndexSet (aIdx, pOrder, dim, output);

  unsigned nve1d = 9u;
  unsigned nve = static_cast< unsigned > (pow (nve1d, dim));
// unsigned nel1d = nve1d - 1u;
// unsigned nel = static_cast< unsigned > (pow (nel1d, dim));

  double Xv[9] = {0., 0.186, 0.409, 0.567, 0.8235, 1.03, 1.2222, 1.4345, 1.6};

  double scale = 0.25;

  std::vector < double > sMax (nve1d);
  std::vector < double > sMin (nve1d);

  for (int i = 0; i < nve1d; i++) {
    int deltai = (i % 2 == 0) ? pOrder + nonLocal : pOrder - (!nonLocal);

    int im = (i - deltai >= 0) ?  i - deltai : 0;
    int ip = (i + deltai < nve) ? i + deltai : nve - 1;

    sMin[i] = Xv[im] - Xv[i];
    sMax[i] = Xv[ip] - Xv[i];

    if(nonLocal){
      if( i - deltai < 0 ) sMin[i] += Xv[0] - Xv[1];
      if( i + deltai >= nve1d) sMax[i] += Xv[nve1d - 1] - Xv[nve1d - 2];
    }
   
  }

  std::vector< std::vector < std::vector< double> > > XVprint (1);
  XVprint[0].resize (nve);
  for (unsigned i = 0; i < XVprint[0].size(); i++) {
    XVprint[0][i].assign (dim, 0.);
  }

  for (unsigned i = 0; i < XVprint[0].size(); i++) {
    for (unsigned d = 0; d < dim; d++) {
      XVprint[0][i][d] = Xv[i];
    }
  }

  //PrintLine ("./output/", XVprint, false, 1);

  unsigned Np = 500;

  unsigned maxNumberOfNodes = (2u * (pOrder + 1u) < nve) ? 2u * (pOrder + 1u) : nve;

  std::vector < GMPM *> gmpm (Np);
  for (unsigned p = 0; p < gmpm.size(); p++) {
    gmpm[p] = new GMPM (dim, maxNumberOfNodes);
  }

  std::vector <double> distance (dim);

  maxNumberOfNodes = 0;

  for (unsigned p = 0; p < gmpm.size(); p++) {
    for (unsigned d = 0; d < dim; d++) {
      gmpm[p]->_xp[d] = Xv[0] + (Xv[nve - 1] - Xv[0]) / (Np - 1) * (p);
    }
    for (unsigned j = 0; j < nve; j++) {
      bool particleIsWhithin = false;
      for (unsigned d = 0; d < dim; d++) {
        distance[d] = gmpm[p]->_xp[d] - Xv[j];
        if ( (distance[d] >= 0 && distance[d] <= sMax[j]) || (distance[d] <= 0 && distance[d] >= sMin[j])) {
          particleIsWhithin = true;
        }
        //if (fabs (distance[d]) <= distanceMax[j] )  particleIsWhithin = true;
      }
      if (particleIsWhithin) {
        unsigned size = gmpm[p]->_node.size();
        gmpm[p]->_node.resize (size + 1);
        gmpm[p]->_s.resize (size + 1);

        gmpm[p]->_node[size] = j;
        gmpm[p]->_s[size].resize (dim);
        for (unsigned d = 0; d < dim; d++) {
          gmpm[p]->_s[size][d] = distance[d];
        }
      }
    }
    if (maxNumberOfNodes <  gmpm[p]->_node.size()) {
      maxNumberOfNodes = gmpm[p]->_node.size();
    }
  }

  std::cout << "MaxNumberOfNOdes = " << maxNumberOfNodes << std::endl;

  std::vector< std::vector < std::vector< double> > > XPprint (1);
  XPprint[0].resize (gmpm.size());
  for (unsigned i = 0; i < gmpm.size(); i++) {
    XPprint[0][i].resize (dim);
  }
  for (unsigned p = 0; p < gmpm.size(); p++) {
    for (unsigned d = 0; d < dim; d++) {
      XPprint[0][p][d] = gmpm[p]->_xp[d];
    }
  }
  //PrintLine ("./output/", XPprint, false, 0);



  std::vector < unsigned > pOrderTest (dim);
  for (unsigned d = 0; d < dim; d++) {
    pOrderTest[d] = pOrder / dim;
  }
  pOrderTest[0] += pOrder % dim;
  std::cout << "Testing Polynomial \nPn = ";
  for (unsigned d = 0 ; d < dim; d++) {
    std::cout << "x" << d << "^" << pOrderTest[d] << " * ";
  }
  std::cout << "\b\b  " << std::endl;
  std::vector < double > Ur (Np, 0.);
  std::vector < double > dUr (Np, 0.);
  std::vector < double > Ue (Np, 0.);

  std::ofstream fout;
  std::ofstream fouti;
  std::vector < std::ofstream > dfout (dim);
  std::vector < std::ofstream > dfouti (dim);
  fout.open ("./output/phiSum.txt");
  for (unsigned d = 0; d < dim; d++) {
    std::ostringstream stream;
    stream << "./output/dphi" << "dx" << d << "Sum.txt";
    dfout[d].open (stream.str().c_str());
  }


  for (unsigned i = 0; i < nve; i++) {
    std::ostringstream stream;
    stream << "./output/phi" << i << ".txt";
    fouti.open (stream.str().c_str());
    fouti.close();

    for (unsigned d = 0; d < dim; d++) {
      std::ostringstream stream;
      stream << "./output/dphi" << i << "dx" << d << ".txt";
      dfouti[d].open (stream.str().c_str());
      dfouti[d].close();
    }
  }
  double ymin = 0.;
  double ymax = 0.;

  std::vector < double > dymin (dim, 0.);
  std::vector < double > dymax (dim, 0.);

  std::vector< double> alpha (aIdx.size());
  std::vector< std::vector< double> > dalpha (dim);
  for (unsigned d = 0; d < dim; d++) {
    dalpha[d].resize (aIdx.size());
  }

  std::vector< double> b (aIdx.size());
  std::vector< unsigned> pivotIndex (aIdx.size());

  std::vector < std::vector< double> > Mp (aIdx.size());
  for (unsigned k = 0; k < aIdx.size(); k++) {
    Mp[k].resize (aIdx.size() + 1);
    //Mp[k].resize (aIdx.size());
  }

  std::vector < std::vector < std::vector< double> > > dMp (dim);
  for (unsigned d = 0; d < dim; d++) {
    dMp[d].resize (aIdx.size());
    for (unsigned k = 0; k < aIdx.size(); k++) {
      dMp[d][k].resize (aIdx.size() + 1);
    }
  }


  // array of matrices
  std::vector < std::vector < std::vector < double > > > T;
  std::vector < std::vector < std::vector < double > > > dT;
  T.resize (maxNumberOfNodes);
  dT.resize (maxNumberOfNodes);
  for (unsigned i = 0; i < maxNumberOfNodes; i++) {
    T[i].resize (dim);
    dT[i].resize (dim);
  }

  std::vector < double > weight (maxNumberOfNodes);
  std::vector < std::vector < double > > dweight (maxNumberOfNodes);
  for (unsigned i = 0; i < maxNumberOfNodes; i++) {
    dweight[i].resize (dim);
  }

  std::vector < double > T0;
  std::vector < double > dT0;
  //GetChebyshev (T0, dT0, pOrder, 0., true);
  GetPolynomial (T0, dT0, pOrder, 0., false);

  for (unsigned p = 0; p < Np; p++) { // particle loop

    for (unsigned k = 0; k < aIdx.size(); k++) {
      Mp[k].assign (aIdx.size() + 1, 0.);
      //Mp[k].assign (aIdx.size(), 0.);
    }

    for (unsigned d = 0; d < dim; d++) {
      for (unsigned k = 0; k < aIdx.size(); k++) {
        dMp[d][k].assign (aIdx.size() + 1, 0.);
      }
    }

    for (unsigned i = 0; i <  gmpm[p]->_node.size(); i++) {

      unsigned inode = gmpm[p]->_node[i];
      
//       double dist2 = 0.;
//       for (unsigned d = 0 ; d < dim; d++) { // multidimensional loop
//         dist2 += gmpm[p]->_s[i][d] * gmpm[p]->_s[i][d];
//       }
//       dist2 /= distanceMax[inode] * distanceMax[inode];
//       if (dist2 > 1.) {
//         weight[i] = 0.;
//         dweight[i].assign (dim, 0.);
//       }
//       else {
//         weight[i] = pow (1. - dist2, 4);
//         for (unsigned d = 0 ; d < dim; d++) {
//           dweight[i][d] = 4. * pow (1. - dist2, 3) * (-2.) * gmpm[p]->_s[i][d] / (distanceMax[inode] * distanceMax[inode]);
//         }
//       }
      if(nonLocal){
        double s = gmpm[p]->_s[i][0];
        weight[i] = pow ( (1. - s / sMax[inode]) * (1. - s / sMin[inode]), 4.);
        for (unsigned d = 0 ; d < dim; d++) {
          dweight[i][d] = 4. * pow ( (1. - s / sMax[inode]) * (1. - s / sMin[inode]), 3.) * 
                          ( (- 1. / sMax[inode]) * (1. - s / sMin[inode]) + 
                            (1. - s / sMax[inode]) * (- 1. / sMin[inode]) );
        }
      }
      else{
        weight[i] =1.;
        for (unsigned d = 0 ; d < dim; d++) {
          dweight[i][d] = 0.; 
        }
      }


      //std::cout << i << " " << inode << " " << weight[i] << std::endl ;

      if (weight[i] > 0.) { // take only contribution form the nodes whose weight function overlap with xp
        for (unsigned d = 0 ; d < dim; d++) { // multidimensional loop
          //GetChebyshev (T[i][d], dT[i][d], pOrder, gmpm[p]->_s[i][d] / scale, false);  //1D Chebyshev
          GetPolynomial (T[i][d], dT[i][d], pOrder, gmpm[p]->_s[i][d] / scale, false);  //1D Polynomials
        }
        for (unsigned k = 0; k < aIdx.size(); k++) {
          for (unsigned l = 0; l < aIdx.size(); l++) {
            double TkTl = 1.;
            //std::vector< double > dTkTl (3, 1.);
            for (unsigned d = 0 ; d < dim; d++) {
              TkTl *= T[i][d][aIdx[k][d]] * T[i][d][aIdx[l][d]]; //alpha * beta multidimendional product
//               for (unsigned d2 = 0 ; d2 < dim; d2++) {
//                 if (d == d2) {
//                   dTkTl[d] *= (dT[i][d][aIdx[k][d]] * T[i][d][aIdx[l][d]] + T[i][d][aIdx[k][d]] * dT[i][d][aIdx[l][d]]) ;
//                 }
//                 else {
//                   dTkTl[d] *= T[i][d2][aIdx[k][d2]] * T[i][d2][aIdx[l][d2]];
//                 }
//               }
            }
            Mp[k][l] +=  weight[i] * TkTl;
            for (unsigned d = 0 ; d < dim; d++) {
              //dMp[d][k][l] += weight[i] * dTkTl[d];
              dMp[d][k][l] += dweight[i][d] * TkTl;
            }
          }
        }
      }
      else {
        T[i].resize (0);
      }
    }

    for (unsigned j = 0; j < aIdx.size(); j++) {
      double rhs = 1.;
      for (unsigned d = 0 ; d < dim; d++) {
        rhs *= T0[ aIdx[j][d] ];
      }
      Mp[j][aIdx.size()] = rhs;
      b[j] = rhs;
    }

    //GaussianEleminationWithPivoting (Mp, alpha, false);
    LUsolve (Mp, pivotIndex, b, alpha, false);

    for (unsigned d = 0; d < dim; d++) {
      b.assign (aIdx.size(), 0.);
      b[1] = -1. / scale;
      for (unsigned j = 0; j < aIdx.size(); j++) {
        for (unsigned k = 0; k < aIdx.size(); k++) {
          b[j] -= dMp[d][j][k] * alpha[k];
        }
      }
      LUbackward (Mp, pivotIndex, b, dalpha[d], false);
    }


    double phiSum = 0.;
    std::vector < double > dphiSum (dim, 0.);

    for (unsigned i = 0; i <  gmpm[p]->_node.size(); i++) {

      unsigned inode = gmpm[p]->_node[i];
      std::ostringstream stream;
      stream << "./output/phi" << inode << ".txt";
      fouti.open (stream.str().c_str(), std::ios_base::app);

      for (unsigned d = 0 ; d < dim; d++) {
        std::ostringstream stream;
        stream << "./output/dphi" << inode << "dx" << d << ".txt";
        dfouti[0].open (stream.str().c_str(), std::ios_base::app);
      }


      double P = 1.;
      for (unsigned d = 0 ; d < dim; d++) {
        P *= pow (Xv[inode] , pOrderTest[d]);
      }

      if (weight[i] > 0.) {
        double sumAlphaT = 0.;
        std::vector < double > sumdAlphaT (dim, 0.);
        std::vector < double > sumAlphadT (dim, 0.);
        for (unsigned k = 0; k < aIdx.size(); k++) {
          double Tk = 1;
          std::vector < double > dTk (dim, 1.);
          for (unsigned d = 0 ; d < dim; d++) {
            Tk *= T[i][d][aIdx[k][d]];
            for (unsigned d2 = 0 ; d2 < dim; d2++) {
              if (d == d2) {
                dTk[d] *= dT[i][d][aIdx[k][d]];
              }
              else {
                dTk[d] *= T[i][d2][aIdx[k][d2]];
              }
            }
          }

          sumAlphaT += alpha[k] * Tk;
          for (unsigned d = 0; d < dim; d++) {
            sumdAlphaT[d] += dalpha[d][k] * Tk;
            sumAlphadT[d] += alpha[k] * dTk[d];
          }
        }
        double phi = weight[i] * sumAlphaT;

        std::vector < double > dphi (dim, 0.);
        for (unsigned d = 0; d < dim; d++) {
          dphi[d] = weight[i] * sumdAlphaT[d] + dweight[i][d] * sumAlphaT; //+ weight[i] * sumAlphadT[d];
        }

        Ur[p] += phi * P;
        dUr[p] += dphi[0] * P;

        fouti << gmpm[p]->_xp[0] << " " << phi << std::endl;

        for (unsigned d = 0 ; d < dim; d++) {
          dfouti[0] << gmpm[p]->_xp[0] << " " << dphi[0] << std::endl;
        }

        ymin = (ymin < phi) ? ymin : phi;
        ymax = (ymax > phi) ? ymax : phi;
        phiSum += phi;
        for (unsigned d = 0 ; d < dim; d++) {

          dymin[d] = (dymin[d] < dphi[d]) ? dymin[d] : dphi[d];
          dymax[d] = (dymax[d] > dphi[d]) ? dymax[d] : dphi[d];

          dphiSum[d] += dphi[d];
        }
      }
      fouti.close();
      dfouti[0].close();
      for (unsigned d = 0 ; d < dim; d++) {
        dfouti[d].close();
      }
    }
    fout << gmpm[p]->_xp[0] << " " << phiSum << std::endl;
    for (unsigned d = 0 ; d < dim; d++) {
      dfout[d] << gmpm[p]->_xp[0] << " " << dphiSum[d] << std::endl;
    }
  }
  fout.close();
  for (unsigned d = 0 ; d < dim; d++) {
    dfout[d].close();
  }


  PrintGnuplotScript (Xv[0], Xv[nve - 1], ymin, ymax, nve);
  bool printDerivative = true;
  PrintGnuplotScript (Xv[0], Xv[nve - 1], dymin[0], dymax[0], nve, printDerivative);

  bool test_passed = true;
  for (unsigned p = 0; p < Np; p++) { // particle loop
    Ue[p] = 1.;
    double dUe = 1.;
    for (unsigned d = 0 ; d < dim; d++) {
      Ue[p] *= pow (gmpm[p]->_xp[d], pOrderTest[d]);
    }
    dUe *= pOrderTest[0] * pow (gmpm[p]->_xp[0], pOrderTest[0] - 1);
    if (fabs (Ue[p] - Ur[p]) > 1.0e-6) {
      test_passed = false;
      std::cout << "Error at node = " << p << " exact value = " << Ue[p] << " reconstructed value = " << Ur[p] << std::endl;
      std::cout << "Error at node = " << p << " derivative exact value = " << dUe << " reconstructed value = " << dUr[p] << std::endl;
    }
    if (fabs (dUe - dUr[p]) > 1.0e-6) {
      test_passed = false;
      std::cout << "Error at node = " << p << " exact value = " << Ue[p] << " reconstructed value = " << Ur[p] << std::endl;
      std::cout << "Error at node = " << p << " derivative exact value = " << dUe << " reconstructed value = " << dUr[p] << std::endl;
    }


  }
  if (test_passed == true) std::cout << "Test passed";
  std::cout << std::endl;

  for (unsigned p = 0; p < gmpm.size(); p++) {
    delete gmpm[p] ;
  }

}

void GaussianEleminationWithPivoting (std::vector<std::vector < double > > & A, std::vector < double> &x, const bool &output) {
  // A is nx(n+1) augmented matrix
  unsigned n  = A.size();

  std::vector < unsigned > idx (n);
  for (unsigned i = 0; i < n; i++) {
    idx[i] = i;
  }

  if (output) {
    std::cout << "Before pivoting, A = " << std::endl;
    for (unsigned i = 0; i < n; i++) {
      for (unsigned j = 0; j <= n; j++) {
        std::cout << A[idx[i]][j] << " " ;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  const double tol = 1e-12;
  double max_row = 0.0, abs;
  for (unsigned i = 0; i < n; i++) {
    for (unsigned k = i + 1; k < n; k++) {
      max_row = fabs (A[idx[i]][i]);
      if ( (abs = fabs (A[idx[k]][i])) > max_row) {
        max_row = abs;
        unsigned idxi = idx[i];
        idx[i] = idx[k];
        idx[k] = idxi;
      }
    }
    if (max_row < tol) {
      std::cout << "Degenerate matrix" << std::endl;
      //exit (0);
    }

    for (unsigned j = i + 1; j < n; j++) {
      double mji = A[idx[j]][i] / A[idx[i]][i];
      for (unsigned k = i ; k < n + 1; k++) {
        A[idx[j]][k] -= mji * A[idx[i]][k];
      }
    }
  }
  if (fabs (A[idx[n - 1]][n - 1] / A[idx[n - 1]][n]) < tol) {
    std::cout << "Zero row! Matrix is singular" << std::endl;
    for (unsigned i = 0; i < n; i++) {
      for (unsigned j = 0; j <= n; j++) {
        std::cout << A[idx[i]][j] << " " ;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
    //exit (0);
  }

  x[n - 1] = A[idx[n - 1]][n] / A[idx[n - 1]][n - 1];
  for (unsigned ip1 = n - 1; ip1 > 0; ip1--) {
    unsigned i = ip1 - 1u;
    x[i] = A[idx[i]][n];
    for (int j = ip1; j < n; j++) {
      x[i] -= A[idx[i]][j] * x[j];
    }
    x[i] /= A[idx[i]][i];
  }

  if (output) {
    std::cout << "After pivoting, A = " << std::endl;
    for (unsigned i = 0; i < n; i++) {
      for (unsigned j = 0; j <= n; j++) {
        std::cout << A[idx[i]][j] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl << std::endl;
    std::cout << "Solution is:" << std::endl;
    for (unsigned i = 0; i < n; i++) {
      std::cout << x[i] << std::endl;
    }
  }
}

void LUwithPivoting (std::vector<std::vector < double > > & A, std::vector < unsigned > &idx, const bool &output) {
  unsigned n  = A.size();
  idx.resize (n);
  for (unsigned i = 0; i < n; i++) {
    idx[i] = i;
  }

  if (output) {
    std::cout << "Before pivoting, A = " << std::endl;
    for (unsigned i = 0; i < n; i++) {
      for (unsigned j = 0; j < n; j++) {
        std::cout << A[idx[i]][j] << " " ;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  const double tol = 1e-12;
  double max_row = 0.0, abs;
  for (unsigned i = 0; i < n; i++) {
    for (unsigned k = i + 1; k < n; k++) {
      max_row = fabs (A[idx[i]][i]);
      if ( (abs = fabs (A[idx[k]][i])) > max_row) {
        max_row = abs;
        unsigned idxi = idx[i];
        idx[i] = idx[k];
        idx[k] = idxi;
      }
    }
    if (max_row < tol) {
      std::cout << "Degenerate matrix.";
      //exit (0);
    }

    for (unsigned j = i + 1; j < n; j++) {
      double mji = A[idx[j]][i] / A[idx[i]][i];
      A[idx[j]][i] = mji;
      for (unsigned k = i + 1; k < n; k++) {
        A[idx[j]][k] -= mji * A[idx[i]][k];
      }
    }
  }
  if (fabs (A[idx[n - 1]][n - 1] / A[idx[n - 1]][n]) < tol) {
    std::cout << "Zero row! Matrix is singular" << std::endl;
    for (unsigned i = 0; i < n; i++) {
      for (unsigned j = 0; j < n; j++) {
        std::cout << A[idx[i]][j] << " " ;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    //exit (0);
  }

  if (output) {
    std::cout << "After pivoting, A = " << std::endl;
    for (unsigned i = 0; i < n; i++) {
      for (unsigned j = 0; j < n; j++) {
        std::cout << A[idx[i]][j] << " " ;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}


void LUbackward (const std::vector<std::vector < double > > & A, const std::vector < unsigned > &idx,
                 std::vector < double > &b, std::vector < double > &x, const bool &output) {

  unsigned n  = A.size();

  if (output) {
    std::cout << "RHS before pivoting:" << std::endl;
    for (unsigned i = 0; i < n; i++) {
      std::cout << b[idx[i]] << std::endl;
    }
    std::cout << std::endl;
  }
  for (unsigned i = 0; i < n; i++) {
    for (unsigned j = i + 1; j < n; j++) {
      b[idx[j]] -= A[idx[j]][i] * b[idx[i]];
    }
  }

  if (output) {
    std::cout << "RHS after pivoting:" << std::endl;
    for (unsigned i = 0; i < n; i++) {
      std::cout << b[idx[i]] << std::endl;
    }
    std::cout << std::endl;
  }

  x[n - 1] = b[idx[n - 1]] / A[idx[n - 1]][n - 1];
  for (unsigned ip1 = n - 1; ip1 > 0; ip1--) {
    unsigned i = ip1 - 1u;
    x[i] = b[idx[i]];
    for (int j = ip1; j < n; j++) {
      x[i] -= A[idx[i]][j] * x[j];
    }
    x[i] /= A[idx[i]][i];
  }

  if (output) {
    std::cout << "Solution is:" << std::endl;
    for (unsigned i = 0; i < n; i++) {
      std::cout << x[i] << std::endl;
    }
    std::cout << std::endl;
  }
}

void LUsolve (std::vector<std::vector < double > > & A, std::vector < unsigned > &idx,
              std::vector < double > &b, std::vector < double > &x, const bool &output) {
  LUwithPivoting (A, idx, output);
  LUbackward (A, idx, b, x, output);
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
          std::cout << "alpha[" << index << "][" << j << "]= " << Jp[index][j] << " ";
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


void GetChebyshev (std::vector<double> &T, std::vector<double> &dT, const unsigned &n, const double &x, const bool &output) {
  T.resize (n + 1);
  T[0] = 1;
  T[1] = x;
  for (unsigned i = 2; i < n + 1; i++) {
    T[i] = 2 * x * T[ i - 1] - T[ i - 2];
  }
  if (output) {
    std::cout << "Chebyshev Polynomials at x = " << x << std::endl;
    for (unsigned i = 0; i < n + 1; i++) {
      std::cout << "T" << i << " [x] = " << T[i] << std::endl;
    }
    std::cout << std::endl;
  }

  dT.resize (n + 1);
  dT[0] = 1;
  dT[1] = 2 * x;
  for (unsigned i = 2; i < n + 1; i++) {
    dT[i] = 2 * x * T[ i - 1] - T[ i - 2];
  }

  for (unsigned i = n; i > 0; i--) {
    dT[i] = i * dT[i - 1];
  }
  dT[0] = 0.;

  if (output) {
    std::cout << "Chebyshev Polynomial derivatives at x = " << x << std::endl;
    for (unsigned i = 0; i < n + 1; i++) {
      std::cout << "dT" << i << " [x] = " << dT[i] << std::endl;
    }
    std::cout << std::endl;
  }


}


void GetPolynomial (std::vector<double> &P, std::vector<double> &dP, const unsigned &n, const double &x, const bool &output) {
  P.resize (n + 1);
  P[0] = 1.;

  for (unsigned i = 1; i < n + 1; i++) {
    P[i] = P[i - 1] * x;
  }

  if (output) {
    std::cout << "Chebyshev Polynomials at x = " << x << std::endl;
    for (unsigned i = 0; i < n + 1; i++) {
      std::cout << "P" << i << " [x] = " << P[i] << std::endl;
    }
    std::cout << std::endl;
  }

  dP.resize (n + 1);
  dP[0] = 0.;
  for (unsigned i = 1; i < n + 1; i++) {
    dP[i] = P[i - 1] * i;
  }

  if (output) {
    std::cout << "Polynomial derivatives at x = " << x << std::endl;
    for (unsigned i = 0; i < n + 1; i++) {
      std::cout << "dP" << i << " [x] = " << dP[i] << std::endl;
    }
    std::cout << std::endl;
  }



}


void GetMultiIndex (std::vector <unsigned> &idx, const unsigned &dim, const unsigned& n, const unsigned &i) {
  idx.resize (dim);
  for (unsigned d = 0; d < dim; d++) {
    idx[dim - 1u - d] = (i % static_cast < unsigned > (pow (n, dim - d))) / static_cast < unsigned > (pow (n, dim - 1 - d));
  }
}


void SetElementDofs (std::vector <unsigned> &elementDof, const std::vector < unsigned > & idx, const unsigned & nve1d) {

  unsigned dim = idx.size();
  unsigned nDofs = static_cast < unsigned > (pow (2, dim));

  elementDof.assign (nDofs, 0);
  unsigned sizeHalf = nDofs / 2u;

  unsigned jj;
  for (unsigned d = 0; d < dim; d++) {
    for (unsigned j = 0; j < nDofs; j++) {
      jj = j;
      while (jj >= (2u * sizeHalf)) {
        jj -= 2u * sizeHalf;
      }
      jj /= sizeHalf;

      //elementDof[j] += ( idx[d] + jj ) * ( static_cast <unsigned> ( pow( nve1d , dim - 1u - d)) );
      elementDof[j] += (idx[d] + jj) * (static_cast <unsigned> (pow (nve1d , d)));

    }
    sizeHalf /= 2;
  }
}

void PrintSolution (const std::vector <double> &Ue, const std::vector <double> &Ur, const double *x, const unsigned &dim, const unsigned &n) {
  if (dim > 3) {
    std::cout << "VTK Output: dimension greater than 3 is not supported" << std::endl;
    return;
  }
  std::ofstream fout;
  fout.open ("./output/IMPM.vtk");

  fout << "# vtk DataFile Version 2.0" << std::endl;
  fout << "IMPM example" << std::endl;
  fout << "ASCII" << std::endl;
  fout << "DATASET RECTILINEAR_GRID " << std::endl;
  fout << "DIMENSIONS ";
  unsigned totalNodes = 1;
  for (unsigned k = 0u; k < dim; k++) {
    fout << n << " ";
    //totalNodes *= n;
  }
  if (dim == 2) fout << "1 ";

  if (dim == 1) fout << "1 1 ";
  fout << std::endl;

  fout << "X_COORDINATES " << n << " float " << std::endl;

  for (unsigned i = 0u; i < n; i++) {
    fout << x[i] << " ";
  }
  fout << std::endl;

  fout << "Y_COORDINATES ";
  if (dim > 1) {
    fout << n << " float " << std::endl;
    for (unsigned i = 0u; i < n; i++) {
      fout << x[i] << " ";
    }
  }
  else {
    fout << 1 << " float " << std::endl << 0.;
  }
  fout << std::endl;

  fout << "Z_COORDINATES ";
  if (dim > 2) {
    fout << n << " float " << std::endl;
    for (unsigned i = 0u; i < n; i++) {
      fout << x[i] << " ";
    }
  }
  else {
    fout << 1 << " float " << std::endl << 0.;
  }
  fout << std::endl;


  fout << "POINT_DATA " << Ur.size() << std::endl;
  fout << "SCALARS Ue float 1" << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;
  for (unsigned i = 0u; i < Ue.size(); i++) {
    fout << Ue[i] << " ";
  }

  fout << std::endl;
  fout << "SCALARS Ur float 1" << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;
  for (unsigned i = 0u; i < Ur.size(); i++) {
    fout << Ur[i] << " ";
  }
  fout << std::endl;

  fout << std::endl;
  fout << "SCALARS Ue_minus_Ur float 1" << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;
  for (unsigned i = 0u; i < Ue.size(); i++) {
    fout << Ue[i] - Ur[i] << " ";
  }
  fout << std::endl;

  fout << std::endl;
  fout.close();
}

void PrintGnuplotScript (const double & xmin, const double & xmax, const double & ymin, const double & ymax, const unsigned &nve, const bool &printDerivative) {
  std::ofstream fout;
  if (printDerivative) {
    fout.open ("./output/gnuScriptdPhi.txt");
  }
  else {
    fout.open ("./output/gnuScriptPhi.txt");
  }
  fout << "set xrange[" << xmin << ":" << xmax << "]" << std::endl;
  fout << "set yrange[" << ymin - 0.02 << ":" << ymax + 0.02 << "]" << std::endl;
  fout << "plot ";
  for (unsigned i = 0; i < nve; i++) {
    //if(i==0 || i==nve-1){

    std::ostringstream stream;
    if (printDerivative) {
      stream << "\"dphi" << i << "dx0.txt\"";
      fout << stream.str().c_str() << " u 1:2 title \"d{/Symbol f}_{" << i << "}\" with line,";
    }
    else {
      stream << "\"phi" << i << ".txt\"";
      fout << stream.str().c_str() << " u 1:2 title \"{/Symbol f}_{" << i << "}\" with line,";
    }
  }
  if (printDerivative) {
    fout << "\"dphidx0Sum.txt\" u 1:2 title \"{/Symbol S}_{i}d{/Symbol f}_{i}\" with line,";
  }
  else {
    fout << "\"phiSum.txt\" u 1:2 title \"{/Symbol S}_{i}{/Symbol f}_{i}\" with line,";
  }
  fout << "\npause -1 " << std::endl;
  fout.close();

}


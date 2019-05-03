#include "FemusInit.hpp"

#include "./include/gmpm.hpp"



//double NeumannFactor = 0.;

using namespace femus;

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

  double Xv[9] = {0., 0.15, 0.43, 0.58, 0.81, 0.98, 1.25, 1.425, 1.6};

  double scale = 0.25;

  std::vector < double > sMax (nve1d);
  std::vector < double > sMin (nve1d);

  for (int i = 0; i < nve1d; i++) {
    int deltai = (i % 2 == 0) ? pOrder + nonLocal : pOrder - (!nonLocal);

    int im = (i - deltai >= 0) ?  i - deltai : 0;
    int ip = (i + deltai < nve) ? i + deltai : nve - 1;

    sMin[i] = Xv[im] - Xv[i];
    sMax[i] = Xv[ip] - Xv[i];

    if (nonLocal) {
      if (i - deltai < 0) sMin[i] += Xv[0] - Xv[1];
      if (i + deltai >= nve1d) sMax[i] += Xv[nve1d - 1] - Xv[nve1d - 2];
    }

  }

  unsigned Np = 500;
  //double L = ( Xv[nve1d - 1] - Xv[0]);
  double L = (0.5 * (Xv[nve1d - 3] + Xv[nve1d - 3]) + 0.0025 - Xv[0]);

  double DX = L / Np;

  unsigned maxNumberOfNodes = (2u * (pOrder + 1u) < nve) ? 2u * (pOrder + 1u) : nve;

  std::vector < GMPM *> gmpm (Np);
  for (unsigned p = 0; p < gmpm.size(); p++) {
    gmpm[p] = new GMPM (dim, maxNumberOfNodes);
  }

  std::vector <double> distance (dim);

  maxNumberOfNodes = 0;
  double EA = 10.e06;
  
  for (unsigned p = 0; p < gmpm.size(); p++) {
    for (unsigned d = 0; d < dim; d++) {
      gmpm[p]->_xp[d] = Xv[0] + DX * (0.5 + p);
    }
    for (unsigned j = 0; j < nve; j++) {
      bool particleIsWhithin = false;
      for (unsigned d = 0; d < dim; d++) {
        distance[d] = gmpm[p]->_xp[d] - Xv[j];
        if ( (distance[d] >= 0 && distance[d] <= sMax[j]) || (distance[d] <= 0 && distance[d] >= sMin[j])) {
          particleIsWhithin = true;
        }
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
    
  GMPM gmpmR(dim, maxNumberOfNodes);
  {
    for (unsigned d = 0; d < dim; d++) {
      gmpmR._xp[d] = gmpm[Np - 1]->_xp[d] + DX * 0.5;
    }
    for (unsigned j = 0; j < nve; j++) {
      bool particleIsWhithin = false;
      for (unsigned d = 0; d < dim; d++) {
        distance[d] = gmpmR._xp[d] - Xv[j];
        if ( (distance[d] >= 0 && distance[d] <= sMax[j]) || (distance[d] <= 0 && distance[d] >= sMin[j])) {
          particleIsWhithin = true;
        }
      }
      if (particleIsWhithin) {
        unsigned size = gmpmR._node.size();
        gmpmR._node.resize (size + 1);
        gmpmR._s.resize (size + 1);
        
        gmpmR._node[size] = j;
        gmpmR._s[size].resize (dim);
        for (unsigned d = 0; d < dim; d++) {
          gmpmR._s[size][d] = distance[d];
        }
      }
    }
  
  }
  

  std::cout << "MaxNumberOfNOdes = " << maxNumberOfNodes << std::endl;


  std::vector < std::vector < double > > K (nve);
  for (unsigned i = 0; i < nve; i++) {
    K[i].assign (nve, 0.);
  }
  std::vector < double > F (nve, 0);
  std::vector < double > U (nve, 0);

  std::vector< double> alpha (aIdx.size());
  std::vector< std::vector< double> > dalpha (dim);
  for (unsigned d = 0; d < dim; d++) {
    dalpha[d].resize (aIdx.size());
  }

  std::vector< double> b (aIdx.size());
  std::vector< unsigned> pivotIndex (aIdx.size());

  std::vector < std::vector< double> > Mp (aIdx.size());
  for (unsigned k = 0; k < aIdx.size(); k++) {
    Mp[k].resize (aIdx.size());
  }

  std::vector < std::vector < std::vector< double> > > dMp (dim);
  for (unsigned d = 0; d < dim; d++) {
    dMp[d].resize (aIdx.size());
    for (unsigned k = 0; k < aIdx.size(); k++) {
      dMp[d][k].resize (aIdx.size());
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
  GetPolynomial (T0, dT0, pOrder, 0., false);

  std::vector < double > phi;
  phi.reserve (maxNumberOfNodes);
  std::vector < std::vector < double > > dphi;
  dphi.reserve (maxNumberOfNodes);


  for (unsigned p = 0; p < Np; p++) { // particle loop

    gmpm[p]->GetTestFunction(aIdx, Mp, dMp, nonLocal, weight, sMax, sMin, 
                             dweight, T, dT, T0, dT0, pOrder, scale,
                             b, pivotIndex, alpha, dalpha, phi, dphi);

    for (unsigned i = 0; i <  gmpm[p]->_node.size(); i++) {
      if (weight[i] > 0.) {
        unsigned inode = gmpm[p]->_node[i];
        for (unsigned j = 0; j <  gmpm[p]->_node.size(); j++) {
          if (weight[j] > 0.) {
            unsigned jnode = gmpm[p]->_node[j];
            for (unsigned d = 0; d < dim; d++) {
              K[inode][jnode] += EA * dphi[i][d] * dphi[j][d] * DX;
            }
          }
        }
      }
    }
  }

  double P = -1;

  K[0].assign (nve, 0.);
  K[0][0] = 1;
  
  gmpmR.GetTestFunction(aIdx, Mp, dMp, nonLocal, weight, sMax, sMin, 
                        dweight, T, dT, T0, dT0, pOrder, scale,
                        b, pivotIndex, alpha, dalpha, phi, dphi);
  for (unsigned i = 0; i <  gmpmR._node.size(); i++) {
    unsigned inode = gmpmR._node[i];
    F[inode] -= phi[i];
  }
  
  LUsolve (K, pivotIndex, F, U, true);

  double ULP = 0.;
  

  
  for (unsigned i = 0; i <  gmpmR._node.size(); i++) {
    unsigned inode = gmpmR._node[i];
    ULP += phi[i] * U[inode];
    //std::cout << i << " " << inode << " " << phiLP[i] << " " << U[inode] << std::endl;
  }

  std::cout << "Calculated value = " << ULP << std::endl;
  std::cout << "Exact value = " << P * L / EA << std::endl;

  std::cout << "Relative Error = " << fabs ( (ULP - P * L / EA) / (P * L / EA)) << std::endl;

 

  for (unsigned p = 0; p < Np; p++) { // particle loop
    double Up = 0.;
    gmpm[p]->GetTestFunction(aIdx, Mp, dMp, nonLocal, weight, sMax, sMin, 
                             dweight, T, dT, T0, dT0, pOrder, scale,
                             b, pivotIndex, alpha, dalpha, phi, dphi);
    
    for (unsigned i = 0; i <  gmpm[p]->_node.size(); i++) {
      unsigned inode = gmpm[p]->_node[i];
      Up += phi[i] * U[inode];
    }
    
    std::cout << gmpm[p]->_xp[0] << " " << Up << std::endl;
    
  }
  
  
  
  
  
  
  for (unsigned p = 0; p < gmpm.size(); p++) {
    delete gmpm[p] ;
  }
}





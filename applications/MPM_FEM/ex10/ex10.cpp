#include "FemusInit.hpp"

#include "./include/gmpm.hpp"

using namespace femus;

int main (int argc, char** args) {

  FemusInit (argc, args, MPI_COMM_WORLD);

  bool nonLocal = true;
  bool output = true;
  unsigned Cn = 1;
  unsigned Np = 201;
  unsigned pOrder = 1;
  unsigned dim = 1;
  double scale = 0.25;
   
  std::vector < std::vector <unsigned> > aIdx;
  ComputeIndexSet (aIdx, pOrder, dim, output);

  unsigned nel1d = 4;
  unsigned nel =  static_cast< unsigned > (pow (nel, dim));
  std::vector < double > Xel {0., 0.43, 0.81, 1.25, 1.6};
  
  unsigned nve1d =  nel1d * pOrder + 1u;
  unsigned nve = static_cast< unsigned > (pow (nve1d, dim));

  std::vector < double > Xv (nve1d);

  unsigned counter = 0;
  Xv[counter] = Xel[0];
  counter++;
  for (unsigned i = 0; i < nel1d; i++) {
    double Dx = (Xel[i + 1] - Xel[i]) / pOrder;
    for (unsigned j = 0; j < pOrder; j++) {
      Xv[counter] = Xv[counter - 1] + Dx;
      counter++;
    }
  }
  
  std::vector < double > sMax (nve1d);
  std::vector < double > sMin (nve1d);

  for (int i = 0; i < nve1d; i++) {
    int deltaip = nonLocal * pOrder + ( pOrder  - i % pOrder);
    int deltaim = (i % pOrder == 0) ? nonLocal * pOrder + pOrder : nonLocal * pOrder + i % pOrder;

    int im = (i - deltaim >= 0) ?  i - deltaim : 0;
    int ip = (i + deltaip < nve1d) ? i + deltaip : nve1d - 1;

    sMin[i] = Xv[im] - Xv[i];
    sMax[i] = Xv[ip] - Xv[i];
    
    if (nonLocal) {
      if (i - deltaim < 0) sMin[i] += Xv[0] - Xv[1];
      if (i + deltaip >= nve1d) sMax[i] += Xv[nve1d - 1] - Xv[nve1d - 2];
    }
  }
  
  double L = (Xv[nve1d - pOrder - 1]  - Xv[0]);

  double DX = L / (Np - 1.);
  unsigned maxNumberOfNodes = ( (3u * pOrder + 1u) < nve) ? (3u * pOrder + 1u) : nve;

  std::vector < GMPM *> gmpm (Np);
  for (unsigned p = 0; p < gmpm.size(); p++) {
    gmpm[p] = new GMPM (dim, maxNumberOfNodes);
  }

  std::vector < double > sMaxR = sMax;
  std::vector < double > sMinR = sMin;
  std::vector < double > XvR = Xv;

  sMaxR.resize (nve1d - pOrder);
  sMinR.resize (nve1d - pOrder);
  XvR.resize (nve1d - pOrder);

  maxNumberOfNodes = 0;
  for (unsigned p = 0; p < gmpm.size(); p++) {
    for (unsigned d = 0; d < dim; d++) {
      gmpm[p]->_xp[d] = Xv[0] + DX * p;
    }
    gmpm[p]->CheckIfParticleIsWhithin (XvR, sMinR, sMaxR);
    
    if (maxNumberOfNodes <  gmpm[p]->_node.size()) {
      maxNumberOfNodes = gmpm[p]->_node.size();
    }
  }

  gmpm[0]->SetVolume (DX / 2.);
  for (unsigned p = 1; p < Np - 1; p++) {
    gmpm[p]->SetVolume (DX);
  }
  gmpm[Np - 1]->SetVolume (DX / 2.);

  std::cout << "MaxNumberOfNOdes = " << maxNumberOfNodes << std::endl;

  std::vector < std::vector < double > > K (nve);
  for (unsigned i = 0; i < nve; i++) {
    K[i].assign (nve, 0.);
  }
  std::vector < double > F (nve, 0);

  std::vector < double > phi (maxNumberOfNodes);
  std::vector < std::vector < double > > dphi (maxNumberOfNodes);

  double EA = 1.e07;
  double P = -1.;


  if (output) {
    std::ofstream fout;
    for (unsigned i = 0; i < XvR.size(); i++) {
      std::ostringstream stream;
      stream << "./output/phi" << i << ".txt";
      fout.open (stream.str().c_str());
      fout.close();

      fout.open ("./output/phiSum.txt");
      fout.close();

      for (unsigned d = 0; d < dim; d++) {
        std::ostringstream stream;
        stream << "./output/dphi" << i << "dx" << d << ".txt";
        fout.open (stream.str().c_str());
        fout.close();

        stream << "./output/dphi" << "dx" << d << "Sum.txt";
        fout.open (stream.str().c_str());
        fout.close();
      }
    }
    PrintGrid (Xv);
    PrintGnuplotScript (Xv[0], Xv[nve - 1], -0.5, 1.5, XvR.size());
    bool printDerivative = true;
    PrintGnuplotScript (Xv[0], Xv[nve - 1], -10., 10., XvR.size(), printDerivative);
  }

  double weight;
  for (unsigned p = 0; p < Np; p++) { // particle loop

    WindowFunction wf;
    wf.BuildWeight (XvR, pOrder, gmpm[p]->_xp[0], nonLocal, Cn);

    gmpm[p]->GetTestFunction (aIdx, nonLocal, XvR, sMaxR, sMinR, pOrder, scale, wf, phi, dphi, weight);

    for (unsigned i = 0; i <  gmpm[p]->_node.size(); i++) {
      unsigned inode = gmpm[p]->_node[i];
      for (unsigned j = 0; j <  gmpm[p]->_node.size(); j++) {
        unsigned jnode = gmpm[p]->_node[j];
        for (unsigned d = 0; d < dim; d++) {
          K[inode][jnode] += EA * dphi[i][d] * dphi[j][d] * weight;
        }
      }
      if (p == Np - 1) { //boundary force
        F[inode] = P * phi[i];
      }
    }

    if (output) {
      std::ofstream fout;
      double sumPhi = 0.;
      std::vector<double> sumdPhi (dim, 0.);
      for (unsigned i = 0; i <  gmpm[p]->_node.size(); i++) {
        unsigned inode = gmpm[p]->_node[i];
        std::ostringstream stream;

        stream << "./output/phi" << inode << ".txt";
        fout.open (stream.str().c_str(), std::ios_base::app);
        fout << gmpm[p]->_xp[0] << " " << phi[i] << std::endl;
        fout.close();

        sumPhi += phi[i];
        for (unsigned d = 0 ; d < dim; d++) {
          std::ostringstream stream;
          stream << "./output/dphi" << inode << "dx" << d << ".txt";
          fout.open (stream.str().c_str(), std::ios_base::app);
          fout << gmpm[p]->_xp[0] << " " << dphi[i][d] << std::endl;
          fout.close();
          sumdPhi[d] += 0;//dphi[i][d];
        }
      }
      fout.open ("./output/phiSum.txt", std::ios_base::app);
      fout << gmpm[p]->_xp[0] << " " << sumPhi << std::endl;
      fout.close();
      for (unsigned d = 0 ; d < dim; d++) {
        std::ostringstream stream;
        stream << "./output/dphi" << "dx" << d << "Sum.txt";
        fout.open (stream.str().c_str(), std::ios_base::app);
        fout << gmpm[p]->_xp[0] << " " << sumdPhi[d] << std::endl;
        fout.close();
      }
    }
  }

  K[0].assign (nve, 0.);
  K[0][0] = 1;
  for (unsigned i = 0; i < nve; i++) {
    if (K[i][i] == 0.)  {
      K[i][i] = 1.;
      F[i] = 0.;
    }
  }

  std::vector < double > U (nve);
  std::vector< unsigned> pivotIndex;
  LUsolve (K, pivotIndex, F, U, true);

  double Ur = 0.;

  for (unsigned i = 0; i <  gmpm[Np - 1]->_node.size(); i++) {
    unsigned inode = gmpm[Np - 1]->_node[i];
    Ur += phi[i] * U[inode];
  }

  std::cout << "Calculated value = " << Ur << std::endl;
  std::cout << "Exact value = " << P * L / EA << std::endl;
  std::cout << "Relative Error = " << fabs ( (Ur - P * L / EA) / (P * L / EA)) << std::endl;
  if (output) {
    std::ofstream fout;
    fout.open ("./output/solution.txt");
    for (unsigned p = 0; p < Np; p++) { // particle loop
      double Up = 0.;

      WindowFunction wf;
      wf.BuildWeight (XvR, pOrder, gmpm[p]->_xp[0], nonLocal,  Cn );

      gmpm[p]->GetTestFunction (aIdx, nonLocal, XvR, sMaxR, sMinR, pOrder, scale, wf, phi, dphi, weight);

      for (unsigned i = 0; i <  gmpm[p]->_node.size(); i++) {
        unsigned inode = gmpm[p]->_node[i];
        Up += phi[i] * U[inode];
      }
      fout << gmpm[p]->_xp[0] << " " << Up << std::endl;
    }
    fout.close();
    fout.open ("./output/gnuSolutionScript.txt");
    fout << "set xrange[" << Xv[0] << ":" << Xv[nve1d - 1] << "]" << std::endl;
    fout << "set yrange[ 0 :" << 1.2 * P * L / EA << "]" << std::endl;
    fout << "plot \"solution.txt\" u 1:2 with line" << std::endl;
    fout << "pause -1 " << std::endl;
    fout.close();
  }

  for (unsigned p = 0; p < gmpm.size(); p++) {
    delete gmpm[p] ;
  }

}





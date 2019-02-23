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

bool sparse = true;
//bool sparse = false;

unsigned alpha = 7;
unsigned M = pow (10, alpha);   //number of samples
unsigned N = 3; //dimension of the parameter space (each of the M samples has N entries)
unsigned L = 9; //max refinement level
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

void PrintBaseCoefficients (const std::vector < double > & c, const std::vector <unsigned> &I, const std::vector < unsigned > &n);
unsigned GetBaseIndex (const std::vector <unsigned> &IL, const std::vector <unsigned> &I);
unsigned GetBaseIndex (const std::vector <unsigned> &IL, const std::vector < std::vector <unsigned> > &iiL);


void GetBaseIndexes (const std::vector <unsigned> &IL, const unsigned &i, std::vector <unsigned> &I);


unsigned GetCoarseBaseIndex (const std::vector <unsigned> &IL, const std::vector <unsigned> &JL, std::vector <unsigned> &I);


void IncreaseLevelIndex (const unsigned &Lm1, const unsigned &N, std::vector <unsigned> &IL);

unsigned SumIndexes (std::vector <unsigned> &IL) {
  unsigned sum = 0u;
  for (unsigned k = 0u; k < IL.size(); k++) {
    sum += IL[k];
  }
  return sum;
}

void PrintVTKHistogram (const std::vector <double> &C, const std::vector <double> &Cr, const std::vector <unsigned > &IL, const std::vector < unsigned > &n);

int main (int argc, char** argv) {

  //BEGIN construction of the sample set
  std::vector < std::vector < double > >  samples;
  samples.resize (M);

  for (unsigned m = 0u; m < M; m++) {
    samples[m].resize (N);
  }

  for (unsigned m = 0u; m < M; m++) {
    for (unsigned n = 0u; n < N; n++) {
      double var = var_nor();
      double varunif = var_unif();
      double U = var_unif1();
      samples[m][n] = var;
      //samples[m][n] = varunif;
    }
  }
  //END

  const unsigned Lm1 = L - 1u;
  std::vector< unsigned > n1D (L);
  n1D[0u] = 1u;
  for (unsigned i = 1u; i < n1D.size(); i++) {
    n1D[i] = n1D[i - 1u] * 2u;
  }

  std::vector <unsigned> IL (N, Lm1);

  //Build Histogram

  clock_t time0 = clock();

  std::vector <double> cI (static_cast< unsigned > (pow (n1D[Lm1], N)), 0.);
  std::vector < unsigned > I (N);

  double h = (xmax - xmin) / n1D[Lm1];
  for (unsigned m = 0u; m < M; m++) {
    for (unsigned k = 0u; k < N; k++) {
      double x = (samples[m][k] - xmin) / h;
      if (x < 0.) I[k] = 0u;
      else if (x >= n1D[Lm1]) I[k] = n1D[Lm1] - 1u;
      else I[k] = static_cast < unsigned > (floor (x));
    }
    cI[GetBaseIndex (IL, I)]++;
  }

  double vol = pow (h, N);
  vol *= M;
  for (unsigned i = 0u; i < cI.size(); i++) {
    cI[i] = cI[i] / vol;
  }

  std::cout << std::endl << "HISTOGRAM TIME:\t" << static_cast<double> ( (clock() - time0)) / CLOCKS_PER_SEC << std::endl;

//   std::cout << "Histogram " << std::endl;
//   PrintBaseCoefficients (cI, IL, n1D);

  ///////////////////////////////////////////////////////////////////////////////

  time0 = clock();

  //Build Hierarchical Bases

  std::vector <std::vector <double> > CH (static_cast< unsigned > (pow (L, N)));
  std::vector <std::vector <unsigned> > CHu (static_cast< unsigned > (pow (L, N)));
  unsigned cnt = 0u;

  std::vector<std::vector<unsigned>> JT (static_cast< unsigned > (pow (L, N)));

  IL.assign (N, 0u);
  for (unsigned l = 0u; l < CHu.size(); l++) {
    unsigned sum = SumIndexes (IL);
    if (!sparse || sum < L) {
      JT[cnt] = IL;
      unsigned spaceSize = 1u;
      for (unsigned k = 0u; k < N; k++) {
        spaceSize <<= IL[k];
      }
      CH[cnt].assign (spaceSize, 0.);
      CHu[cnt].assign (spaceSize, 0u);
      cnt++;
    }
    IncreaseLevelIndex (Lm1, N, IL);
  }
  JT.resize (cnt);
  CH.resize (cnt);
  CHu.resize (cnt);
  std::cout << std::endl << "Total number of Hierarchical spaces: " << JT.size() << std::endl;

  double H = (xmax - xmin);
  std::vector < std::vector <unsigned> > iiL (N);
  for (unsigned k = 0u; k < N; k++) {
    iiL[k].assign (L, 0u);
  }

  unsigned nL = n1D[Lm1];
  double nLdH = n1D[Lm1] / H;

  for (unsigned m = 0u; m < M; m++) {
    for (unsigned k = 0u; k < N; k++) {
      iiL[k][Lm1] = static_cast < unsigned > (floor ( (samples[m][k] - xmin) * nLdH));
      if (iiL[k][Lm1] >= nL) iiL[k][Lm1] = nL - 1u;
      for (unsigned j = Lm1; j > 1u; j--) {
        iiL[k][j - 1] = (iiL[k][j] >> 1);
      }
    }
    for (unsigned l = 0u; l < JT.size(); l++) {
      CHu[l][GetBaseIndex (JT[l], iiL)]++;
    }
  }

  for (unsigned l = 0u; l < JT.size(); l++) {
    double vol = 1.;
    for (unsigned k = 0u; k < N; k++) {
      vol *= H / n1D[JT[l][k]];
    }
    vol *= M;
    for (unsigned i = 0u; i < CH[l].size(); i++) {
      CH[l][i] = CHu[l][i] / vol;
    }
  }

  std::cout << std::endl << "HIERARCHICAL BASE STAGE 1 TIME:\t" << static_cast<double> ( (clock() - time0)) / CLOCKS_PER_SEC << std::endl;
  clock_t time1 = clock();

  for (unsigned l = 0u; l < JT.size(); l++) {
    for (unsigned m = 0u; m < l; m++) {
      bool subSpace = true;
      for (unsigned k = 0u; k < N; k++) {
        if (JT[m][k] > JT[l][k]) {
          subSpace = false;
          break;
        }
      }
      if (subSpace) {
        for (unsigned  i = 0u; i < CH[l].size(); i++) {
          GetBaseIndexes (JT[l], i, I);
          CH[l][i] -= CH[m][ GetCoarseBaseIndex (JT[l], JT[m], I) ];
        }
      }
    }
  }
  std::cout << std::endl << "HIERARCHICAL BASE STAGE 2 TIME:\t" << static_cast<double> ( (clock() - time1)) / CLOCKS_PER_SEC << std::endl;
  std::cout << std::endl << "HIERARCHICAL BASE TOTAL TIME:\t" << static_cast<double> ( (clock() - time0)) / CLOCKS_PER_SEC << std::endl;

//   std::cout << std::endl;
//   std::cout << "Hierarchical bases " << std::endl;
//   for (unsigned l = 0u; l < JT.size(); l++) {
//   PrintBaseCoefficients (CH[l], JT[l], n1D);
//   }

//////////////////////////////////////////////////////////////////////

//Reconstruct Histogram from Hierarchical Bases

  time0 = clock();
  std::vector <double> cIr (static_cast< unsigned > (pow (n1D[Lm1], N)), 0);
  IL.assign (N, Lm1);
  for (unsigned i = 0u; i < cIr.size(); i++) {
    GetBaseIndexes (IL, i, I);
    for (unsigned l = 0u; l < JT.size(); l++) {
      cIr[i] += CH[l][GetCoarseBaseIndex (IL, JT[l], I)];
    }
  }
  std::cout << std::endl << "RECONSTRUCTION TIME:\t" << static_cast<double> ( (clock() - time0)) / CLOCKS_PER_SEC << std::endl;

//   std::cout << "Histogram reconstructed from Hierarchical Bases" << std::endl;
//   IL.assign (N, Lm1);
//   PrintBaseCoefficients (cIr, IL, n1D);
//
//   std::cout << "Difference between Histogram and Histogram reconstructed from Hierarchical Bases" << std::endl;
//   for (unsigned i = 0u; i < cIr.size(); i++) cIr[i] -= cI[i];
//   PrintBaseCoefficients (cIr, IL, n1D);

  double error = 0.;
  for (unsigned i = 0u; i < cIr.size(); i++)
    error += fabs (cIr[i] - cI[i]);

  std::cout << "Difference between Histogram and Histogram reconstructed from Hierarchical Bases = " << error << std::endl;

  PrintVTKHistogram (cI, cIr, IL, n1D);


  return 0u;

} //end main


void PrintBaseCoefficients (const std::vector < double > & c, const std::vector <unsigned> &I, const std::vector < unsigned > &n) {

  for (unsigned k = 0; k < N; k++) {
    std::cout << "I[" << k << "] = " << I[k] << ", ";
  }
  std::cout << std::endl;

  unsigned N = I.size();
  unsigned Nm1 = N - 1u;
  unsigned Nm2 = N - 2u;
  unsigned nNm1 = n[ I[Nm1] ];
  for (unsigned i = 0u; i < c.size(); i++) {
    std::cout << c[i] << " ";
    unsigned size = nNm1;
    unsigned ip1 = i + 1u;
    for (unsigned k = 0u; k < N; k++) {
      if (ip1 % size == 0u) std::cout << std::endl;
      if (k < Nm1) size *= n[ I[Nm2 - k]];
    }
  }
}

unsigned GetBaseIndex (const std::vector <unsigned> &IL, const std::vector <unsigned> &I) {
  unsigned i = I[0u];
  for (unsigned k = 1u; k < I.size(); k++) {
    i = (i << IL[k]) + I[k];
  }
  return i;
}

unsigned GetBaseIndex (const std::vector <unsigned> &IL, const std::vector < std::vector <unsigned> > &iiL) {
  unsigned i = iiL[0u][IL[0u]];
  unsigned N = IL.size();
  for (unsigned k = 1u; k < N; k++) {
    i = (i << IL[k]) + iiL[k][IL[k]];
  }
  return i;
}


void GetBaseIndexes (const std::vector <unsigned> &IL, const unsigned &i, std::vector <unsigned> &I) {
  unsigned N = IL.size();
  unsigned ii = i, jj;
  for (unsigned k = 0u; k < N; k++) {
    jj = ii >> IL[N - 1u - k]  ;
    I[N - 1u - k] = ii - (jj << IL[N - 1u - k]);
    ii = jj;
  }
}

unsigned GetCoarseBaseIndex (const std::vector <unsigned> &IL, const std::vector <unsigned> &JL, std::vector <unsigned> &I) {

  unsigned i = (I[0] >> (IL[0] - JL[0]));
  for (unsigned k = 1u; k < I.size(); k++) {
    i = (i << JL[k]) + (I[k] >> (IL[k] - JL[k]));
  }
  return i;
}

void PrintVTKHistogram (const std::vector <double> &C, const std::vector <double> &Cr, const std::vector <unsigned > &IL, const std::vector < unsigned > &n) {
  if (IL.size() > 3) {
    std::cout << "VTK Output: dimension greater than 3 is not supported" << std::endl;
    return;
  }
  std::ofstream fout;
  fout.open ("./output/Histogram.vtk");

  fout << "# vtk DataFile Version 2.0" << std::endl;
  fout << "Volume example" << std::endl;
  fout << "ASCII" << std::endl;
  fout << "DATASET STRUCTURED_POINTS " << std::endl;
  fout << "DIMENSIONS ";
  unsigned totalNodes = 1;
  for (unsigned k = 0u; k < IL.size(); k++) {
    fout << n[IL[k]] << " ";
    totalNodes *= n[IL[k]];
  }
  if (IL.size() == 2) fout << "1 ";
  if (IL.size() == 1) fout << "1 1 ";
  fout << std::endl;
  fout << "ASPECT_RATIO 1 1 1" << std::endl;
  fout << "ORIGIN 0 0 0" << std::endl;
  fout << "POINT_DATA " << totalNodes << std::endl;
  fout << "SCALARS Histogram float 1" << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;
  for (unsigned i = 0u; i < C.size(); i++) {
    fout << C[i] << " ";
  }

  fout << std::endl;
  fout << "SCALARS Reconstructed float 1" << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;
  for (unsigned i = 0u; i < C.size(); i++) {
    fout << Cr[i] << " ";
  }
  fout << std::endl;

  fout << std::endl;
  fout << "SCALARS Difference float 1" << std::endl;
  fout << "LOOKUP_TABLE default" << std::endl;
  for (unsigned i = 0u; i < C.size(); i++) {
    fout << C[i] - Cr[i] << " ";
  }
  fout << std::endl;

  fout << std::endl;
  fout.close();
}

void IncreaseLevelIndex (const unsigned &Lm1, const unsigned &N, std::vector <unsigned> &IL) {
  unsigned* ptr = &IL[N - 1u];
  for (int k = 0; k < N; k++, ptr--) {
    if ( (*ptr) == Lm1) {
      (*ptr) = 0u;
    }
    else {
      (*ptr) += 1u;
      return;
    }
  }
}



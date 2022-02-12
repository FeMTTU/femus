
#include <iostream>
#include <fstream>

#include <boost/math/special_functions/factorials.hpp>

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

class WindowFunction {
  public:
    void BuildWeight (const std::vector <double> &Xv, const unsigned &pOrder, const double &x, const bool &nonLocal, const unsigned &n);

    double x0;
    double x1;

    double w0;
    double w1;

    double dw0;
    double dw1;

    const static double ak[5][5];
};

const double WindowFunction::ak[5][5] = {{1.}, {1., 2.}, {1., 3., 6.}, {1., 4., 10., 20.}, {1., 5., 15., 35., 70.}};

void WindowFunction::BuildWeight (const std::vector <double> &Xv, const unsigned &pOrder, const double &x, const bool &nonLocal, const unsigned &Cn) {

  unsigned i = pOrder;
  while (x > Xv[i]) i += pOrder;
  x0 = Xv[i - pOrder];
  x1 = Xv[i];

  if (nonLocal) {
      
    double l0 = (x - x1) / (x0 - x1);
    double l1 = (x - x0) / (x1 - x0);

    double dl0 = 1. / (x0 - x1);
    double dl1 = -dl0;

    unsigned n = Cn + 1;

    double Pl0 = ak[Cn][0];
    double Pl1 = ak[Cn][0];
    double dPl0 = 0.;
    double dPl1 = 0.;

    double l0k = 1.;
    double l1k = 1.;
    for (unsigned k = 1; k < n; k++) {
      dPl0 += k * ak[Cn][k] * l0k * dl0;
      dPl1 += k * ak[Cn][k] * l1k * dl1;
      l0k *= l0;
      l1k *= l1;
      Pl0 += ak[Cn][k] * l0k;
      Pl1 += ak[Cn][k] * l1k;
    }
    dw0 = l0k * (n * dl0 * Pl1 + l0 * dPl1);
    dw1 = l1k * (n * dl1 * Pl0 + l1 * dPl0);

    l0k *= l0;
    l1k *= l1;
    w0 =  l0k * Pl1;
    w1 =  l1k * Pl0;


    // std::cout << l0 << " " << l1 << " " << dl0 << " " <<dl1 << std::endl;
    // std::cout << w0 << " " << w1 << " " << dw0 << " " <<dw1 << std::endl;

  }
  else {
    w0 = w1 = dw0 = dw1 = 0.;
  }
}

class GMPM {
  public:
    GMPM (const unsigned &dim, const unsigned &Nv) {
      _dim = dim;
      _xp.resize (dim);
      _node.reserve (Nv);
      _s.reserve (Nv);
    };

    void CheckIfParticleIsWhithin (const std::vector < double > &Xv,
                                   const std::vector < double > &sMin,
                                   const std::vector < double > &sMax);

    void GetTestFunction (const std::vector < std::vector <unsigned> > &aIdx,
                          const bool &nonLocal,
                          const std::vector < double > &Xv,
                          const std::vector < double > &sMax,
                          const std::vector < double > &sMin,
                          const unsigned &pOrder,
                          const double &scale,
                          const WindowFunction &w,
                          std::vector <double > &phi,
                          std::vector < std::vector < double > > &dphi,
                          double & weight
                         );

    unsigned _dim;
    std::vector < unsigned > _node;
    std::vector < double > _xp;
    std::vector < std::vector < double > > _s;
    void SetVolume (const double &volume) {
      _volume = volume;
    }
  private:
    double _volume;
    static std::vector < double > _T0;
    static std::vector < double > _dT0;
    static std::vector < std::vector< double> > _Mp;
    static std::vector < std::vector < std::vector< double> > > _dMp;
    static std::vector < double > _weight;
    static std::vector < std::vector < double > > _dweight;
    static std::vector < std::vector < std::vector < double > > > _T;
    static std::vector < std::vector < std::vector < double > > > _dT;
    static std::vector< double> _b;
    static std::vector< unsigned> _pivotIndex;
    static std::vector< double> _alpha;
    static std::vector< std::vector< double> > _dalpha;
};

std::vector < double > GMPM::_T0;
std::vector < double > GMPM::_dT0;
std::vector < std::vector< double> > GMPM::_Mp;
std::vector < std::vector < std::vector< double> > > GMPM::_dMp;
std::vector < double > GMPM::_weight;
std::vector < std::vector < double > > GMPM::_dweight;
std::vector < std::vector < std::vector < double > > > GMPM::_T;
std::vector < std::vector < std::vector < double > > > GMPM::_dT;
std::vector< double> GMPM::_b;
std::vector< unsigned> GMPM::_pivotIndex;
std::vector< double> GMPM::_alpha;
std::vector< std::vector< double> > GMPM::_dalpha;

void GMPM::CheckIfParticleIsWhithin (const std::vector < double > &Xv,
                                     const std::vector < double > &sMin,
                                     const std::vector < double > &sMax) {
  std::vector < double > s (_dim);
  for (unsigned j = 0; j < Xv.size(); j++) {
    bool particleIsWhithin = true;
    for (unsigned d = 0; d < _dim; d++) {
      s[d] = _xp[d] - Xv[j];
      if ( (s[d] >= 0 && s[d] > sMax[j]) || (s[d] <= 0 && s[d] < sMin[j])) {
        particleIsWhithin = false;
        break;
      }
    }
    if (particleIsWhithin) {
      unsigned size = _node.size();
      _node.resize (size + 1);
      _s.resize (size + 1);

      _node[size] = j;
      _s[size].resize (_dim);
      for (unsigned d = 0; d < _dim; d++) {
        _s[size][d] = s[d];
      }
    }
  }
}


void GMPM::GetTestFunction (const std::vector < std::vector <unsigned> > &aIdx,
                            const bool &nonLocal,
                            const std::vector < double > &Xv,
                            const std::vector < double > &sMax,
                            const std::vector < double > &sMin,
                            const unsigned &pOrder,
                            const double &scale,
                            const WindowFunction &w,
                            std::vector <double > &phi,
                            std::vector < std::vector < double > > &dphi,
                            double & weight
                           ) {

  weight = _volume;
  GetPolynomial (_T0, _dT0, pOrder, 0., false);

  _Mp.resize (aIdx.size());
  for (unsigned k = 0; k < aIdx.size(); k++) {
    _Mp[k].assign (aIdx.size(), 0.);
  }

  _dMp.resize (_dim);
  for (unsigned d = 0; d < _dim; d++) {
    _dMp[d].resize (aIdx.size());
    for (unsigned k = 0; k < aIdx.size(); k++) {
      _dMp[d][k].assign (aIdx.size(), 0.);
    }
  }

  _alpha.resize (aIdx.size());
  _dalpha.resize (_dim);
  for (unsigned d = 0; d < _dim; d++) {
    _dalpha[d].resize (aIdx.size());
  }

  _b.resize (aIdx.size());
  _pivotIndex.resize (aIdx.size());

  _T.resize (_node.size());
  _dT.resize (_node.size());
  for (unsigned i = 0; i < _node.size(); i++) {
    _T[i].resize (_dim);
    _dT[i].resize (_dim);
  }

  _weight.assign (_node.size(), 0.);
  _dweight.resize (_node.size());
  for (unsigned i = 0; i < _node.size(); i++) {
    _dweight[i].resize (_dim);
  }

  bool pIsNode = false;
  unsigned pnode = 0;
  if(nonLocal){
    for (unsigned i = 0; i < _node.size(); i++) {
        unsigned inode = _node[i];
        if (inode % pOrder == 0) {
        if (fabs (_xp[0] - Xv[inode]) < 1.e-14) {
            pIsNode = true;
            pnode = inode;
        }
      }
    }
  }

  for (unsigned i = 0; i < _node.size(); i++) {
    unsigned inode = _node[i];
    //if (nonLocal) {

    double x = Xv[inode];
    if (pIsNode == false) {
    
      if (x < w.x0 - 1.e-14) {
        _weight[i] = w.w0;
        _dweight[i][0] = w.dw0;
      }
      else if (x <= w.x1 + 1.e-14) {
        _weight[i] = 1.;
        _dweight[i][0] = 0.;
      }
      else {
        _weight[i] = w.w1;
        _dweight[i][0] = w.dw1;
      }
    }
    else {
      if (pnode > inode) {
        _weight[i] = (pnode - inode > pOrder) ? 0. : 1.;
      }
      else {
        _weight[i] = (inode - pnode > pOrder) ? 0. : 1.;
      }
      _dweight[i][0] = 0;
    }

//     if (fabs (_xp[0] - 0.43) < 1.0e-8) {
//       std::cout << x << " " << w.x0 << " " << w.x1 << " " << w.w0 << " " << w.w1 << " " << w.dw0 << " " << w.dw1 << std::endl;
//       std::cout << inode << " " << _weight[i] << " " << _dweight[i][0] << std::endl << std::flush;
//     }


//       double s = _s[i][0];
//       _weight[i] = pow ( (1. - s / sMax[inode]) * (1. - s / sMin[inode]), 4.);
//       for (unsigned d = 0 ; d < _dim; d++) {
//         _dweight[i][d] = 4. * pow ( (1. - s / sMax[inode]) * (1. - s / sMin[inode]), 3.) *
//                          ( (- 1. / sMax[inode]) * (1. - s / sMin[inode]) +
//                            (1. - s / sMax[inode]) * (- 1. / sMin[inode]));
//       }
    // }
    // else {
    //   _weight[i] = 1.;
    //   for (unsigned d = 0 ; d < _dim; d++) {
    //     _dweight[i][d] = 0.;
    //   }
    // }

    //std::cout << inode << " " << _weight[i] <<std::endl;

    if (_weight[i] > 0.) { // take only contribution form the nodes whose weight function overlap with xp
      for (unsigned d = 0 ; d < _dim; d++) { // multi_dimensional loop
        GetPolynomial (_T[i][d], _dT[i][d], pOrder, _s[i][d] / scale, false); //1D Polynomials
      }
      for (unsigned k = 0; k < aIdx.size(); k++) {
        for (unsigned l = 0; l < aIdx.size(); l++) {
          double TkTl = 1.;
          for (unsigned d = 0 ; d < _dim; d++) {
            TkTl *= _T[i][d][aIdx[k][d]] * _T[i][d][aIdx[l][d]]; //alpha * beta multi_dimendional product
          }
          _Mp[k][l] +=  _weight[i] * TkTl;
          for (unsigned d = 0 ; d < _dim; d++) {
            _dMp[d][k][l] += _dweight[i][d] * TkTl;
          }
        }
      }
    }
    else {
      _T[i].resize (0);
    }
  }

  for (unsigned j = 0; j < aIdx.size(); j++) {
    double rhs = 1.;
    for (unsigned d = 0 ; d < _dim; d++) {
      rhs *= _T0[ aIdx[j][d] ];
    }
    _b[j] = rhs;
  }
  LUsolve (_Mp, _pivotIndex, _b, _alpha, false);

  for (unsigned d = 0; d < _dim; d++) {
    _b.assign (aIdx.size(), 0.);
    _b[1] = -1. / scale;
    for (unsigned j = 0; j < aIdx.size(); j++) {
      for (unsigned k = 0; k < aIdx.size(); k++) {
        _b[j] -= _dMp[d][j][k] * _alpha[k];
      }
    }
    LUbackward (_Mp, _pivotIndex, _b, _dalpha[d], false);
  }

  phi.assign (_node.size(), 0.);
  dphi.resize (_node.size());

  for (unsigned i = 0; i < _node.size(); i++) {

    dphi[i].assign (_dim, 0.);

    unsigned inode = _node[i];

    if (_weight[i] > 0.) {
      double sumAlphaT = 0.;
      std::vector < double > sumdAlphaT (_dim, 0.);
      std::vector < double > sumAlphadT (_dim, 0.);
      for (unsigned k = 0; k < aIdx.size(); k++) {
        double Tk = 1;
        std::vector < double > dTk (_dim, 1.);
        for (unsigned d = 0 ; d < _dim; d++) {
          Tk *= _T[i][d][aIdx[k][d]];
          for (unsigned d2 = 0 ; d2 < _dim; d2++) {
            if (d == d2) {
              dTk[d] *= _dT[i][d][aIdx[k][d]];
            }
            else {
              dTk[d] *= _T[i][d2][aIdx[k][d2]];
            }
          }
        }

        sumAlphaT += _alpha[k] * Tk;
        for (unsigned d = 0; d < _dim; d++) {
          sumdAlphaT[d] += _dalpha[d][k] * Tk;
          sumAlphadT[d] += _alpha[k] * dTk[d];
        }
      }
      phi[i] = _weight[i] * sumAlphaT;
      for (unsigned d = 0; d < _dim; d++) {
        dphi[i][d] = _weight[i] * sumdAlphaT[d] + _dweight[i][d] * sumAlphaT;
      }
    }
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
      std::cout << b[i] << std::endl;
    }
    std::cout << std::endl;
  }
  for (unsigned i = 0; i < n; i++) {
    for (unsigned j = i + 1; j < n; j++) {
      b[idx[j]] -= A[idx[j]][i] * b[idx[i]];
    }
  }

  const double tol = 1e-12;
  if (fabs (A[idx[n - 1]][n - 1] / b[idx[n - 1]]) < tol) {
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
  fout << "\"grid.txt\" u 1:2 title \"\" with linespoints,";
  fout << "\npause -1 " << std::endl;
  fout.close();

}

void PrintGrid (const std::vector<double> &Xv) {
  std::ofstream fout;
  fout.open ("./output/grid.txt");
  for (unsigned i = 0; i < Xv.size(); i++) {
    fout << Xv[i] << " " << 0. << std::endl;
  }
  fout.close();
}




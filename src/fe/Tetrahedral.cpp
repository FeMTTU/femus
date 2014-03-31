/*=========================================================================

 Program: FEMUS
 Module: Tetrahedral
 Authors: Eugenio Aulisa
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------

#include "Basis.hpp"


double tet1::eval_phi(const int *I,const double* x) const {
  return tet_1(x[0],x[1],x[2],I[0],I[1],I[2]);
}

double tet1::eval_dphidx(const int *I,const double* x) const {
  return dtet_1dx(x[0],x[1],x[2],I[0],I[1],I[2]);
}

double tet1::eval_dphidy(const int *I,const double* x) const {
  return dtet_1dy(x[0],x[1],x[2],I[0],I[1],I[2]);
}

double tet1::eval_dphidz(const int *I,const double* x) const {
  return dtet_1dz(x[0],x[1],x[2],I[0],I[1],I[2]);
}

double tet1::tet_1(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const {
  return (!i*!j*!k)*(1.-x-y-z)+ !(i-2)*x + !(j-2)*y + !(k-2)*z;
}

double tet1::dtet_1dx(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const {
  return -(!i*!j*!k) + !(i-2);
}

double tet1::dtet_1dy(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const {
  return -(!i*!j*!k) + !(j-2);
}

double tet1::dtet_1dz(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const {
  return -(!i*!j*!k) + !(k-2);
}

//************************************************************

double  tet2::eval_phi(const int *I,const double* x) const {
  return tet_2(x[0],x[1],x[2],I[0],I[1],I[2]);
}

double  tet2::eval_dphidx(const int *I,const double* x) const {
  return dtet_2dx(x[0],x[1],x[2],I[0],I[1],I[2]);
}

double  tet2::eval_dphidy(const int *I,const double* x) const {
  return dtet_2dy(x[0],x[1],x[2],I[0],I[1],I[2]);
}


double  tet2::eval_dphidz(const int *I,const double* x) const {
  return dtet_2dz(x[0],x[1],x[2],I[0],I[1],I[2]);
}

double tet2::tet_2(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const {
  double f=0.;
  double t=x+y+z;
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=(1.-t)*(1.-2.*t);
        break;
      case 1:
        f=4.*z*(1.-t);
        break;
      case 2:
        f=(-z+2.*z*z);
        break;
      }
      break;
    case 1:
      switch (k) {
      case 0:
        f=4.*y*(1.-t);
        break;
      case 1:
        f=4.*y*z;
        break;
      }
      break;
    case 2:
      f=(-y+2.*y*y);
      break;
    }
    break;
  case 1:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=4.*x*(1.-t);
        break;
      case 1:
        f=4.*x*z;
        break;
      }
      break;
    case 1:
      f=4.*x*y;
      break;
    }
    break;
  case 2:
    f=(-x+2.*x*x);
    break;
  }
  return f;
}

double tet2::dtet_2dx(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const {
  double f=0.;
  double t=x+y+z;
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=(-3.+4.*t);
        break;
      case 1:
        f=-4.*z;
        break;
      case 2:
        f=0.;
        break;
      }
      break;
    case 1:
      switch (k) {
      case 0:
        f=-4.*y;
        break;
      case 1:
        f=0.;
        break;
      }
      break;
    case 2:
      f=0;
      break;
    }
    break;
  case 1:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=4.*(1.-t-x);
        break;
      case 1:
        f=4.*z;
        break;
      }
      break;
    case 1:
      f=4.*y;
      break;
    }
    break;
  case 2:
    f=(-1.+4.*x);
    break;
  }
  return f;
}

double tet2::dtet_2dy(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const {
  double f=0.;
  double t=x+y+z;
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=(-3.+4.*t);
        break;
      case 1:
        f=-4.*z;
        break;
      case 2:
        f=0;
        break;
      }
      break;
    case 1:
      switch (k) {
      case 0:
        f=4.*(1.-t-y);
        break;
      case 1:
        f=4.*z;
        break;
      }
      break;
    case 2:
      f=(-1.+4.*y);
      break;
    }
    break;
  case 1:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=-4.*x;
        break;
      case 1:
        f=0;
        break;
      }
      break;
    case 1:
      f=4.*x;
      break;
    }
    break;
  case 2:
    f=0.;
    break;
  }
  return f;
}

double tet2::dtet_2dz(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const {
  double f=0.;
  double t=x+y+z;
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=(-3.+4.*t);
        break;
      case 1:
        f=4.*(1.-t-z);
        break;
      case 2:
        f=(-1.+4.*z);
        break;
      }
      break;
    case 1:
      switch (k) {
      case 0:
        f=-4.*y;
        break;
      case 1:
        f=4.*y;
        break;
      }
      break;
    case 2:
      f=0.;
      break;
    }
    break;
  case 1:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=-4.*x;
        break;
      case 1:
        f=4.*x;
        break;
      }
      break;
    case 1:
      f=0.;
      break;
    }
    break;
  case 2:
    f=0.;
    break;
  }
  return f;
}

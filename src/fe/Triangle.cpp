#include "Basis.hpp"
//************************************************************

double tri1::eval_phi(const int *I,const double* x) const {
  return tri1a(x[0],x[1],I[0],I[1]);
}

double tri1::eval_dphidx(const int *I,const double* x) const {
  return dtri1dx(x[0],x[1],I[0],I[1]);
}

double tri1::eval_dphidy(const int *I,const double* x) const {
  return dtri1dy(x[0],x[1],I[0],I[1]);
}

double tri1::eval_dphidz(const int *I,const double* x) const {
  return 0;
}

double tri1::tri1a(const double& x, const double& y, const int& i,const int& j) const {
  return (!i*!j)*(1.-x-y)+ !(i-2)*x + !(j-2)*y;
}

double tri1::dtri1dx(const double& x, const double& y, const int& i,const int& j) const {
  return -(!i*!j) + !(i-2);
}

double tri1::dtri1dy(const double& x, const double& y, const int& i,const int& j) const {
  return -(!i*!j) + !(j-2);
}

//************************************************************

double tri2::eval_phi(const int *I,const double* x) const {
  return tri2a(x[0],x[1],I[0],I[1]);
}

double tri2::eval_dphidx(const int *I,const double* x) const {
  return dtri2dx(x[0],x[1],I[0],I[1]);
}

double tri2::eval_dphidy(const int *I,const double* x) const {
  return dtri2dy(x[0],x[1],I[0],I[1]);
}

double tri2::eval_dphidz(const int *I,const double* x) const {
  return 0;
}

double tri2::tri2a(const double& x, const double& y, const int& i,const int& j) const {
  double f=0.;
  switch (i) {
  case(0):
    switch (j) {
    case(0):
      f=(1.-x-y)*(1.-2.*x-2.*y);
      break;
    case(1):
      f=4.*y*(1.-x-y);
      break;
    case(2):
      f=(-y+2.*y*y);
      break;
    }
    break;
  case(1):
    switch (j) {
    case(0):
      f=4.*x*(1.-x-y);
      break;
    case(1):
      f=4.*x*y;
      break;
    }
    break;
  case(2):
    f=(-x+2.*x*x);
    break;
  }
  return f;
}

double tri2::dtri2dx(const double& x, const double& y, const int& i,const int& j) const {
  double f=0.;
  switch (i) {
  case(0):
    switch (j) {
    case(0):
      f=(-3.+4.*x+4.*y);
      break;
    case(1):
      f=-4.*y;
      break;
    case(2):
      f=0;
      break;
    }
    break;
  case(1):
    switch (j) {
    case(0):
      f=4.*(1-2.*x-y);
      break;
    case(1):
      f=4.*y;
      break;
    }
    break;
  case(2):
    f=-1.+4.*x;
    break;
  }
  return f;
}


double tri2::dtri2dy(const double& x, const double& y, const int& i,const int& j) const {
  double f=0.;
  switch (i) {
  case(0):
    switch (j) {
    case(0):
      f=(-3.+4.*x+4.*y);
      break;
    case(1):
      f=4.*(1-x-2.*y);
      break;
    case(2):
      f=-1.+4.*y;
      break;
    }
    break;
  case(1):
    switch (j) {
    case(0):
      f=-4.*x;
      break;
    case(1):
      f=4.*x;
      break;
    }
    break;
  case(2):
    f=0.;
    break;
  }
  return f;
}


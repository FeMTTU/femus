#include "Basis.hpp"
//************************************************************

double wedge1::eval_phi(const int *I,const double* x) const {
  return tri1(x[0],x[1],I[0],I[1])*lag1(x[2],I[2]);
}

double wedge1::eval_dphidx(const int *I,const double* x) const {
  return dtri1dx(x[0],x[1],I[0],I[1])*lag1(x[2],I[2]);
}

double wedge1::eval_dphidy(const int *I,const double* x) const {
  return dtri1dy(x[0],x[1],I[0],I[1])*lag1(x[2],I[2]);
}

double wedge1::eval_dphidz(const int *I,const double* x) const {
  return tri1(x[0],x[1],I[0],I[1])*dlag1(x[2],I[2]);
}

double wedge1::tri1(const double& x, const double& y, const int& i,const int& j) const {
  return (!i*!j)*(1.-x-y)+ !(i-2)*x + !(j-2)*y;
}

double wedge1::dtri1dx(const double& x, const double& y, const int& i,const int& j) const {
  return -(!i*!j) + !(i-2);
}

double wedge1::dtri1dy(const double& x, const double& y, const int& i,const int& j) const {
  return -(!i*!j) + !(j-2);
}

double wedge1::lag1(const double& x, const int& i) const {
  return (!i)*0.5*(1.-x)+!(i-2)*0.5*(1.+x);
}

double wedge1::dlag1(const double& x, const int& i) const {
  return (!i)*(-0.5)+!(i-2)*0.5;
}

//************************************************************

double wedge2::eval_phi(const int *I,const double* x) const {
  return tri2(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
}

double wedge2::eval_dphidx(const int *I,const double* x) const {
  return dtri2dx(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
}

double wedge2::eval_dphidy(const int *I,const double* x) const {
  return dtri2dy(x[0],x[1],I[0],I[1])*lag2(x[2],I[2]);
}

double wedge2::eval_dphidz(const int *I,const double* x) const {
  return tri2(x[0],x[1],I[0],I[1])*dlag2(x[2],I[2]);
}

double wedge2::lag2(const double& x, const int& i) const {
  return !i*0.5*x*(x-1.) + !(i-1)*(1.-x)*(1.+x) + !(i-2)*0.5*x*(1.+x);
}

double wedge2::dlag2(const double& x, const int& i) const {
  return !i*(x-0.5) + !(i-1)*(-2.*x) + !(i-2)*(x+0.5);
}

double wedge2::tri2(const double& x, const double& y, const int& i,const int& j) const {
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

double wedge2::dtri2dx(const double& x, const double& y, const int& i,const int& j) const {
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


double wedge2::dtri2dy(const double& x, const double& y, const int& i,const int& j) const {
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


//************************************************************
double  wedgeth::eval_phi(const int *I,const double* x) const {
  return wed_th(x[0],x[1],x[2],I[0],I[1],I[2]);
}

double  wedgeth::eval_dphidx(const int *I,const double* x) const {
  return dwed_thdx(x[0],x[1],x[2],I[0],I[1],I[2]);
}

double  wedgeth::eval_dphidy(const int *I,const double* x) const {
  return dwed_thdy(x[0],x[1],x[2],I[0],I[1],I[2]);
}

double  wedgeth::eval_dphidz(const int *I,const double* x) const {
  return dwed_thdz(x[0],x[1],x[2],I[0],I[1],I[2]);
}

double wedgeth::wed_th(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const {
  double f=0.;
  double t=x+y;
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=(1-t)*(-2*t-z)*(1-z)*0.5;
        break;
      case 1:
        f=(1-t)*(1-z*z);
        break;
      case 2:
        f=(1-t)*(-2*t+z)*(1+z)*0.5;
        break;
      }
      break;
    case 1:
      switch (k) {
      case 0:
        f=2*y*(1-z)*(1-t);
        break;
      case 2:
        f=2*y*(1+z)*(1-t);
        break;
      }
      break;
    case 2:
      switch (k) {
      case 0:
        f=y*(-2+2*y-z)*(1-z)*0.5;
        break;
      case 1:
        f=y*(1-z*z);
        break;
      case 2:
        f=y*(-2+2*y+z)*(1+z)*0.5;
        break;
      }
      break;
    }
    break;
  case 1:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=2*x*(1-z)*(1-t);
        break;
      case 2:
        f=2*x*(1+z)*(1-t);
        break;
      }
      break;
    case 1:
      switch (k) {
      case 0:
        f=2*x*y*(1-z);
        break;
      case 2:
        f=2*x*y*(1+z);
        break;
      }
      break;
    }
    break;
  case 2:
    switch (k) {
    case 0:
      f=x*(-2+2*x-z)*(1-z)*0.5;
      break;
    case 1:
      f=x*(1-z*z);
      break;
    case 2:
      f=x*(-2+2*x+z)*(1+z)*0.5;
      break;
    }
    break;
  }
  return f;
}


double wedgeth::dwed_thdx(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const {
  double f=0.;
  double t=x+y;
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=(-1.+2.*t+0.5*z)*(1.-z);
        break;
      case 1:
        f=(z*z-1.);
        break;
      case 2:
        f=(-1.+2.*t-0.5*z)*(1.+z);
        break;
      }
      break;
    case 1:
      switch (k) {
      case 0:
        f=-2.*y*(1-z);
        break;
      case 2:
        f=-2.*y*(1+z);
        break;
      }
      break;
    case 2:
      switch (k) {
      case 0:
        f=0.;
        break;
      case 1:
        f=0.;
        break;
      case 2:
        f=0.;
        break;
      }
      break;
    }
    break;
  case 1:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=2*(1.-z)*(1.-2.*x-y);
        break;
      case 2:
        f=2*(1.+z)*(1.-2.*x-y);
        break;
      }
      break;
    case 1:
      switch (k) {
      case 0:
        f=2.*y*(1.-z);
        break;
      case 2:
        f=2.*y*(1.+z);
        break;
      }
      break;
    }
    break;
  case 2:
    switch (k) {
    case 0:
      f=(-2.+4.*x-z)*(1.-z)*0.5;
      break;
    case 1:
      f=(1.-z*z);
      break;
    case 2:
      f=(-2.+4.*x+z)*(1+z)*0.5;
      break;
    }
    break;
  }
  return f;
}


double wedgeth::dwed_thdy(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const {
  double f=0.;
  double t=x+y;
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=(-1.+2.*t+0.5*z)*(1.-z);
        break;
      case 1:
        f=(z*z-1.);
        break;
      case 2:
        f=(-1.+2.*t-0.5*z)*(1.+z);
        break;
      }
      break;
    case 1:
      switch (k) {
      case 0:
        f=2.*(1-z)*(1.-x-2.*y);
        break;
      case 2:
        f=2.*(1+z)*(1.-x-2.*y);
        break;
      }
      break;
    case 2:
      switch (k) {
      case 0:
        f=(-2.+4.*y-z)*(1.-z)*0.5;
        break;
      case 1:
        f=(1.-z*z);
        break;
      case 2:
        f=(-2.+4.*y+z)*(1.+z)*0.5;
        break;
      }
      break;
    }
    break;
  case 1:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=-2.*x*(1.-z);
        break;
      case 2:
        f=-2.*x*(1.+z);
        break;
      }
      break;
    case 1:
      switch (k) {
      case 0:
        f=2.*x*(1.-z);
        break;
      case 2:
        f=2.*x*(1.+z);
        break;
      }
      break;
    }
    break;
  case 2:
    switch (k) {
    case 0:
      f=0.;
      break;
    case 1:
      f=0.;
      break;
    case 2:
      f=0.;
      break;
    }
    break;
  }
  return f;
}


double wedgeth::dwed_thdz(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const {
  double f=0.;
  double t=1.-(x+y);
  switch (i) {
  case 0:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=t*(0.5-t+z);
        break;
      case 1:
        f=-2.*t*z;
        break;
      case 2:
        f=t*(-0.5+t+z);
        break;
      }
      break;
    case 1:
      switch (k) {
      case 0:
        f=-2.*y*t;
        break;
      case 2:
        f=2.*y*t;
        break;
      }
      break;
    case 2:
      switch (k) {
      case 0:
        f=y*(0.5-y+z);
        break;
      case 1:
        f=-2.*y*z;
        break;
      case 2:
        f=y*(-0.5+y+z);
        break;
      }
      break;
    }
    break;
  case 1:
    switch (j) {
    case 0:
      switch (k) {
      case 0:
        f=-2.*x*t;
        break;
      case 2:
        f=2.*x*t;
        break;
      }
      break;
    case 1:
      switch (k) {
      case 0:
        f=-2.*x*y;
        break;
      case 2:
        f=2.*x*y;
        break;
      }
      break;
    }
    break;
  case 2:
    switch (k) {
    case 0:
      f=x*(0.5-x+z);
      break;
    case 1:
      f=-2.*x*z;
      break;
    case 2:
      f=x*(-0.5+x+z);
      break;
    }
    break;
  }
  return f;
}



// double wedgeth::wed_dthdz(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const{
//   double f;
//   double t=x+y;
//   switch(i){
//   case 0:
//     switch(j){
//     case 0:
//       switch(k){
//       case 0:
// 	f=(1-t)*(-0.5+t+z);
// 	break;
//       case 1:
// 	f=-(1-t)*2*z;
// 	break;
//       case 2:
// 	f=(1-t)*(0.5-t+z);
// 	break;
//       }
//       break;
//     case 1:
//       switch(k){
//       case 0:
// 	f=-2*y*(1-t);
// 	break;
//       case 2:
// 	f=2*y*(1-t);
// 	break;
//       }
//       break;
//     case 2:
//       switch(k){
//       case 0:
// 	f=y*(0.5-y+z);
// 	break;
//       case 1:
// 	f=-y*2*z;
// 	break;
//       case 2:
// 	f=y*(-2+2*y+z)*(1+z)*0.5;
// 	break;
//       }
//       break;
//     }
//     break;
//   case 1:
//     switch(j){
//     case 0:
//       switch(k){
//       case 0:
// 	f=2*x*(1-z)*(1-t);
// 	break;
//       case 2:
// 	f=2*x*(1+z)*(1-t);
// 	break;
//       }
//       break;
//     case 1:
//       switch(k){
//       case 0:
// 	f=2*x*y*(1-z);
// 	break;
//       case 2:
// 	f=2*x*y*(1+z);
// 	break;
//       }
//       break;
//     }
//     break;
//   case 2:
//     switch(k){
//     case 0:
//       f=x*(-2+2*x-z)*(1-z)*0.5;
//       break;
//     case 1:
//       f=x*(1-z*z);
//       break;
//     case 2:
//       f=x*(-2+2*x+z)*(1+z)*0.5;
//       break;
//     }
//     break;
//   }
//   return f;
// }

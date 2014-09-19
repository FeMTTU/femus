/*=========================================================================

 Program: FEMUS
 Module: basis
 Authors: Eugenio Aulisa
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
 * This class contains the fe basis function and their derivatives
*/


#ifndef __basis_h___
#define __basis_h___


namespace femus {
  
  
class basis {
public:
  virtual double eval_phi(const int *I,const double* x) const = 0;
  virtual double eval_dphidx(const int *I,const double* x) const = 0;
  virtual double eval_dphidy(const int *I,const double* x) const = 0;
  virtual double eval_dphidz(const int *I,const double* x) const = 0;
};

//************************************************************

class hex1: public basis {
private:
  double lag1(const double& x, const int& i) const;
  double dlag1(const double& x, const int& i) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class hex2: public basis {
private:
  double lag2(const double& x, const int& i) const;
  double dlag2(const double& x, const int& i) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class hexth: public basis {
private:
  double th2(const double& x,  const int& i) const;
  double dth2(const double& x,  const int& i) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

class hex0: public basis {
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//******************************************************************************

class hexpwl: public basis {
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************


class wedge1: public basis {
  double lag1(const double& x, const int& i) const;
  double dlag1(const double& x, const int& i) const;
  double tri1(const double& x,const double& y, const int& j,const int& i) const;
  double dtri1dx(const double& x,const double& y, const int& j,const int& i) const;
  double dtri1dy(const double& x,const double& y, const int& j,const int& i) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class wedge2: public basis {
  double lag2(const double& x, const int& i) const;
  double dlag2(const double& x, const int& i) const;
  double tri2(const double& x,const double& y, const int& j,const int& i) const;
  double dtri2dx(const double& x,const double& y, const int& j,const int& i) const;
  double dtri2dy(const double& x,const double& y, const int& j,const int& i) const;
public:
  virtual double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class wedgeth: public basis {
  double wed_th(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const;
  double dwed_thdx(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const;
  double dwed_thdy(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const;
  double dwed_thdz(const double& x, const double& y, const double &z, const int& i,const int& j, const int & k) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class tet0: public basis {
private:
  double tet_0(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const;
  double dtet_0dx(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const;
  double dtet_0dy(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const;
  double dtet_0dz(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class tet1: public basis {
private:
  double tet_1(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const;
  double dtet_1dx(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const;
  double dtet_1dy(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const;
  double dtet_1dz(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class tet2: public basis {
private:
  double tet_2(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const;
  double dtet_2dx(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const;
  double dtet_2dy(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const;
  double dtet_2dz(const double& x, const double& y, const double& z, const int & i,const int & j,const int &k) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//******************************************************************************

class quad0: public basis {
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//******************************************************************************

class quadpwl: public basis {
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class quad1: public basis {
private:
  double lag1(const double& x, const int& i) const;
  double dlag1(const double& x, const int& i) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class quad2: public basis {
private:
  double lag2(const double& x, const int& i) const;
  double dlag2(const double& x, const int& i) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class quadth: public basis {
private:
  double th2(const double& x,  const int& i) const;
  double dth2(const double& x,  const int& i) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class tri0: public basis {
  double tri0a(const double& x,const double& y, const int& j,const int& i) const;
  double dtri0dx(const double& x,const double& y, const int& j,const int& i) const;
  double dtri0dy(const double& x,const double& y, const int& j,const int& i) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class tri1: public basis {
  double tri1a(const double& x,const double& y, const int& j,const int& i) const;
  double dtri1dx(const double& x,const double& y, const int& j,const int& i) const;
  double dtri1dy(const double& x,const double& y, const int& j,const int& i) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class tri2: public basis {
  double tri2a(const double& x,const double& y, const int& j,const int& i) const;
  double dtri2dx(const double& x,const double& y, const int& j,const int& i) const;
  double dtri2dy(const double& x,const double& y, const int& j,const int& i) const;
public:
  virtual double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class line0: public basis {
  double lag0(const double& x, const int& i) const;
  double dlag0(const double& x, const int& i) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class line1: public basis {
  double lag1(const double& x, const int& i) const;
  double dlag1(const double& x, const int& i) const;
public:
  double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};

//************************************************************

class line2: public basis {
  double lag2(const double& x, const int& i) const;
  double dlag2(const double& x, const int& i) const;
public:
  virtual double eval_phi(const int *I,const double* x) const;
  double eval_dphidx(const int *I,const double* x) const;
  double eval_dphidy(const int *I,const double* x) const;
  double eval_dphidz(const int *I,const double* x) const;
};







} //end namespace femus


#endif

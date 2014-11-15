/*=========================================================================

  Program: FEMUS
  Module: basis
  Authors: Eugenio Aulisa, Giorgio Bornia
 
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

#include <iostream>
#include <stdlib.h>  

namespace femus {
    
  class basis {
  public:
    
    const int _nc,_nf,_ncf0,_ncf1,_ncf2;
    const unsigned _type;
    basis(const unsigned &type, const int &nc,const int &nf,const int &ncf0, const int &ncf1, const int &ncf2):
      _type(type),
      _nc(nc),
      _nf(nf),
      _ncf0(ncf0),
      _ncf1(ncf1),
      _ncf2(ncf2){ };
           
    virtual void PrintType() const = 0 ;
    virtual double eval_phi(const int *I,const double* x) const {
      std::cout<<"Error this phi is not available for this element \n"; abort();
    };
    virtual double eval_dphidx(const int *I,const double* x) const {
      std::cout<<"Error this dphidx is not available for this element dimension\n"; abort();
    };
    virtual double eval_dphidy(const int *I,const double* x) const {
      std::cout<<"Error this dphidy is not available for this element dimension\n"; abort();
    };
    virtual double eval_dphidz(const int *I,const double* x) const {
      std::cout<<"Error this dphidz is not available for this element dimension\n"; abort();
    };
    virtual double eval_d2phidx2(const int *I,const double* x) const {
      std::cout<<"Error this d2phix2 is not available for this element dimension\n"; abort();
    };
    virtual double eval_d2phidy2(const int *I,const double* x) const {
      std::cout<<"Error this d2phidy2 is not available for this element dimension\n"; abort();
    };
    virtual double eval_d2phidz2(const int *I,const double* x) const {
      std::cout<<"Error this d2phidz2 is not available for this element dimension\n"; abort();
    };
    virtual double eval_d2phidxdy(const int *I,const double* x) const {
      std::cout<<"Error this d2phidxdy is not available for this element dimension\n"; abort();
    };
    virtual double eval_d2phidydz(const int *I,const double* x) const {
      std::cout<<"Error this d2phidydz is not available for this element dimension\n"; abort();
    };
    virtual double eval_d2phidzdx(const int *I,const double* x) const {
      std::cout<<"Error this d2phidydz is not available for this element dimension\n"; abort();
    };
    
    virtual const double* getX(const int &i) const = 0;
    virtual const int* getIND(const int &i) const = 0;
    virtual const int* getKVERT_IND(const int &i) const = 0;
   
    
    
  protected:  
    //1D basis
    // linear lagrangian
    inline double lag1(const double& x, const int& i) const {
      return (!i)*0.5*(1.-x)+!(i-2)*0.5*(1.+x);
    }
    inline double dlag1(const double& x, const int& i) const {
      return (!i)*(-0.5)+!(i-2)*0.5;
    }
  
    //quadratic lagrangian  
    inline double th2(const double& x, const int& i) const {
      return !i*(0.5)*(1.-x) + !(i-1)*(1.-x)*(1.+x) + !(i-2)*(0.5)*(1.+x);
    }

    inline double dth2(const double& x, const int& i) const {
      return (!i)*(-0.5) + !(i-1)*(-2.*x) + !(i-2)*(0.5);
    }

    inline double d2th2(const double& x, const int& i) const {
      return !(i-1)*(-2.);
    }
   
    //bi-quadratic lagrangian  
    inline double lag2(const double& x, const int& i) const {
      return !i*0.5*x*(x-1.) + !(i-1)*(1.-x)*(1.+x) + !(i-2)*0.5*x*(1.+x);
    }

    inline double dlag2(const double& x, const int& i) const {
      return !i*(x-0.5) + !(i-1)*(-2.*x) + !(i-2)*(x+0.5);
    }

    inline double d2lag2(const double& x, const int& i) const{
      return !i + !(i-1)*(-2.) + !(i-2);
    }
  
    //2D basis 
    // linear triangle 
    inline double triangle1(const double& x, const double& y, const int& i,const int& j) const {
      return (!i*!j)*(1.-x-y)+ !(i-2)*x + !(j-2)*y;
    }

    inline double dtriangle1dx(const double& x, const double& y, const int& i,const int& j) const {
      return -(!i*!j) + !(i-2);
    }

    inline double dtriangle1dy(const double& x, const double& y, const int& i,const int& j) const {
      return -(!i*!j) + !(j-2);
    }
  
    // quadratic triangle 
    inline double triangle2(const double& x, const double& y, const int& i,const int& j) const {
      return 
        !i     * (!j* (1.-x-y)*(1.-2.*x-2.*y) + !(j-1)* 4.*y*(1.-x-y) + !(j-2)*(-y+2.*y*y)) +
	!(i-1) * (!j* 4.*x*(1.-x-y) + !(j-1) * 4.*x*y) +
	!(i-2) * (!j*(-x+2.*x*x));  
    }

    inline double dtriangle2dx(const double& x, const double& y, const int& i,const int& j) const {
      return 
        !i     * (!j* (-3.+4.*x+4.*y) + !(j-1)*y*(-4.) ) +
	!(i-1) * (!j* 4.*(1.-2.*x-y)  + !(j-1)*y*(4.)) +
	!(i-2) * (!j*(-1 + 4.*x));
    }


    inline double dtriangle2dy(const double& x, const double& y, const int& i,const int& j) const {
      return 
        !j     * (!i* (-3.+4.*y+4.*x) + !(i-1)*x*(-4.) ) +
	!(j-1) * (!i* 4.*(1.-2.*y-x)  + !(i-1)*x*(4.)) +
	!(j-2) * (!i*(-1 + 4.*y));
    }

    inline double d2triangle2dx2(const double& x, const double& y, const int& i,const int& j) const {
      return !j*( (!i)*4. +!(i-1)*(-8.) + !(i-2)*4. );
    }

    inline double d2triangle2dy2(const double& x, const double& y, const int& i,const int& j) const {
      return !i*( (!j)*4. +!(j-1)*(-8.) + !(j-2)*4. );
    }

    inline double d2triangle2dxdy(const double& x, const double& y, const int& i,const int& j) const {
      return ( (!i)*(!j) + !(i-1)*!(j-1) )*4. + ( !(i-1)*(!j) + (!i)*!(j-1) )*(-4.);
    }
   
  };

  //************************************************************
    
  class hex_lag : public basis {  
  public:  
    hex_lag(const unsigned int& type, const int& nc, const int& nf):
      basis(type, nc, nf, 8, 20, 27){ };
     
    const double* getX(const int &i) const{return X[i];};
    const int* getIND(const int &i) const{return IND[i];};
    const int* getKVERT_IND(const int &i) const {return KVERT_IND[i];};
    
  protected: 
    static const double X[125][3];
    static const int IND[27][3];
    static const int KVERT_IND[125][2];
  };
  
  class hex1: public hex_lag {
  public:
    hex1(): hex_lag(0, 8, 27) {}; 
    void PrintType() const { std::cout<<" hex1 ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const;
    double eval_d2phidydz(const int *I,const double* x) const;
    double eval_d2phidzdx(const int *I,const double* x) const;
  };

  //************************************************************

  class hexth: public hex_lag { 
  public:
    hexth(): hex_lag(1, 20, 81) {};    
    void PrintType() const { std::cout<<" hexth ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const;
    double eval_d2phidy2(const int *I,const double* x) const;
    double eval_d2phidz2(const int *I,const double* x) const;
    double eval_d2phidxdy(const int *I,const double* x) const;
    double eval_d2phidydz(const int *I,const double* x) const;
    double eval_d2phidzdx(const int *I,const double* x) const;
  
  };

  //************************************************************

  class hex2: public hex_lag {
  public:
     hex2(): hex_lag(2, 27, 125) {};    
    void PrintType() const { std::cout<<" hex2 ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const;
    double eval_d2phidy2(const int *I,const double* x) const;
    double eval_d2phidz2(const int *I,const double* x) const;
    double eval_d2phidxdy(const int *I,const double* x) const;
    double eval_d2phidydz(const int *I,const double* x) const;
    double eval_d2phidzdx(const int *I,const double* x) const;
  
  };

  //************************************************************

  class hex_const : public basis {
  public: 
    hex_const(const unsigned int& type, const int& nc, const int& nf):
      basis(type, nc, nf, 8, 20, 27){ };
    
    const double* getX(const int &i) const{return X[i];};
    const int* getIND(const int &i) const{return IND[i];};
    const int* getKVERT_IND(const int &i) const {return KVERT_IND[i];};
    
  protected: 
    static const double X[32][3];
    static const int IND[4][3];
    static const int KVERT_IND[32][2];
  };
  
  
  class hex0: public hex_const {
  public:
    hex0(): hex_const(17, 1, 8) {};    
    void PrintType() const { std::cout<<" hex0 ";};
    
    double eval_phi(const int *I,const double* x) const{ return 1.; };
    double eval_dphidx(const int *I,const double* x) const{ return 0.; };
    double eval_dphidy(const int *I,const double* x) const{ return 0.; };
    double eval_dphidz(const int *I,const double* x) const{ return 0.; };
    
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
  };

  //******************************************************************************

  class hexpwl: public hex_const {
  public:
    hexpwl(): hex_const(18, 4, 32) {};     
    void PrintType() const { std::cout<<" hexpwl ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
  
  };

  //************************************************************
  
  class wedge_lag : public basis{
  public:
    wedge_lag(const unsigned int& type, const int& nc, const int& nf):
      basis(type, nc, nf, 6, 15, 18){ };
    
    const double* getX(const int &i) const{return X[i];};
    const int* getIND(const int &i) const{return IND[i];};
    const int* getKVERT_IND(const int &i) const {return KVERT_IND[i];};
    
  protected: 
    static const double X[75][3];
    static const int IND[18][3];
    static const int KVERT_IND[75][2];
  };

  class wedge1: public wedge_lag {
  public:
    wedge1(): wedge_lag(3, 6, 18) {};    
    void PrintType() const { std::cout<<" wedge1 ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
  
  };

  //************************************************************

  class wedgeth: public wedge_lag {
  public:
    wedgeth(): wedge_lag(4, 15, 57) {};    
    void PrintType() const { std::cout<<" wedgeth ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const;
    double eval_d2phidy2(const int *I,const double* x) const;
    double eval_d2phidz2(const int *I,const double* x) const;
    double eval_d2phidxdy(const int *I,const double* x) const;
    double eval_d2phidydz(const int *I,const double* x) const;
    double eval_d2phidzdx(const int *I,const double* x) const;
  
  };

  //************************************************************

  class wedge2: public wedge_lag {
  
  public:
    wedge2(): wedge_lag(5, 18, 75) {};    
    void PrintType() const { std::cout<<" wedge2 ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const;
    double eval_d2phidy2(const int *I,const double* x) const;
    double eval_d2phidz2(const int *I,const double* x) const;
    double eval_d2phidxdy(const int *I,const double* x) const;
    double eval_d2phidydz(const int *I,const double* x) const;
    double eval_d2phidzdx(const int *I,const double* x) const;
  
  };

  class wedge0: public basis {
  public:
    wedge0(): basis(20, 1, 8, 6, 15, 18){ };    
    void PrintType() const { std::cout<<" wedge0 ";};
    
    double eval_phi(const int *I,const double* x) const{ return 1.; };
    double eval_dphidx(const int *I,const double* x) const{ return 0.; };
    double eval_dphidy(const int *I,const double* x) const{ return 0.; };
    double eval_dphidz(const int *I,const double* x) const{ return 0.; };
    
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
    
    const double* getX(const int &i) const{std::cout<<"Abort in Wedge0\n"; abort();};
    const int* getIND(const int &i) const{std::cout<<"Abort in Wedge0\n"; abort();};
    const int* getKVERT_IND(const int &i) const {std::cout<<"Abort in Wedge0\n"; abort();};
  };
  
  
  //************************************************************

  class tet_lag : public basis{
  public:
    tet_lag(const unsigned int& type, const int& nc, const int& nf):
      basis(type, nc, nf, 4, 10, 10){ };
    const double* getX(const int &i) const{return X[i];};
    const int* getIND(const int &i) const{return IND[i];};
    const int* getKVERT_IND(const int &i) const {return KVERT_IND[i];};
  
  protected: 
    static const double X[35][3];
    static const int IND[10][3];
    static const int KVERT_IND[35][2];
  };
  
  
  
  class tet1: public tet_lag {
  public:
    tet1(): tet_lag(6, 4, 10) {};    
    void PrintType() const { std::cout<<" tet1 ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
  };

  //************************************************************

  class tet2: public tet_lag {
  public:
    tet2(): tet_lag(7, 10, 35) {};
    void PrintType() const { std::cout<<" tet2 ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const;
    double eval_d2phidy2(const int *I,const double* x) const;
    double eval_d2phidz2(const int *I,const double* x) const;
    double eval_d2phidxdy(const int *I,const double* x) const;
    double eval_d2phidydz(const int *I,const double* x) const;
    double eval_d2phidzdx(const int *I,const double* x) const;
  };

  //************************************************************

  class tet0: public basis {
  public:
    tet0(): basis(30, 1, 8, 4, 10, 10){ };
    void PrintType() const { std::cout<<" tet0 ";};  
    
    double eval_phi(const int *I,const double* x) const{ return 1.; };
    double eval_dphidx(const int *I,const double* x) const{ return 0.; };
    double eval_dphidy(const int *I,const double* x) const{ return 0.; };
    double eval_dphidz(const int *I,const double* x) const{ return 0.; };
    
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
        
    const double* getX(const int &i) const{std::cout<<"Abort in tet0\n"; abort();};
    const int* getIND(const int &i) const{std::cout<<"Abort in tet0\n"; abort();};
    const int* getKVERT_IND(const int &i) const {std::cout<<"Abort in tet0\n"; abort();};
  };

  //************************************************************

  class quad_lag : public basis{
  public:
    quad_lag(const unsigned int& type, const int& nc, const int& nf):
      basis(type, nc, nf, 4, 8, 9){ };
    const double* getX(const int &i) const{return X[i];};
    const int* getIND(const int &i) const{return IND[i];};
    const int* getKVERT_IND(const int &i) const {return KVERT_IND[i];};
  protected: 
    static const double X[25][3];
    static const int IND[9][3];
    static const int KVERT_IND[25][2];
    
  };
  
  
  class quad1: public quad_lag {
  public:
    quad1(): quad_lag(8, 4, 9) {};
    void PrintType() const { std::cout<<" quad1 ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; }
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; }
    double eval_d2phidxdy(const int *I,const double* x) const;
  };

  //************************************************************

  class quadth: public quad_lag {
  public:
    quadth(): quad_lag(9, 8, 21) {};
    void PrintType() const { std::cout<<" quadth ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const;
    double eval_d2phidy2(const int *I,const double* x) const;
    double eval_d2phidxdy(const int *I,const double* x) const;
  };

  //************************************************************

  class quad2: public quad_lag {
  public:
    quad2(): quad_lag(10, 9, 25) {};
    void PrintType() const { std::cout<<" quad2 ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const;
    double eval_d2phidy2(const int *I,const double* x) const;
    double eval_d2phidxdy(const int *I,const double* x) const;
  };

  //******************************************************************************

  class quad_const : public basis{
  public:
    quad_const(const unsigned int& type, const int& nc, const int& nf):
      basis(type, nc, nf, 4, 8, 9){ };
    const double* getX(const int &i) const{return X[i];};
    const int* getIND(const int &i) const{return IND[i];};
    const int* getKVERT_IND(const int &i) const {return KVERT_IND[i];};
  
  protected: 
    static const double X[12][3];
    static const int IND[3][3];
    static const int KVERT_IND[12][2];
  };
  
  
  class quad0: public quad_const {
  public:
    quad0(): quad_const(15, 1, 4) {};
    void PrintType() const { std::cout<<" quad0 ";};
    
    double eval_phi(const int *I,const double* x) const{ return 1.; };
    double eval_dphidx(const int *I,const double* x) const{ return 0.; };
    double eval_dphidy(const int *I,const double* x) const{ return 0.; };
    
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
    
  };

  //******************************************************************************

  class quadpwl: public quad_const {
  public:
    quadpwl(): quad_const(16, 3, 12) {};
    void PrintType() const { std::cout<<" quadpwl ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
  };

  //************************************************************

  class tri_lag : public basis{
  public: 
    tri_lag(const unsigned int& type, const int& nc, const int& nf):
      basis(type, nc, nf, 3, 6, 6){ };
    const double* getX(const int &i) const{return X[i];};
    const int* getIND(const int &i) const{return IND[i];};
    const int* getKVERT_IND(const int &i) const {return KVERT_IND[i];};
    
  protected:
    static const double X[15][3];
    static const int IND[6][3];
    static const int KVERT_IND[15][2];
  };
  
  
  class tri1: public tri_lag {
  public:
    tri1(): tri_lag(11, 3, 6) {}; 
    void PrintType() const { std::cout<<" tri1 ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
  };

  //************************************************************

  class tri2: public tri_lag {
  public:
    tri2(): tri_lag(12, 6, 15) {}; 
    void PrintType() const { std::cout<<" tri2 ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
  
    double eval_d2phidx2(const int *I,const double* x) const;
    double eval_d2phidy2(const int *I,const double* x) const;
    double eval_d2phidxdy(const int *I,const double* x) const;
  
  };

  //************************************************************

  class tri0: public basis {
  public:
    tri0(): basis(40, 1, 4, 3, 6, 6){ };    
    void PrintType() const { std::cout<<" tri0 ";};
        
    double eval_phi(const int *I,const double* x) const{ return 1.; };
    double eval_dphidx(const int *I,const double* x) const{ return 0.; };
    double eval_dphidy(const int *I,const double* x) const{ return 0.; };
    
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
        
    const double* getX(const int &i) const{std::cout<<"Abort in tri0\n"; abort();};
    const int* getIND(const int &i) const{std::cout<<"Abort in tri0\n"; abort();};
    const int* getKVERT_IND(const int &i) const {std::cout<<"Abort in tri0\n"; abort();};
  };

  //************************************************************

  class line_lag : public basis{ 
  public:
    line_lag(const unsigned int& type, const int& nc, const int& nf):
      basis(type, nc, nf, 2, 3, 3){ };
    const double* getX(const int &i) const{return X[i];};
    const int* getIND(const int &i) const{return IND[i];};
    const int* getKVERT_IND(const int &i) const {return KVERT_IND[i];};
    
  protected:
    static const double X[5][3];
    static const int IND[3][3];
    static const int KVERT_IND[5][2];
  };
  
  class line1: public line_lag {
  public:
    line1(): line_lag(13, 2, 3) {}; 
    void PrintType() const { std::cout<<" line1 ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
  
  };

  //************************************************************
  class line2: public line_lag {
  public:
    line2(): line_lag(14, 3, 5) {}; 
    void PrintType() const { std::cout<<" line2 ";};
    
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_d2phidx2(const int *I,const double* x) const;
  };

  //************************************************************

  class line0: public basis {
  public:
    line0(): basis(50, 1, 2, 2, 3, 3){ };  
    void PrintType() const { std::cout<<" line0 ";};
    
    double eval_phi(const int *I,const double* x) const{ return 1.; };
    double eval_dphidx(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    
    const double* getX(const int &i) const{std::cout<<"Abort in line0\n"; abort();};
    const int* getIND(const int &i) const{std::cout<<"Abort in line0\n"; abort();};
    const int* getKVERT_IND(const int &i) const {std::cout<<"Abort in line0\n"; abort();};
    
  };

}//end namespace femus


#endif

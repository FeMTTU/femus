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

#include <iostream>
#include <stdlib.h>  

namespace femus {
    
  class basis {
  public:
    virtual void PrintType() const = 0 ;
    virtual double eval_phi(const int *I,const double* x) const{ return 1.; };
    virtual double eval_dphidx(const int *I,const double* x) const{ return 0.; };
    virtual double eval_dphidy(const int *I,const double* x) const{ return 0.; };
    virtual double eval_dphidz(const int *I,const double* x) const{ return 0.; };
       
    virtual double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    virtual double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    virtual double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    virtual double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
    virtual double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    virtual double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
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


  class hex1: public basis {
  public:
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

  class hexth: public basis {
  
  public:
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

  class hex2: public basis {
  private:
  
  public:
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

  class hex0: public basis {
  public:
    void PrintType() const { std::cout<<" hex0 ";};
  };

  //******************************************************************************

  class hexpwl: public basis {
  public:
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


  class wedge1: public basis {
  public:
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

  class wedgeth: public basis {
  public:
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

  class wedge2: public basis {
  
  public:
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
    void PrintType() const { std::cout<<" wedge0 ";};
  };
  
  
  //************************************************************

  class tet1: public basis {
  public:
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

  class tet2: public basis {

  public:
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
    void PrintType() const { std::cout<<" tet0 ";};
  };

  //************************************************************

  class quad1: public basis {
  public:
    void PrintType() const { std::cout<<" quad1 ";};
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const{ return 0.; };
  
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; }
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; }
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const;
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
  };

  //************************************************************

  class quadth: public basis {
  private:
    double th2(const double& x,  const int& i) const;
    double dth2(const double& x,  const int& i) const;
    double d2th2(const double& x,  const int& i) const;
  public:
    void PrintType() const { std::cout<<" quadth ";};
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const{ return 0.; };
  
    double eval_d2phidx2(const int *I,const double* x) const;
    double eval_d2phidy2(const int *I,const double* x) const;
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const;
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
  };

  //************************************************************

  class quad2: public basis {
  public:
    void PrintType() const { std::cout<<" quad2 ";};
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const{ return 0.; };
  
    double eval_d2phidx2(const int *I,const double* x) const;
    double eval_d2phidy2(const int *I,const double* x) const;
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const;
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
  
  };

  //******************************************************************************

  class quad0: public basis {
  public:
    void PrintType() const { std::cout<<" quad0 ";};
  };

  //******************************************************************************

  class quadpwl: public basis {
  public:
    void PrintType() const { std::cout<<" quadpwl ";};
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const{ return 0.; };
  
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
  };

  //************************************************************

  class tri1: public basis {
  public:
    void PrintType() const { std::cout<<" tri1 ";};
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
    double eval_dphidz(const int *I,const double* x) const{ return 0.; };
  
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
  };

  //************************************************************

  class tri2: public basis {
  public:
    void PrintType() const { std::cout<<" tri2 ";};
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const;
  
    double eval_dphidz(const int *I,const double* x) const{ return 0.; };
  
    double eval_d2phidx2(const int *I,const double* x) const;
    double eval_d2phidy2(const int *I,const double* x) const;
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const;
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
  
  };

  //************************************************************

  class tri0: public basis {
  public:
    void PrintType() const { std::cout<<" tri0 ";};
  };

  //************************************************************

  class line1: public basis {
  public:
    void PrintType() const { std::cout<<" line1 ";};
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const{ return 0.; };
    double eval_dphidz(const int *I,const double* x) const{ return 0.; };
  
    double eval_d2phidx2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
  
  };

  //************************************************************
  class line2: public basis {
  public:
    void PrintType() const { std::cout<<" line2 ";};
    double eval_phi(const int *I,const double* x) const;
    double eval_dphidx(const int *I,const double* x) const;
    double eval_dphidy(const int *I,const double* x) const{ return 0.; };
    double eval_dphidz(const int *I,const double* x) const{ return 0.; };
  
    double eval_d2phidx2(const int *I,const double* x) const;
    double eval_d2phidy2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidz2(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidxdy(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidydz(const int *I,const double* x) const{ return 0.; };
    double eval_d2phidzdx(const int *I,const double* x) const{ return 0.; };
  };

  //************************************************************

  class line0: public basis {
  public:
    void PrintType() const { std::cout<<" line0 ";};
  };

} //end namespace femus


#endif

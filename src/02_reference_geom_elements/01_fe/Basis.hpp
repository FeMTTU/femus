/*=========================================================================

  Program: FEMUS
  Module: basis
  Authors: Eugenio Aulisa, Sara Calandrini and Giorgio Bornia

  Copyright (c) FEMTTU
  All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.  See the above copyright notice for more information.

  =========================================================================*/

/**
 * This class contains the fe basis function and their derivatives
 */


#ifndef __femus_fe_Basis_hpp__
#define __femus_fe_Basis_hpp__

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>



namespace femus {
    
    

  class basis {
      
    public:

      /**
       * _nc = number of dofs of 1 element;  
       * _nf = number of dofs in that element after refinement; 
        _nlag[0] = number of linear dofs in 1 element;
        _nlag[1] = number of serendipity dofs in 1 element; 
        _nlag[2] = number of tensor-product quadratic dofs in 1 element; 
        _nlag[3] = number of tensor-product quadratic dofs in that element after 1 refinement; 
      */
      const int _nc, _nf, _nlag0, _nlag1, _nlag2, _nlag3;
      
      /** this is only needed in 3d to handle faces of different types (wedges, pyramides, ...) */
      int faceNumber[3];

      basis(const int &nc, const int &nf, const int &nlag0, const int &nlag1, const int &nlag2, const int &nlag3,
	    const int  &faceNumber0,const int  &faceNumber1, const int  &faceNumber2):
        _nc(nc),
        _nf(nf),
        _nlag0(nlag0),
        _nlag1(nlag1),
        _nlag2(nlag2),
        _nlag3(nlag3)
        {
	  faceNumber[0] = faceNumber0;
	  faceNumber[1] = faceNumber1;
	  faceNumber[2] = faceNumber2;
	  }

      virtual void PrintType() const = 0;
      
      double eval_phi(const unsigned &j, const std::vector < double > &x) const {
        return eval_phi(this->GetIND(j), &x[0]);
      }
      
      double eval_dphidx(const unsigned &j, const std::vector < double > &x) const {
        return eval_dphidx(this->GetIND(j), &x[0]);
      }
      
      double eval_dphidy(const unsigned &j, const std::vector < double > &x) const {
        return eval_dphidy(this->GetIND(j), &x[0]);
      }
      
      double eval_dphidz(const unsigned &j, const std::vector < double > &x) const {
        return eval_dphidz(this->GetIND(j), &x[0]);
      }
      
      /** Evaluate the derivatives either in the x, y or z direction **/
      double eval_dphidxyz(const unsigned int dim, const int* j, const double* x) const {
          
        assert(dim < 3); //0, 1, 2

        switch(dim) {
            case(0): { return eval_dphidx(j, x); }
            case(1): { return eval_dphidy(j, x); }
            case(2): { return eval_dphidz(j, x); }
            default: {std::cout << "Only up to dim 3" << std::endl; abort(); }
        }
        
    }
      
      double eval_d2phidx2(const unsigned &j, const std::vector < double > &x) const {
        return eval_d2phidx2(this->GetIND(j), &x[0]);
      }
      
      double eval_d2phidy2(const unsigned &j, const std::vector < double > &x) const {
        return eval_d2phidy2(this->GetIND(j), &x[0]);
      }
      
      double eval_d2phidz2(const unsigned &j, const std::vector < double > &x) const {
        return eval_d2phidz2(this->GetIND(j), &x[0]);
      }
      
      double eval_d2phidxdy(const unsigned &j, const std::vector < double > &x) const {
        return eval_d2phidxdy(this->GetIND(j), &x[0]);
      }
      
      double eval_d2phidydz(const unsigned &j, const std::vector < double > &x) const {
        return eval_d2phidydz(this->GetIND(j), &x[0]);
      }
      
      double eval_d2phidzdx(const unsigned &j, const std::vector < double > &x) const {
        return eval_d2phidzdx(this->GetIND(j), &x[0]);
      }
      
      
      
      virtual double eval_phi(const int *I, const double* x) const {
        std::cout << "Error this phi is not available for this element \n";
        abort();
      }
      
      virtual double eval_dphidx(const int *I, const double* x) const {
        std::cout << "Error this dphidx is not available for this element dimension\n";
        abort();
      }
      
      virtual double eval_dphidy(const int *I, const double* x) const {
        std::cout << "Error this dphidy is not available for this element dimension\n";
        abort();
      }
      
      virtual double eval_dphidz(const int *I, const double* x) const {
        std::cout << "Error this dphidz is not available for this element dimension\n";
        abort();
      }
      
      virtual double eval_d2phidx2(const int *I, const double* x) const {
        std::cout << "Error this d2phix2 is not available for this element dimension\n";
        abort();
      }
      
      virtual double eval_d2phidy2(const int *I, const double* x) const {
        std::cout << "Error this d2phidy2 is not available for this element dimension\n";
        abort();
      }
      
      virtual double eval_d2phidz2(const int *I, const double* x) const {
        std::cout << "Error this d2phidz2 is not available for this element dimension\n";
        abort();
      }
      
      virtual double eval_d2phidxdy(const int *I, const double* x) const {
        std::cout << "Error this d2phidxdy is not available for this element dimension\n";
        abort();
      }
      
      virtual double eval_d2phidydz(const int *I, const double* x) const {
        std::cout << "Error this d2phidydz is not available for this element dimension\n";
        abort();
      }
      
      virtual double eval_d2phidzdx(const int *I, const double* x) const {
        std::cout << "Error this d2phidydz is not available for this element dimension\n";
        abort();
      }

      virtual const double* GetX(const int &i) const = 0;

      virtual const double* GetXcoarse(const int &i) const {
        std::cout << "Warning this function in not yet implemented for this element type" << std::endl;
        return NULL;
      }

      virtual void SetX(const unsigned &i, const unsigned &j, const double &value) {
        std::cout << "Warning this function in not yet implemented for this element type" << std::endl;
      }

      virtual const int* GetIND(const int &i) const = 0;

      virtual const int* GetKVERT_IND(const int &i) const = 0;

      virtual const unsigned GetFine2CoarseVertexMapping(const int &i, const unsigned &j) const {
        std::cout << "Warning this function in not implemented for const element type" << std::endl;
        return 0u;
      }

      virtual const unsigned GetFaceDof(const unsigned &i, const unsigned &j) const {
	std::cout << "Warning AAA this function in not yet implemented for this element type" << std::endl;
    return 0u;
      }

    protected:
        
      //1D basis - BEGIN
      
      // linear lagrangian
      inline double lagLinear(const double& x, const int& i) const {
        return (!i) * 0.5 * (1. - x) + !(i - 2) * 0.5 * (1. + x);
      }
      
      inline double dlagLinear(const double& x, const int& i) const {
        return (!i) * (-0.5) + !(i - 2) * 0.5;
      }

      //quadratic lagrangian
      inline double lagQuadratic(const double& x, const int& i) const {
        return !i * (0.5) * (1. - x) + !(i - 1) * (1. - x) * (1. + x) + !(i - 2) * (0.5) * (1. + x);
      }

      inline double dlagQuadratic(const double& x, const int& i) const {
        return (!i) * (-0.5) + !(i - 1) * (-2.*x) + !(i - 2) * (0.5);
      }

      inline double d2lagQuadratic(const double& x, const int& i) const {
        return !(i - 1) * (-2.);
      }

      //bi-quadratic lagrangian
      inline double lagBiquadratic(const double& x, const int& i) const {
        return !i * 0.5 * x * (x - 1.) + !(i - 1) * (1. - x) * (1. + x) + !(i - 2) * 0.5 * x * (1. + x);
      }

      inline double dlagBiquadratic(const double& x, const int& i) const {
        return !i * (x - 0.5) + !(i - 1) * (-2.*x) + !(i - 2) * (x + 0.5);
      }

      inline double d2lagBiquadratic(const double& x, const int& i) const {
        return !i + !(i - 1) * (-2.) + !(i - 2);
      }
      //1D basis - END


      //2D basis - BEGIN
      
      // linear triangle
      inline double triangleLinear(const double& x, const double& y, const int& i, const int& j) const {
        return (!i * !j) * (1. - x - y) + !(i - 2) * x + !(j - 2) * y;
      }

      inline double dtriangleLineardx(const double& x, const double& y, const int& i, const int& j) const {
        return -(!i * !j) + !(i - 2);
      }

      inline double dtriangleLineardy(const double& x, const double& y, const int& i, const int& j) const {
        return -(!i * !j) + !(j - 2);
      }

      // quadratic triangle
      inline double triangleQuadratic(const double& x, const double& y, const int& i, const int& j) const {
        return
          !i     * (!j * (1. - x - y) * (1. - 2.*x - 2.*y) + !(j - 1) * 4.*y * (1. - x - y) + !(j - 2) * (-y + 2.*y * y)) +
          !(i - 1) * (!j * 4.*x * (1. - x - y) + !(j - 1) * 4.*x * y) +
          !(i - 2) * (!j * (-x + 2.*x * x));
      }

      inline double dtriangleQuadraticdx(const double& x, const double& y, const int& i, const int& j) const {
        return
          !i     * (!j * (-3. + 4.*x + 4.*y) + !(j - 1) * y * (-4.)) +
          !(i - 1) * (!j * 4.*(1. - 2.*x - y)  + !(j - 1) * y * (4.)) +
          !(i - 2) * (!j * (-1 + 4.*x));
      }

      inline double dtriangleQuadraticdy(const double& x, const double& y, const int& i, const int& j) const {
        return
          !j     * (!i * (-3. + 4.*y + 4.*x) + !(i - 1) * x * (-4.)) +
          !(j - 1) * (!i * 4.*(1. - 2.*y - x)  + !(i - 1) * x * (4.)) +
          !(j - 2) * (!i * (-1 + 4.*y));
      }

      inline double d2triangleQuadraticdx2(const double& x, const double& y, const int& i, const int& j) const {
        return !j * ((!i) * 4. + !(i - 1) * (-8.) + !(i - 2) * 4.);
      }

      inline double d2triangleQuadraticdy2(const double& x, const double& y, const int& i, const int& j) const {
        return !i * ((!j) * 4. + !(j - 1) * (-8.) + !(j - 2) * 4.);
      }

      inline double d2triangleQuadraticdxdy(const double& x, const double& y, const int& i, const int& j) const {
        return ((!i) * (!j) + !(i - 1) * !(j - 1)) * 4. + (!(i - 1) * (!j) + (!i) * !(j - 1)) * (-4.);
      }

      //biquadratic triangle

      inline double triangleBiquadratic(const double& x, const double& y, const int& i, const int& j) const {
        return
          !i     * (!j * ((1. - x - y) * (1. - 2.*x - 2.*y) + 3.*x * y * (1 - x - y)) +
                    !(j - 1) * 4.*(y * (1. - x - y) - 3.*x * y * (1 - x - y))  +
                    !(j - 2) * (-y + 2.*y * y + 3.*x * y * (1 - x - y))) +
          !(i - 1) * (!j * 4.*(x * (1. - x - y) - 3.*x * y * (1 - x - y)) +
                      !(j - 1) * 4.*(x * y - 3.*x * y * (1 - x - y))) +
          !(i - 2) * (!j * (-x + 2.*x * x + 3.*x * y * (1 - x - y))) +
          !(i - 7) * (!(j - 7) * 27.*x * y * (1 - x - y));
      }

      inline double dtriangleBiquadraticdx(const double& x, const double& y, const int& i, const int& j) const {
        return
          !i     * (!j * (-3. + 4.*x + 4.*y + 3.*(y - 2.*x * y - y * y)) +
                    !(j - 1) * 4.*(-y - 3.*(y - 2.*x * y - y * y))  +
                    !(j - 2) * 3.*(y - 2.*x * y - y * y)) +
          !(i - 1) * (!j * 4.*(1. - 2.*x - y - 3.*(y - 2.*x * y - y * y))  +
                      !(j - 1) * 4.*(y - 3.*(y - 2.*x * y - y * y))) +
          !(i - 2) * (!j * (-1 + 4.*x + 3.*(y - 2.*x * y - y * y))) +
          !(i - 7) * (!(j - 7) * 27.*(y - 2.*x * y - y * y));
      }

      inline double dtriangleBiquadraticdy(const double& x, const double& y, const int& i, const int& j) const {
        return
          !j     * (!i * (-3. + 4.*y + 4.*x + 3.*(x - x * x - 2.*x * y)) +
                    !(i - 1) * 4.*(-x - 3.*(x - x * x - 2.*x * y))  +
                    !(i - 2) * 3.*(x - x * x - 2.*x * y)) +
          !(j - 1) * (!i * 4.*(1. - 2.*y - x - 3.*(x - x * x - 2.*x * y))  +
                      !(i - 1) * 4.*(x - 3.*(x - x * x - 2.*x * y))) +
          !(j - 2) * (!i * (-1 + 4.*y + 3.*(x - x * x - 2.*x * y))) +
          !(j - 7) * (!(i - 7) * 27.*(x - x * x - 2.*x * y));
      }

      inline double d2triangleBiquadraticdx2(const double& x, const double& y, const int& i, const int& j) const {
        return
          !i     * (!j * (4. - 6.*y) +
                    !(j - 1) * 4.*(6.*y)  +
                    !(j - 2) * (-6.*y)) +
          !(i - 1) * (!j * 4.*(-2. + 6.*y)  +
                      !(j - 1) * 4.*(6.*y)) +
          !(i - 2) * (!j * (4. - 6.*y)) +
          !(i - 7) * (!(j - 7) * (-54.*y));
      }

      inline double d2triangleBiquadraticdy2(const double& x, const double& y, const int& i, const int& j) const {
        return
          !j     * (!i * (4. - 6.*x) +
                    !(i - 1) * 4.*(6.*x) +
                    !(i - 2) * (-6.*x)) +
          !(j - 1) * (!i * 4.*(-2. + 6.*x)  +
                      !(i - 1) * 4.*(6.*x)) +
          !(j - 2) * (!i * (4. - 6.*x)) +
          !(j - 7) * (!(i - 7) * (-54.*x));
      }

      inline double d2triangleBiquadraticdxdy(const double& x, const double& y, const int& i, const int& j) const {
        return
          !j     * (!i * (4. + 3.*(1 - 2.*x - 2.*y)) +
                    !(i - 1) * 4.*(-1. - 3.*(1 - 2.*x - 2.*y))  +
                    !(i - 2) * 3.*(1 - 2.*x - 2.*y)) +
          !(j - 1) * (!i * 4.*(-1. - 3.*(1 - 2.*x - 2.*y))  +
                      !(i - 1) * 4.*(1. - 3.*(1 - 2.*x - 2.*y))) +
          !(j - 2) * (!i * (3.*(1 - 2.*x - 2.*y))) +
          !(j - 7) * (!(i - 7) * 27.*(1 - 2.*x - 2.*y));
      }
      //2D basis - END
      

  };



}//end namespace femus



#endif

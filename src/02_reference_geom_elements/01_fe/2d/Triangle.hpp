#ifndef __femus_fe_Triangle_hpp__
#define __femus_fe_Triangle_hpp__



#include "Basis.hpp"

#include <iostream>
#include <vector>


//maximum is Tri7
//3, 6, 7
#define LAGRANGE_TRIANGLE_NDOFS_MAXIMUM_FINE  19
#define DISCPOLY_TRIANGLE_NDOFS_MAXIMUM_FINE  12


namespace femus {
    

  //******** TRIANGLE - BEGIN ****************************************************
  
  
  //******** C0 LAGRANGE - BEGIN ****************************************************
  

  class tri_lag : public basis {
      
    public:
        
      tri_lag(const int& nc, const int& nf):
        basis(nc, nf, LAGRANGE_TRIANGLE_NDOFS_MAXIMUM_FINE, 0, 3, 3) { }
        
      const double* GetX(const int &i) const {
        return X[i];
      }
      
      void SetX(const unsigned &i, const unsigned &j, const double &value) {
        X[i][j] = value;
      }
      
      const double* GetXcoarse(const int &i) const {
        return Xc[i];
      }
      
      const int* GetIND(const int &i) const {
        return IND[i];
      }
      
      const int* GetKVERT_IND(const int &i) const {
        return KVERT_IND[i];
      }

      const unsigned GetFine2CoarseVertexMapping(const int &i, const unsigned &j)  const {
        return fine2CoarseVertexMapping[i][j];
      }
      
      const unsigned GetFaceDof(const unsigned &i, const unsigned &j) const {
        return faceDofs[i][j];
      }
      
      
 // === FE Shape Functions, Service for other elements (bit of repetition) - 2D basis, triangle BEGIN ============
      
    public:
      
      // linear triangle
      static inline double triangleLinear(const double& x, const double& y, const int& i, const int& j) {
        return (!i * !j) * (1. - x - y) + !(i - 2) * x + !(j - 2) * y;
      }

      static inline double dtriangleLineardx(const double& x, const double& y, const int& i, const int& j) {
        return -(!i * !j) + !(i - 2);
      }

      static inline double dtriangleLineardy(const double& x, const double& y, const int& i, const int& j) {
        return -(!i * !j) + !(j - 2);
      }

      // quadratic triangle
      static inline double triangleQuadratic(const double& x, const double& y, const int& i, const int& j) {
        return
          !i     * (!j * (1. - x - y) * (1. - 2.*x - 2.*y) + !(j - 1) * 4.*y * (1. - x - y) + !(j - 2) * (-y + 2.*y * y)) +
          !(i - 1) * (!j * 4.*x * (1. - x - y) + !(j - 1) * 4.*x * y) +
          !(i - 2) * (!j * (-x + 2.*x * x));
      }

      static inline double dtriangleQuadraticdx(const double& x, const double& y, const int& i, const int& j) {
        return
          !i     * (!j * (-3. + 4.*x + 4.*y) + !(j - 1) * y * (-4.)) +
          !(i - 1) * (!j * 4.*(1. - 2.*x - y)  + !(j - 1) * y * (4.)) +
          !(i - 2) * (!j * (-1 + 4.*x));
      }

      static inline double dtriangleQuadraticdy(const double& x, const double& y, const int& i, const int& j) {
        return
          !j     * (!i * (-3. + 4.*y + 4.*x) + !(i - 1) * x * (-4.)) +
          !(j - 1) * (!i * 4.*(1. - 2.*y - x)  + !(i - 1) * x * (4.)) +
          !(j - 2) * (!i * (-1 + 4.*y));
      }

      static inline double d2triangleQuadraticdx2(const double& x, const double& y, const int& i, const int& j) {
        return !j * ((!i) * 4. + !(i - 1) * (-8.) + !(i - 2) * 4.);
      }

      static inline double d2triangleQuadraticdy2(const double& x, const double& y, const int& i, const int& j) {
        return !i * ((!j) * 4. + !(j - 1) * (-8.) + !(j - 2) * 4.);
      }

      static inline double d2triangleQuadraticdxdy(const double& x, const double& y, const int& i, const int& j) {
        return ((!i) * (!j) + !(i - 1) * !(j - 1)) * 4. + (!(i - 1) * (!j) + (!i) * !(j - 1)) * (-4.);
      }

      //biquadratic triangle

      static inline double triangleBiquadratic(const double& x, const double& y, const int& i, const int& j) {
        return
          !i     * (!j * ((1. - x - y) * (1. - 2.*x - 2.*y) + 3.*x * y * (1 - x - y)) +
                    !(j - 1) * 4.*(y * (1. - x - y) - 3.*x * y * (1 - x - y))  +
                    !(j - 2) * (-y + 2.*y * y + 3.*x * y * (1 - x - y))) +
          !(i - 1) * (!j * 4.*(x * (1. - x - y) - 3.*x * y * (1 - x - y)) +
                      !(j - 1) * 4.*(x * y - 3.*x * y * (1 - x - y))) +
          !(i - 2) * (!j * (-x + 2.*x * x + 3.*x * y * (1 - x - y))) +
          !(i - 7) * (!(j - 7) * 27.*x * y * (1 - x - y));
      }

      static inline double dtriangleBiquadraticdx(const double& x, const double& y, const int& i, const int& j) {
        return
          !i     * (!j * (-3. + 4.*x + 4.*y + 3.*(y - 2.*x * y - y * y)) +
                    !(j - 1) * 4.*(-y - 3.*(y - 2.*x * y - y * y))  +
                    !(j - 2) * 3.*(y - 2.*x * y - y * y)) +
          !(i - 1) * (!j * 4.*(1. - 2.*x - y - 3.*(y - 2.*x * y - y * y))  +
                      !(j - 1) * 4.*(y - 3.*(y - 2.*x * y - y * y))) +
          !(i - 2) * (!j * (-1 + 4.*x + 3.*(y - 2.*x * y - y * y))) +
          !(i - 7) * (!(j - 7) * 27.*(y - 2.*x * y - y * y));
      }

      static inline double dtriangleBiquadraticdy(const double& x, const double& y, const int& i, const int& j) {
        return
          !j     * (!i * (-3. + 4.*y + 4.*x + 3.*(x - x * x - 2.*x * y)) +
                    !(i - 1) * 4.*(-x - 3.*(x - x * x - 2.*x * y))  +
                    !(i - 2) * 3.*(x - x * x - 2.*x * y)) +
          !(j - 1) * (!i * 4.*(1. - 2.*y - x - 3.*(x - x * x - 2.*x * y))  +
                      !(i - 1) * 4.*(x - 3.*(x - x * x - 2.*x * y))) +
          !(j - 2) * (!i * (-1 + 4.*y + 3.*(x - x * x - 2.*x * y))) +
          !(j - 7) * (!(i - 7) * 27.*(x - x * x - 2.*x * y));
      }

      static inline double d2triangleBiquadraticdx2(const double& x, const double& y, const int& i, const int& j) {
        return
          !i     * (!j * (4. - 6.*y) +
                    !(j - 1) * 4.*(6.*y)  +
                    !(j - 2) * (-6.*y)) +
          !(i - 1) * (!j * 4.*(-2. + 6.*y)  +
                      !(j - 1) * 4.*(6.*y)) +
          !(i - 2) * (!j * (4. - 6.*y)) +
          !(i - 7) * (!(j - 7) * (-54.*y));
      }

      static inline double d2triangleBiquadraticdy2(const double& x, const double& y, const int& i, const int& j) {
        return
          !j     * (!i * (4. - 6.*x) +
                    !(i - 1) * 4.*(6.*x) +
                    !(i - 2) * (-6.*x)) +
          !(j - 1) * (!i * 4.*(-2. + 6.*x)  +
                      !(i - 1) * 4.*(6.*x)) +
          !(j - 2) * (!i * (4. - 6.*x)) +
          !(j - 7) * (!(i - 7) * (-54.*x));
      }

      static inline double d2triangleBiquadraticdxdy(const double& x, const double& y, const int& i, const int& j) {
        return
          !j     * (!i * (4. + 3.*(1 - 2.*x - 2.*y)) +
                    !(i - 1) * 4.*(-1. - 3.*(1 - 2.*x - 2.*y))  +
                    !(i - 2) * 3.*(1 - 2.*x - 2.*y)) +
          !(j - 1) * (!i * 4.*(-1. - 3.*(1 - 2.*x - 2.*y))  +
                      !(i - 1) * 4.*(1. - 3.*(1 - 2.*x - 2.*y))) +
          !(j - 2) * (!i * (3.*(1 - 2.*x - 2.*y))) +
          !(j - 7) * (!(i - 7) * 27.*(1 - 2.*x - 2.*y));
      }
      
 // === FE Shape Functions, Service for other elements (bit of repetition) - 2D basis, triangle END ============
      
      

    private:
        
      static const double Xc[7][2];
      
      static const int IND[7][2];
      
      double X[ LAGRANGE_TRIANGLE_NDOFS_MAXIMUM_FINE ][2];
      static const int KVERT_IND[ LAGRANGE_TRIANGLE_NDOFS_MAXIMUM_FINE ][2];

      static const unsigned fine2CoarseVertexMapping[4][3];
      static const unsigned faceDofs[3][3];

  };


  class TriLinear: public tri_lag {
      
    public:
        
      TriLinear(): tri_lag(3, 6) {}
      
      void PrintType() const {
        std::cout << " TriLinear ";
      }

      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_dphidy(const int *I, const double* x) const;

      double eval_d2phidx2(const int *I, const double* x) const {
        return 0.;
      }
      double eval_d2phidy2(const int *I, const double* x) const {
        return 0.;
      }
      double eval_d2phidxdy(const int *I, const double* x) const {
        return 0.;
      }
      
  };

  //************************************************************

  class TriQuadratic: public tri_lag {
      
    public:
        
      TriQuadratic(): tri_lag(6, 15) {}
      
      void PrintType() const {
        std::cout << " TriQuadratic ";
      }

      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_dphidy(const int *I, const double* x) const;

      double eval_d2phidx2(const int *I, const double* x) const;
      double eval_d2phidy2(const int *I, const double* x) const;
      double eval_d2phidxdy(const int *I, const double* x) const;

  };

  //************************************************************

  class TriBiquadratic: public tri_lag {
      
    public:
        
      TriBiquadratic(): tri_lag(7, 19) {}
      
      void PrintType() const {
        std::cout << " TriBiquadratic ";
      }

      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_dphidy(const int *I, const double* x) const;

      double eval_d2phidx2(const int *I, const double* x) const;
      double eval_d2phidy2(const int *I, const double* x) const;
      double eval_d2phidxdy(const int *I, const double* x) const;

  };

  //******** C0 LAGRANGE - END ****************************************************


  //******** DISCONTINUOUS POLYNOMIAL - BEGIN ****************************************************


  class tri_const : public basis {
      
    public:
        
      tri_const(const int& nc, const int& nf):
        basis(nc, nf, DISCPOLY_TRIANGLE_NDOFS_MAXIMUM_FINE, 0, 3, 3) { }
        
      const double* GetX(const int &i) const {
        return X[i];
      }
      
      const int* GetIND(const int &i) const {
        return IND[i];
      }
      
      const int* GetKVERT_IND(const int &i) const {
        return KVERT_IND[i];
      }

    private:
        
      static const int IND[3][2];
      
      static const double X[ DISCPOLY_TRIANGLE_NDOFS_MAXIMUM_FINE ][2];
      static const int KVERT_IND[  DISCPOLY_TRIANGLE_NDOFS_MAXIMUM_FINE ][2];

  };

  class tri0: public tri_const {
      
    public:
        
      tri0(): tri_const(1, 4) { }
      
      void PrintType() const {
        std::cout << " tri0 ";
      }

      double eval_phi(const int *I, const double* x) const {
        return 1.;
      }
      
      double eval_dphidx(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_dphidy(const int *I, const double* x) const {
        return 0.;
      }

      double eval_d2phidx2(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidy2(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidxdy(const int *I, const double* x) const {
        return 0.;
      }
      
  };

  class tripwLinear: public tri_const {
      
    public:
        
      tripwLinear(): tri_const(3, 12) { }
      
      void PrintType() const {
        std::cout << " tripwLinear ";
      }

      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_dphidy(const int *I, const double* x) const;

      double eval_d2phidx2(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidy2(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidxdy(const int *I, const double* x) const {
        return 0.;
      }
      
  };
  
  
  //******** DISCONTINUOUS POLYNOMIAL - END ****************************************************

  
  
  //******** TRIANGLE - END ****************************************************
  

}

#endif
 

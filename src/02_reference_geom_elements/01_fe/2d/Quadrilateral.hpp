#ifndef __femus_fe_Quadri_hpp__
#define __femus_fe_Quadri_hpp__


#include "Basis.hpp"

#include <iostream>
#include <vector>


namespace femus {
    

  //******** QUADRILATERAL - BEGIN ****************************************************

  
  //******** C0 LAGRANGE - BEGIN ****************************************************
  
  
  class quad_lag : public basis {
      
    public:
        
      quad_lag(const int& nc, const int& nf):
        basis(nc, nf, 4, 8, 9, 25, 0, 4, 4) { }
        
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

    protected:
        
      static const double Xc[9][2];
      
      double X[25][2];
      static const int IND[9][2];
      static const int KVERT_IND[25][2];

      static const unsigned fine2CoarseVertexMapping[4][4];
      static const unsigned faceDofs[4][3];

  };


  class QuadLinear: public quad_lag {
      
    public:
        
      QuadLinear(): quad_lag(4, 9) {}
      
      void PrintType() const {
        std::cout << " QuadLinear ";
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
      double eval_d2phidxdy(const int *I, const double* x) const;
  };

  //************************************************************

  class QuadQuadratic: public quad_lag {
      
    public:
        
      QuadQuadratic(): quad_lag(8, 21) {}
      
      void PrintType() const {
        std::cout << " QuadQuadratic ";
      }

      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_dphidy(const int *I, const double* x) const;

      double eval_d2phidx2(const int *I, const double* x) const;
      double eval_d2phidy2(const int *I, const double* x) const;
      double eval_d2phidxdy(const int *I, const double* x) const;
  };

  //************************************************************

  class QuadBiquadratic: public quad_lag {
      
    public:
        
      QuadBiquadratic(): quad_lag(9, 25) {}
      
      void PrintType() const {
        std::cout << " QuadBiquadratic ";
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

  class quad_const : public basis {
      
    public:
        
      quad_const(const int& nc, const int& nf):
        basis(nc, nf, 4, 8, 9, 25,  0, 4, 4) { }
        
      const double* GetX(const int &i) const {
        return X[i];
      }
      const int* GetIND(const int &i) const {
        return IND[i];
      }
      const int* GetKVERT_IND(const int &i) const {
        return KVERT_IND[i];
      }


    protected:
        
      static const double X[12][2];
      static const int IND[3][2];
      static const int KVERT_IND[12][2];


  };


  class quad0: public quad_const {
      
    public:
        
      quad0(): quad_const(1, 4) {}
      
      void PrintType() const {
        std::cout << " quad0 ";
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

  //******************************************************************************

  class quadpwLinear: public quad_const {
      
    public:
        
      quadpwLinear(): quad_const(3, 12) {}
      
      void PrintType() const {
        std::cout << " quadpwLinear ";
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

  
  //******** QUADRILATERAL - END ****************************************************


}


#endif


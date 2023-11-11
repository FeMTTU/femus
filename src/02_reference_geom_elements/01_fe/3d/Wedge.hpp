#ifndef __femus_fe_Wedge_hpp__
#define __femus_fe_Wedge_hpp__



#include "Basis.hpp"

#include <iostream>
#include <vector>


namespace femus {
    

  //******** WEDGE - BEGIN ****************************************************

  
  //******** C0 LAGRANGE - BEGIN ****************************************************
  
  
  class wedge_lag : public basis {
      
    public:
        
      wedge_lag(const int& nc, const int& nf):
        basis(nc, nf, 6, 15, 21, 95, 0, 3, 5) { }

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
        
      static const double Xc[21][3];
      double X[95][3];
      static const int IND[21][3];
      static const int KVERT_IND[95][2];

      static const unsigned fine2CoarseVertexMapping[8][8];
      static const unsigned faceDofs[5][9];

  };
  
  

  class WedgeLinear: public wedge_lag {
      
    public:
        
      WedgeLinear(): wedge_lag(6, 18) {}
      
      void PrintType() const {
        std::cout << " WedgeLinear ";
      }

      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_dphidy(const int *I, const double* x) const;
      double eval_dphidz(const int *I, const double* x) const;

      double eval_d2phidx2(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidy2(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidz2(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidxdy(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidydz(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidzdx(const int *I, const double* x) const {
        return 0.;
      }

  };

  //************************************************************

  class WedgeQuadratic: public wedge_lag {
      
    public:
        
      WedgeQuadratic(): wedge_lag(15, 57) {};
      
      void PrintType() const {
        std::cout << " WedgeQuadratic ";
      }

      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_dphidy(const int *I, const double* x) const;
      double eval_dphidz(const int *I, const double* x) const;

      double eval_d2phidx2(const int *I, const double* x) const;
      double eval_d2phidy2(const int *I, const double* x) const;
      double eval_d2phidz2(const int *I, const double* x) const;
      double eval_d2phidxdy(const int *I, const double* x) const;
      double eval_d2phidydz(const int *I, const double* x) const;
      double eval_d2phidzdx(const int *I, const double* x) const;

  };
  

  //************************************************************

  class WedgeBiquadratic: public wedge_lag {
    public:
      WedgeBiquadratic(): wedge_lag(21, 95) {}
      
      void PrintType() const {
        std::cout << " WedgeBiquadratic ";
      }

      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_dphidy(const int *I, const double* x) const;
      double eval_dphidz(const int *I, const double* x) const;

      double eval_d2phidx2(const int *I, const double* x) const;
      double eval_d2phidy2(const int *I, const double* x) const;
      double eval_d2phidz2(const int *I, const double* x) const;
      double eval_d2phidxdy(const int *I, const double* x) const;
      double eval_d2phidydz(const int *I, const double* x) const;
      double eval_d2phidzdx(const int *I, const double* x) const;

  };
  
  //******** C0 LAGRANGE - END ****************************************************

  
  //******** DISCONTINUOUS POLYNOMIAL - BEGIN ****************************************************
  
  
  class wedge_const : public basis {
      
    public:
        
      wedge_const(const int& nc, const int& nf):
        basis(nc, nf, 6, 15, 21, 95, 0, 3, 5) { }

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
        
      static const double X[32][3];
      static const int IND[4][3];
      static const int KVERT_IND[32][2];

  };

  
  class wedge0: public wedge_const {
      
  public:
      
      wedge0(): wedge_const(1, 8) { }
      
      void PrintType() const {
        std::cout << " wedge0 ";
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
      
      double eval_dphidz(const int *I, const double* x) const {
        return 0.;
      }

      double eval_d2phidx2(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidy2(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidz2(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidxdy(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidydz(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidzdx(const int *I, const double* x) const {
        return 0.;
      }

  };
  

  class wedgepwLinear: public wedge_const {
      
    public:
        
      wedgepwLinear(): wedge_const(4, 32) { }
      
      void PrintType() const {
        std::cout << " wedgepwLinear ";
      };

      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_dphidy(const int *I, const double* x) const;
      double eval_dphidz(const int *I, const double* x) const;

      double eval_d2phidx2(const int *I, const double* x) const {
        return 0.;
      }
      double eval_d2phidy2(const int *I, const double* x) const {
        return 0.;
      }
      double eval_d2phidz2(const int *I, const double* x) const {
        return 0.;
      }
      double eval_d2phidxdy(const int *I, const double* x) const {
        return 0.;
      }
      double eval_d2phidydz(const int *I, const double* x) const {
        return 0.;
      }
      double eval_d2phidzdx(const int *I, const double* x) const {
        return 0.;
      }

  };
  //******** DISCONTINUOUS POLYNOMIAL - END ****************************************************

  //******** WEDGE - END ****************************************************


}

#endif
 
 

#ifndef __femus_fe_Triangle_hpp__
#define __femus_fe_Triangle_hpp__



#include "Basis.hpp"

#include <iostream>
#include <vector>


namespace femus {
    

  //******** TRIANGLE - BEGIN ****************************************************

  class tri_lag : public basis {
      
    public:
        
      tri_lag(const int& nc, const int& nf):
        basis(nc, nf, 3, 6, 7, 19, 0, 3, 3) { }
        
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
        
      static const double Xc[7][2];
      double X[19][2];
      static const int IND[7][2];
      static const int KVERT_IND[19][2];

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

  //************************************************************

  class tri_const : public basis {
      
    public:
        
      tri_const(const int& nc, const int& nf):
        basis(nc, nf, 3, 6, 7, 19, 0, 3, 3) { }
        
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

  //******** TRIANGLE - END ****************************************************
  

}

#endif
 

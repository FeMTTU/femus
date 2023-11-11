#ifndef __femus_fe_Line_hpp__
#define __femus_fe_Line_hpp__


#include "Basis.hpp"

#include <iostream>
#include <vector>


namespace femus {
    

  
  //******** EDGE - BEGIN ****************************************************
  
  //******** C0 LAGRANGE - BEGIN ****************************************************
  

  class line_lag : public basis {
      
    public:
        
      line_lag(const int& nc, const int& nf):
        basis(nc, nf, 2, 3, 3, 5, 0, 2, 2) { }
        
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

    protected:
        
      static const double Xc[3][1];
      double X[5][2];  ///@todo why does this have 2 in the second component instead of 1?
      static const int IND[3][1];
      static const int KVERT_IND[5][2];
      
      static const unsigned fine2CoarseVertexMapping[2][2];
  };

  class LineLinear: public line_lag {
      
    public:
        
      LineLinear(): line_lag(2, 3) {}
      
      void PrintType() const {
        std::cout << " LineLinear ";
      }

      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_d2phidx2(const int *I, const double* x) const {
        return 0.;
      }

  };

  //************************************************************
  class LineBiquadratic: public line_lag {
      
    public:
        
      LineBiquadratic(): line_lag(3, 5) {}
      
      void PrintType() const {
        std::cout << " LineBiquadratic ";
      }

      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_d2phidx2(const int *I, const double* x) const;
      
  };
  
  
  //******** C0 LAGRANGE - END ****************************************************

  
  
  //******** DISCONTINUOUS POLYNOMIAL - BEGIN ****************************************************

  class line_const : public basis {
      
    public:
        
      line_const(const int& nc, const int& nf):
        basis(nc, nf, 2, 3, 3, 5, 0, 2, 2) { }
        
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
        
      static const double X[4][1];
      static const int IND[2][1];
      static const int KVERT_IND[4][2];
      
  };


  class line0: public line_const {
      
    public:
        
      line0(): line_const(1, 2) { }
      
      void PrintType() const {
        std::cout << " line0 ";
      }

      double eval_phi(const int *I, const double* x) const {
        return 1.;
      }
      
      double eval_dphidx(const int *I, const double* x) const {
        return 0.;
      }
      
      double eval_d2phidx2(const int *I, const double* x) const {
        return 0.;
      }
      

  };


  class linepwLinear: public line_const {
      
    public:
        
      linepwLinear(): line_const(2, 4) { }
      
      void PrintType() const {
        std::cout << " linepwLinear ";
      }

      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_d2phidx2(const int *I, const double* x) const {
        return 0.;
      }

  };

  //******** DISCONTINUOUS POLYNOMIAL - END ****************************************************

  //******** EDGE - END ****************************************************

  
}

#endif

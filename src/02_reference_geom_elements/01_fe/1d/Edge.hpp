#ifndef __femus_fe_Line_hpp__
#define __femus_fe_Line_hpp__


#include "Basis.hpp"
#include "GeomElemEdge.hpp"

#include <iostream>
#include <vector>

//maximum is Edge3
//2, 3, 3,
#define LAGRANGE_EDGE_NDOFS_MAXIMUM_FINE  5
#define DISCPOLY_EDGE_NDOFS_MAXIMUM_FINE  4



namespace femus {
    

  
  //******** EDGE - BEGIN ****************************************************
  
  //******** C0 LAGRANGE - BEGIN ****************************************************
  

  class line_lag : public basis {
      
        
// ===  Constructors / Destructor - BEGIN =================
    public:
      
      line_lag(const int& nc, const int& nf) :
        basis(nc, nf, LAGRANGE_EDGE_NDOFS_MAXIMUM_FINE, new GeomElemEdge()  ) { }   ///@todo consider a smart pointer here
// ===  Constructors / Destructor - END =================


 // ===  FE - BEGIN =================
   public:
     
      const double* GetXcoarse(const int &i) const {
        return Xc[i];
      }
         
    private:
        
      /// DofCarrier of coarse dofs
      static const double Xc[3][1];
      
 // ===  FE - END =================
 
 
 // ===  FE, Shape functions - BEGIN =================
   public:
         
      const int* GetIND(const int &i) const {
        return IND[i];
      }
      
    private:
      
      static const int IND[3][1];
 // ===  FE, Shape functions - END =================

 
 // === FE, Shape Functions, Service for other elements (bit of repetition) - 1D basis BEGIN ============
      
    public:  /*now they are all static and used by the functions needing them */
    
      
      // linear lagrangian
     static inline  double lagLinear(const double& x, const int& i) {
        return (!i) * 0.5 * (1. - x) + !(i - 2) * 0.5 * (1. + x);
      }
      
      static inline double dlagLinear(const double& x, const int& i)  {
        return (!i) * (-0.5) + !(i - 2) * 0.5;
      }

      //quadratic lagrangian
      static inline double lagQuadratic(const double& x, const int& i) {
        return !i * (0.5) * (1. - x) + !(i - 1) * (1. - x) * (1. + x) + !(i - 2) * (0.5) * (1. + x);
      }

      static inline double dlagQuadratic(const double& x, const int& i) {
        return (!i) * (-0.5) + !(i - 1) * (-2.*x) + !(i - 2) * (0.5);
      }

      static inline double d2lagQuadratic(const double& x, const int& i) {
        return !(i - 1) * (-2.);
      }

      //bi-quadratic lagrangian
      static inline double lagBiquadratic(const double& x, const int& i) {
        return !i * 0.5 * x * (x - 1.) + !(i - 1) * (1. - x) * (1. + x) + !(i - 2) * 0.5 * x * (1. + x);
      }

      static inline double dlagBiquadratic(const double& x, const int& i) {
        return !i * (x - 0.5) + !(i - 1) * (-2.*x) + !(i - 2) * (x + 0.5);
      }

      static inline double d2lagBiquadratic(const double& x, const int& i) {
        return !i + !(i - 1) * (-2.) + !(i - 2);
      }


 // === FE, Shape Functions, Service for other elements (bit of repetition) - 1D basis END ============
      

 // ===  FE, Refinement  - BEGIN =================
   public:
         
      const double* GetX(const int &i) const {
        return X[i];
      }
      
      void SetX(const unsigned &i, const unsigned &j, const double &value) {
        X[i][j] = value;
      }
      
      const int* GetKVERT_IND(const int &i) const {
        return KVERT_IND[i];
      }
      
    private:
      
      /// DofCarrier of fine dofs // these are static for discontinuous poly, but here they are filled afterwards
      double X[ LAGRANGE_EDGE_NDOFS_MAXIMUM_FINE ][1];
      static const int KVERT_IND[ LAGRANGE_EDGE_NDOFS_MAXIMUM_FINE ][2];
      
 // ===  FE, Refinement  - END =================
      
      
 // ===  FE, Refinement, only underlying Lagrange linear  - BEGIN =================
   public:
      
      const unsigned GetFine2CoarseVertexMapping(const int &i, const unsigned &j)  const {
        return fine2CoarseVertexMapping[i][j];
      }
 
    private:
      
      static const unsigned fine2CoarseVertexMapping[2][2];
 // ===  FE, Refinement, only underlying Lagrange linear  - END =================

      
  };

  class LineLinear: public line_lag {
      
// ===  Constructors / Destructor - BEGIN =================
    public:
        
      LineLinear(): line_lag(2, 3) {}
// ===  Constructors / Destructor - END =================


 // ===  FE - BEGIN =================
      void PrintType() const {
        std::cout << " LineLinear ";
      }
 // ===  FE - END =================
 
 // ===  FE, Shape functions - BEGIN =================
      double eval_phi(const int *I, const double* x) const;
      double eval_dphidx(const int *I, const double* x) const;
      double eval_d2phidx2(const int *I, const double* x) const {
        return 0.;
      }
 // ===  FE, Shape functions - END =================
 
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
        
// ===  Constructors / Destructor - BEGIN =================
      line_const(const int& nc, const int& nf):
        basis(nc, nf, DISCPOLY_EDGE_NDOFS_MAXIMUM_FINE, new GeomElemEdge() ) { }
// ===  Constructors / Destructor - END =================


 // ===  FE, Shape functions - BEGIN =================
    public:
      
      const int* GetIND(const int &i) const {
        return IND[i];
      }
      
      
    private:
      
      static const int IND[2][1];
 // ===  FE, Shape functions - END =================
      
 // ===  FE, Refinement  - BEGIN =================
      
    public:

      const double* GetX(const int &i) const {
        return X[i];
      }
      

      const int* GetKVERT_IND(const int &i) const {
        return KVERT_IND[i];
      }

    private:
        
      static const double X[ DISCPOLY_EDGE_NDOFS_MAXIMUM_FINE ][1];
      static const int KVERT_IND[  DISCPOLY_EDGE_NDOFS_MAXIMUM_FINE ][2];
      
 // ===  FE, Refinement  - END =================
      
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

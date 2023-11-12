#ifndef __femus_fe_Hexa_hpp__
#define __femus_fe_Hexa_hpp__



#include "Basis.hpp"

#include <iostream>
#include <vector>


namespace femus {
    
  
  //******** HEXAHEDRON - BEGIN ****************************************************

  
  
  //******** C0 LAGRANGE - BEGIN ****************************************************

  class hex_lag : public basis {
    public:
      hex_lag(const int& nc, const int& nf):
        basis(nc, nf, 8, 20, 27, LAGRANGE_HEXAHEDRON_NDOFS_MAXIMUM_FINE, 0, 6, 6) { };

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

    private:
        
      
      static const double Xc[27][3];
      
      static const int IND[27][3];
      
      double X[ LAGRANGE_HEXAHEDRON_NDOFS_MAXIMUM_FINE ][3];
      static const int KVERT_IND[ LAGRANGE_HEXAHEDRON_NDOFS_MAXIMUM_FINE ][2];

      static const unsigned fine2CoarseVertexMapping[8][8];
      static const unsigned faceDofs[6][9];

  };
  
  

  class HexLinear: public hex_lag {
      
    public:
      HexLinear(): hex_lag(8, 27) {};
      
      void PrintType() const {
        std::cout << " HexLinear ";
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
      
      double eval_d2phidxdy(const int *I, const double* x) const;
      double eval_d2phidydz(const int *I, const double* x) const;
      double eval_d2phidzdx(const int *I, const double* x) const;

};

  //************************************************************

  class HexQuadratic: public hex_lag {
      
    public:
        
      HexQuadratic(): hex_lag(20, 81) {};
      void PrintType() const {
        std::cout << " HexQuadratic ";
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

  class HexBiquadratic: public hex_lag {
      
    public:
        
      HexBiquadratic(): hex_lag(27, 125) {};
      void PrintType() const {
        std::cout << " HexBiquadratic ";
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

  class hex_const : public basis {
      
    public:
        
      hex_const(const int& nc, const int& nf):
        basis(nc, nf, 8, 20, 27, DISCPOLY_HEXAHEDRON_NDOFS_MAXIMUM_FINE, 0, 6, 6) { };

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
        
      static const int IND[4][3];
      
      static const double X[ DISCPOLY_HEXAHEDRON_NDOFS_MAXIMUM_FINE ][3];
      static const int KVERT_IND[ DISCPOLY_HEXAHEDRON_NDOFS_MAXIMUM_FINE ][2];


  };


  class hex0: public hex_const {
      
    public:
        
      hex0(): hex_const(1, 8) {}
      
      void PrintType() const {
        std::cout << " hex0 ";
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

  //******************************************************************************

  class hexpwLinear: public hex_const {
      
    public:
        
      hexpwLinear(): hex_const(4, 32) {}
      
      void PrintType() const {
        std::cout << " hexpwLinear ";
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
  //******** DISCONTINUOUS POLYNOMIAL - END ****************************************************

  //******** HEXAHEDRON - END ****************************************************
  
  

}


#endif
 

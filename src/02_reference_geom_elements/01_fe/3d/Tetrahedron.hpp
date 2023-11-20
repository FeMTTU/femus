#ifndef __femus_fe_Tetra_hpp__
#define __femus_fe_Tetra_hpp__



#include "Basis.hpp"
#include "GeomElemTet.hpp"

#include <iostream>
#include <vector>


//maximum is Tet15
// 4, 10, 15
#define LAGRANGE_TETRAHEDRON_NDOFS_MAXIMUM_FINE  67   //8 tetrahedra, ...
#define DISCPOLY_TETRAHEDRON_NDOFS_MAXIMUM_FINE  32



namespace femus {
    


  //******** TETRAHEDRON - BEGIN ****************************************************

  
  //******** C0 LAGRANGE - BEGIN ****************************************************


  class tet_lag : public basis {
      
    public:
        
      tet_lag(const int& nc, const int& nf):
        basis(nc, nf, LAGRANGE_TETRAHEDRON_NDOFS_MAXIMUM_FINE, new GeomElemTet() ) { }
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
        
      static const double Xc[15][3];
      static const int IND[15][3];
      
      double X[ LAGRANGE_TETRAHEDRON_NDOFS_MAXIMUM_FINE ][3];
      static const int KVERT_IND[ LAGRANGE_TETRAHEDRON_NDOFS_MAXIMUM_FINE ][2];

      
      static const unsigned fine2CoarseVertexMapping[8][4];
      static const unsigned faceDofs[4][7];
      

  };



  class TetLinear: public tet_lag {
      
    public:
        
      TetLinear(): tet_lag(4, 10) {}
      
      void PrintType() const {
        std::cout << " TetLinear ";
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

  class TetQuadratic: public tet_lag {
      
    public:
        
      TetQuadratic(): tet_lag(10, 35) {}
      
      void PrintType() const {
        std::cout << " TetQuadratic ";
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

  class TetBiquadratic: public tet_lag {
      
    public:
        
      TetBiquadratic(): tet_lag(15, 67) {}
      
      void PrintType() const {
        std::cout << " TetBiquadratic ";
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

  //************************************************************

  class tet_const : public basis {
      
    public:
        
      tet_const(const int& nc, const int& nf):
        basis(nc, nf, DISCPOLY_TETRAHEDRON_NDOFS_MAXIMUM_FINE, new GeomElemTet() ) { }

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
      
      static const double X[ DISCPOLY_TETRAHEDRON_NDOFS_MAXIMUM_FINE ][3];
      static const int KVERT_IND[ DISCPOLY_TETRAHEDRON_NDOFS_MAXIMUM_FINE ][2];
      
  };

  
  class tet0: public tet_const {
    public:
      tet0(): tet_const(1, 8) { }
      void PrintType() const {
        std::cout << " tet0 ";
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

  
  class tetpwLinear: public tet_const {
      
    public:
        
      tetpwLinear(): tet_const(4, 32) { }
      
      void PrintType() const {
        std::cout << " tetpwLinear ";
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

  //******** TETRAHEDRON - END ****************************************************
  

}


#endif
 

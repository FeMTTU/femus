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



#ifndef __femus_fe_Basis_hpp__
#define __femus_fe_Basis_hpp__


#include "FElemTypeEnum_list.hpp"
#include "GeomElemBase.hpp"

#include <iostream>
#include <cstdlib>
#include <vector>
#include <cassert>











namespace femus {
    


/**
 * This class contains the fe basis function and their derivatives
 */
  class basis {
      
      

// ===  Constructors / Destructor - BEGIN =================
    public:

      basis(const int &nc, const int &nf, 
            const int &nlag3,
            const GeomElemBase * geom_base_in) :
        _geom_elem(geom_base_in),
        faceNumber_offsets{ 0, geom_base_in->num_non_triangular_faces() , geom_base_in->num_non_triangular_faces() + geom_base_in->num_triangular_faces()},
        _nc(nc),
        _nf(nf),
        _n_dofs_max_fine_for_dof_carriers(nlag3)
        {  }
// ===  Constructors / Destructor - END =================

      
 // ===  Geom Elem - BEGIN =================
    public:
      
      const GeomElemBase *  get_geom_elem() const { return  _geom_elem; }
      
      const int n_faces_offset(const unsigned type_in) const { return faceNumber_offsets[type_in]; }
      
      const int n_faces_total() const { return faceNumber_offsets[2]; }
      
      static const unsigned int n_faces_types_max()  {  return  _n_face_types_max; } 
      
    private:
      
      const GeomElemBase *  _geom_elem;
      
      static constexpr const unsigned _n_face_types_max = 2;

      /** this is only needed in 3d to handle faces of different types (wedges, pyramides, ...)
       * 0 = begin non-triangular faces
       * 1 = end non-triangular faces, so that [1] - [0] gives the total of non-triangular faces
       * 1 = begin triangular faces
       * 2 = end triangular faces,   so that [2] - [1] gives the total of non-triangular faces
       */
      const unsigned int faceNumber_offsets[ _n_face_types_max + 1 ];
      
 // ===  Geom Elem - END =================
      
      
 // ===  FE  - BEGIN =================
    public:
      
      virtual void PrintType() const = 0;
      

      const int n_dofs() const { return _nc; }
       
      /// Coordinates of the points that are DofCarrier of the coarse element
      virtual const double* GetXcoarse(const int &i) const {
        std::cout << "Warning this function is not yet implemented for this element type" << std::endl;
        return NULL;
      }


    private:

      /**
       * _nc = number of dofs of 1 element;
      */
      const int _nc;
 // ===  FE  - END =================


 // ===  FE, only Lagrange  - BEGIN =================
    public:
      
      virtual const unsigned GetFaceDof(const unsigned &i, const unsigned &j) const {
	std::cout << "Warning AAA this function is not yet implemented for this element type" << std::endl;
    return 0u;
      }
      
      
      
 // ===  FE, only Lagrange  - END =================

      
 // ===  FE, Shape Functions - BEGIN =================
    public:
     
      virtual const int* GetIND(const int &i) const = 0;
      
      // DERIVATIVES of ORDER 0 - BEGIN ============
      double eval_phi(const unsigned &j, const std::vector < double > &x) const {
        return eval_phi(this->GetIND(j), &x[0]);
      }
      // DERIVATIVES of ORDER 0 - END ============
      
      
      
      // DERIVATIVES of ORDER 1 - BEGIN ============
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
      // DERIVATIVES of ORDER 1 - END ============

      
      // DERIVATIVES of ORDER 2 - BEGIN ============
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
      // DERIVATIVES of ORDER 2 - END ============
      
      
      
      // DERIVATIVES of ORDER 0 - BEGIN ============
      virtual double eval_phi(const int *I, const double* x) const {
        std::cout << "Error this phi is not available for this element \n";
        abort();
      }
      // DERIVATIVES of ORDER 0 - END ============
      
      // DERIVATIVES of ORDER 1 - BEGIN ============
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
      // DERIVATIVES of ORDER 1 - END ============
      
      
      // DERIVATIVES of ORDER 2 - BEGIN ============
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
      // DERIVATIVES of ORDER 2 - END ============
      
 // ===  FE, Shape Functions - END =================
 

 // ===  Refinement  - BEGIN =================

      
 // ===  FE, Refinement  - BEGIN =================
    public:
     
      const unsigned int Get_NDofs_maximum_to_fill_fine_coordinates_of_DofCarriers() const { return _n_dofs_max_fine_for_dof_carriers; }
      
      const int n_dofs_fine() const { return _nf; }
      
      /// Coordinates of the points that are DofCarrier of the fine element
      /** [_nf][_dim] coordinates of the _nf nodes in the refined elements ... @todo in what order? */ 
      virtual const double* GetX(const int &i) const = 0;

      /// Coordinates of the points that are DofCarrier of the fine element
      virtual void SetX(const unsigned &i, const unsigned &j, const double &value) {
        std::cout << "Warning this function in not yet implemented for this element type" << std::endl;
      }

      /** [_nf][2] For each _nf: 0 = id of the subdivision of the fine element, 1 = local id node on the subdivision of the fine element */
      virtual const int* GetKVERT_IND(const int &i) const = 0;

      
    private:
      
      /**
       * _nf = number of dofs in that element after refinement; 
      */
      const int _nf;

     /**
        _n_dofs_max_fine_for_dof_carriers = number of maximum dofs for that FE family in that element after 1 refinement; 
      */
      
      const int _n_dofs_max_fine_for_dof_carriers;
      
 // ===  FE, Refinement  - END =================

      
 // ===  FE, Refinement, only underlying Lagrange linear  - BEGIN =================
    public:
     
      virtual const unsigned GetFine2CoarseVertexMapping(const int &i, const unsigned &j) const {
        std::cout << "Warning this function is only implemented for linear Lagrange element type" << std::endl;
        return 0u;
      }


 // ===  FE, Refinement, only underlying Lagrange linear - END =================
      
 // ===  Refinement  - END =================
 
  };



}//end namespace femus



#endif

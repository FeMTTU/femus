#ifndef __femus_fe_projection_matrices_hpp__
#define __femus_fe_projection_matrices_hpp__


#include "ParallelObject.hpp"

#include "FElemTypeEnum_list.hpp"


 
namespace femus {
    
class elem_type;
class NumericVector;
class SparseMatrix;
class Mesh;


  class FE_Proj_Matrices : public ParallelObject { 
      
public:
    
    FE_Proj_Matrices();
    
    ~FE_Proj_Matrices();
    
    /**  FE: Get the projection matrix between Lagrange FEM at the same level mesh*/
    SparseMatrix* GetQitoQjProjection(const unsigned& itype, const unsigned& jtype, Mesh & mesh_in);

    
private:
    
    /** FE: Build the projection matrix between Lagrange FEM at the same level mesh*/
    void BuildQitoQjProjection(const unsigned& itype, const unsigned& jtype, const Mesh & mesh_in);
    
    /** for solution printing */
    void Get_QitoQjProjection_SparsityPatternSize_OneElement_OneFEFamily_Lagrange_Continuous(const Mesh& Mesh,
                                  const int& iel, 
                                  NumericVector* NNZ_d,
                                  NumericVector* NNZ_o,
                                  const unsigned& itype,
                                  const elem_type * elem_type_in_jtype,
                                  const unsigned ndofs_itype_in) const;

        
    /** for solution printing */
    void Build_QitoQjProjection_OneElement_OneFEFamily_Lagrange_Continuous(const Mesh& mesh,
                             const int& iel,
                             SparseMatrix* Projmat,
                             NumericVector* NNZ_d,
                             NumericVector* NNZ_o,
                             const unsigned& itype,
                             const elem_type * elem_type_in_jtype,
                             const unsigned ndofs_itype_in) const;


    /** FE: The projection matrix between Lagrange FEM at the same level mesh */
    SparseMatrix* _ProjQitoQj[NFE_FAMS_C_ZERO_LAGRANGE][NFE_FAMS_C_ZERO_LAGRANGE];


      
      
  };
    
    
    
    
    
 
 
}
 
 
#endif

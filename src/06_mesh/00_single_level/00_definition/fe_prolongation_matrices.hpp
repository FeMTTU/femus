#ifndef __femus_fe_prolongation_matrices_hpp__
#define __femus_fe_prolongation_matrices_hpp__


#include "ParallelObject.hpp"

#include "FElemTypeEnum_list.hpp"




namespace femus {



class elem_type;
class NumericVector;
class SparseMatrix;
class Mesh;


  class FE_Prolongation_Matrices : public ParallelObject { 

public:
    
    FE_Prolongation_Matrices();
    
    ~FE_Prolongation_Matrices();

    
    /** MESH: Set the coarser mesh from which this mesh is generated */
    void SetCoarseMesh( Mesh* otherCoarseMsh ){
      _coarseMsh = otherCoarseMsh;
    }

    /**  FE: Get the coarse to the fine projection matrix and use it to restrict only on coarse nodes i.e. projection*/
    SparseMatrix* GetCoarseToFineProjectionRestrictionOnCoarse(const unsigned& solType, const Mesh & mesh_in);

    /**  FE: Get the coarse to the fine projection matrix*/
    SparseMatrix* GetCoarseToFineProjection(const unsigned& solType, const Mesh & mesh_in);

private:
  
    /** Pointer to the coarser mesh from which this mesh is generated, it equals NULL if _level = 0 */
    Mesh* _coarseMsh;
    
    /** FE: Build the coarse to the fine projection matrix */
    void BuildCoarseToFineProjection(const unsigned& solType, const char el_dofs[], const Mesh & mesh_in);
    
    void Get_Prolongation_SparsityPatternSize_OneElement_OneFEFamily(const Mesh& meshf,
                                  const Mesh& meshc,
                                  const int& ielc,
                                  NumericVector* NNZ_d,
                                  NumericVector* NNZ_o,
                                  const char is_fine_or_coarse [],
                                  const elem_type * elem_type_in) const;

    void Build_Prolongation_OneElement_OneFEFamily(const Mesh& meshf,
                             const Mesh& meshc,
                             const int& ielc,
                             SparseMatrix* Projmat, 
                             const char is_fine_or_coarse [],
                             const elem_type * elem_type_in) const;

    
    
    /** FE: The coarse to the fine projection matrix */
    SparseMatrix* _ProjCoarseToFine[NFE_FAMS];
    
    
    
};




}


#endif

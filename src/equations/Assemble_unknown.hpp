#ifndef __femus_utils_Assemble_unknown_hpp__
#define __femus_utils_Assemble_unknown_hpp__

#include <vector>
#include <cmath>
#include <string>
#include "FElemTypeEnum.hpp"

#include "MultiLevelSolution.hpp"
#include "System.hpp"

namespace femus {

    
template < class real_num_mov >
 class Phi {
     
 public:
    
  Phi(const unsigned int dim) :
   _dim2( (3 * (dim - 1) + !(dim - 1)) ),
   _max_size_elem_dofs(static_cast< unsigned >(ceil(pow(3, dim)))) {
      
      phi.reserve(_max_size_elem_dofs);
      phi_x.reserve(_max_size_elem_dofs * dim);
      phi_xx.reserve(_max_size_elem_dofs * _dim2);
      
  }   
     
  std::vector < double > phi;
  std::vector < real_num_mov > phi_x;
  std::vector < real_num_mov > phi_xx;
  
  const unsigned int _dim2;        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  const unsigned int _max_size_elem_dofs;

 };
    


//this is based on the AddSolution function in MLSol
 class Unknown {
     
 public:
     
     std::string _name;
     FEFamily _fe_family;     
     FEOrder  _fe_order;
     int _time_order;
     bool _is_pde_unknown;
     
 };
 
 
 class System;
 class MultiLevelSolution;
 
 
template < class real_num >
class UnknownLocal {
     
 public:
     
     UnknownLocal() { }
     
    void initialize(const unsigned int dim, const Unknown & unknown, const MultiLevelSolution * ml_sol, System * mlPdeSys);
    
//  sort of copy constructor, in case you don't have an Unknown
    void initialize(const unsigned int dim, const UnknownLocal & unknown);
    
    void set_elem(const unsigned int iel, const Mesh * msh, const Solution * sol);

    inline const unsigned int & fe_type() const { return SolFEType; }

    inline const unsigned int & sol_index() const { return SolIndex; }

    inline const unsigned int & pde_index() const { return SolPdeIndex; }

    inline const unsigned int & num_elem_dofs() const { return Sol_n_el_dofs; }

    inline const std::vector < real_num > & elem_dofs() const { return Sol_eldofs; }
    
 private:
    
//these do not change with the element
  std::string  Solname;
  unsigned int SolFEType; 
  unsigned int SolPdeIndex;
  unsigned int SolIndex;  
//   _is_pde_unknown
// time order

  
//it will change elem by elem
  unsigned int             Sol_n_el_dofs;
  std::vector < real_num > Sol_eldofs;

 };
 
 
 
 template < class real_num >
 void UnknownLocal< real_num >::initialize(const unsigned int dim, const Unknown & unknown, const MultiLevelSolution * ml_sol, System * mlPdeSys) {
   
     
   Solname     = unknown._name;
   SolPdeIndex = mlPdeSys->GetSolPdeIndex(Solname.c_str());
   SolIndex    = ml_sol->GetIndex        (Solname.c_str());
   SolFEType   = ml_sol->GetSolutionType(SolIndex);
   
   
   const unsigned int max_size_elem_dofs = static_cast< unsigned >(ceil(pow(3, dim)));
   
   Sol_eldofs.reserve(max_size_elem_dofs);
  
  }
  
  
 template < class real_num >
 void UnknownLocal< real_num >::initialize(const unsigned int dim, const UnknownLocal & unknown) {
     
   Solname     =  unknown.Solname;    
   SolFEType   =  unknown.SolFEType; 
   SolPdeIndex =  unknown.SolPdeIndex;    
   SolIndex    =  unknown.SolIndex;
   
   
   const unsigned int max_size_elem_dofs = static_cast< unsigned >(ceil(pow(3, dim)));
   
   Sol_eldofs.reserve(max_size_elem_dofs);
     
} 
 
 
 template < class real_num >
 void UnknownLocal< real_num >::set_elem(const unsigned int iel, const Mesh * msh, const Solution * sol) {
     
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType);
       Sol_n_el_dofs = ndofs_unk;
       Sol_eldofs.resize(ndofs_unk);
    for (unsigned i = 0; i < ndofs_unk; i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType);
       Sol_eldofs[i] = (*sol->_Sol[SolIndex])(solDof);
      }
    
 }
 
 
 
    
}  //end namespace femus



#endif 
 

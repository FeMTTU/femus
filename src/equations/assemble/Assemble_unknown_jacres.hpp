#ifndef __femus_utils_Assemble_unknown_jacres_hpp__
#define __femus_utils_Assemble_unknown_jacres_hpp__


#include "Assemble_unknown.hpp"
#include "Mesh.hpp"
#include "LinearEquationSolver.hpp"



namespace femus {

    
    
    
    
    
template < typename real_num >
class ElementJacRes {
    
public:
    
    ElementJacRes(const unsigned int dim, const std::vector< UnknownLocal < real_num > > & unknowns_in) : 
    unknowns(unknowns_in),
    n_unknowns(unknowns.size()) {
        
        const unsigned int max_size_elem_dofs = static_cast< unsigned >(ceil(pow(3, dim)));
        
        const unsigned total_size_1d = n_unknowns * max_size_elem_dofs;
        Res.reserve( total_size_1d );  
        Jac.reserve( total_size_1d * total_size_1d);

        loc_to_glob_map_all_vars.reserve( total_size_1d );
        
        loc_to_glob_map.resize(n_unknowns);
        for(int i = 0; i < n_unknowns; i++) {    loc_to_glob_map[i].reserve(max_size_elem_dofs); }
        
    }
    
inline void set_loc_to_glob_map(const unsigned int iel, const Mesh * msh,  LinearEquationSolver* pdeSys) {
    
   loc_to_glob_map_all_vars.resize(0);
      for (unsigned  u = 0; u < n_unknowns; u++) {
          
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, unknowns[u].fe_type());
//        Sol_n_el_dofs[u] = ndofs_unk;                              //in unk local
//        SolVAR_eldofs[u].resize(Sol_n_el_dofs[u]);                 //in unk local 
          loc_to_glob_map[u].resize(unknowns[u].num_elem_dofs()); 
    for (unsigned i = 0; i < unknowns[u].num_elem_dofs(); i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, unknowns[u].fe_type());
//        SolVAR_eldofs[u][i] = (*sol->_Sol[SolIndex[u]])(solDof);   //in unk local
          loc_to_glob_map[u][i] = pdeSys->GetSystemDof(unknowns[u].sol_index(), unknowns[u].pde_index(), i, iel);    // global to global mapping between solution node and pdeSys dof
      }
      
   loc_to_glob_map_all_vars.insert(loc_to_glob_map_all_vars.end(),loc_to_glob_map[u].begin(),loc_to_glob_map[u].end());
   
    }
    
    unsigned sum_Sol_n_el_dofs = 0;
    for (unsigned  u = 0; u < n_unknowns; u++) { sum_Sol_n_el_dofs += unknowns[u].num_elem_dofs(); }
    
    Jac.resize(sum_Sol_n_el_dofs * sum_Sol_n_el_dofs);  std::fill(Jac.begin(), Jac.end(), 0.);
    Res.resize(sum_Sol_n_el_dofs);                      std::fill(Res.begin(), Res.end(), 0.);
    
    
  }

  
inline static unsigned int compute_max_n_dofs(const std::vector < unsigned int > & Sol_n_el_dofs) {
    
        unsigned int nDof_max    =  0;
        
        for (unsigned  k = 0; k < Sol_n_el_dofs.size(); k++)     {
            if(Sol_n_el_dofs[k] > nDof_max)    nDof_max = Sol_n_el_dofs[k];
        }
        
      return nDof_max;
}


inline static unsigned int compute_sum_n_dofs(const std::vector < unsigned int > & Sol_n_el_dofs) {
    
        unsigned int sum_Sol_n_el_dofs = 0;
        
        for (unsigned  k = 0; k < Sol_n_el_dofs.size(); k++) {
            sum_Sol_n_el_dofs += Sol_n_el_dofs[k];
        }
        
      return sum_Sol_n_el_dofs;
}

   inline vector < int >      & dof_map() { return loc_to_glob_map_all_vars; }
    
   inline vector < real_num > & res() { return Res; }
    
   inline vector < double >   & jac() { return Jac; }                         
   
   
private:
    
  vector < int >      loc_to_glob_map_all_vars;  
  vector < vector < int > >     loc_to_glob_map;  
  vector < real_num >  Res;                         
  vector < double >    Jac;
  const std::vector< UnknownLocal < real_num > > & unknowns;   
  const unsigned int n_unknowns;
    
};
    
    
    
    
    
    
    
    
    
 
    
    
    
    
    

}










#endif

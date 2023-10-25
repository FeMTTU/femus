#ifndef __femus_Assemble_useful_functions_hpp__
#define __femus_Assemble_useful_functions_hpp__


#include "MultiLevelMesh.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NonLinearImplicitSystem.hpp"

#include "Assemble_jacobian.hpp"

 
namespace femus {

//*******************************************************************************************
//*********************** Domain and Mesh Independent - BEGIN *****************************************
//*******************************************************************************************






 
void el_dofs_quantities_vol(const Solution*                sol,
                        const Mesh * msh,
                        const unsigned int iel,
                        const    std::vector < unsigned > & SolFEType,
                        std::vector < unsigned > & Sol_n_el_dofs 
     ) {
    
        for (unsigned  k = 0; k < Sol_n_el_dofs.size(); k++) {
            unsigned  ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
            Sol_n_el_dofs[k] = ndofs_unk;
        }   
    
}


 
void el_dofs_unknowns_vol(const Solution*                sol,
                      const Mesh * msh,
                      const  LinearEquationSolver* pdeSys,
                      const unsigned int iel,
                      const    std::vector < unsigned > & SolFEType,
                      const std::vector < unsigned > & SolIndex,
                      const std::vector < unsigned > & SolPdeIndex,
                      std::vector < unsigned > & Sol_n_el_dofs, 
                      std::vector < std::vector < double > > & sol_eldofs,  
                      std::vector < std::vector < int > > & L2G_dofmap ) {
    
    assert(Sol_n_el_dofs.size() == sol_eldofs.size());
    
        //all vars###################################################################
        for (unsigned  k = 0; k < Sol_n_el_dofs.size(); k++) {
            unsigned  ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
            Sol_n_el_dofs[k] = ndofs_unk;
            sol_eldofs[k].resize(ndofs_unk);
            L2G_dofmap[k].resize(ndofs_unk);
            for (unsigned i = 0; i < ndofs_unk; i++) {
                unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);
                sol_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);
                L2G_dofmap[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);
            }
        }
        //all vars###################################################################

}




void  print_global_residual_jacobian(const bool print_algebra_global,
                                     const MultiLevelProblem& ml_prob,
                                     const NonLinearImplicitSystem * mlPdeSys,
                                       LinearEquationSolver* pdeSys,
                                     NumericVector*           RES,
                                     SparseMatrix*             JAC,
                                     const unsigned iproc,
                                     const bool assembleMatrix)  {


  if (print_algebra_global) {
    const unsigned nonlin_iter = mlPdeSys->GetNonlinearIt();
    assemble_jacobian< double, double >::print_global_jacobian(assembleMatrix, ml_prob, JAC, nonlin_iter);
//     assemble_jacobian< double, double >::print_global_residual(ml_prob, RES,  mlPdeSys->GetNonlinearIt());

    RES->close();
    std::ostringstream res_out; res_out << ml_prob.GetFilesHandler()->GetOutputPath() << "./" << "res_" << nonlin_iter  << ".txt";
    pdeSys->print_with_structure_matlab_friendly(iproc, res_out.str().c_str(), RES);

     } 
  
  }


  
    
  
  bool check_if_same_elem(const unsigned iel, const unsigned jel) {
      
   return (iel == jel);
      
  }
  
  
  
  bool check_if_same_elem_bdry(const unsigned iel, const unsigned jel, const unsigned iface, const unsigned jface) {

   return (iel == jel && iface == jface);
      
  }


  
  void  set_dense_pattern_for_unknowns(NonLinearImplicitSystem/*WithPrimalDualActiveSetMethod*/  & system, const std::vector < Unknown > unknowns)  {
  
///     system.init();   ///@todo Understand why I cannot put this here but it has to be in the main(), there must be some objects that get destroyed, passing the reference is not enough

  const MultiLevelProblem &  ml_prob = system.GetMLProb();
  const MultiLevelMesh *  ml_mesh = ml_prob.GetMLMesh();
  
  unsigned n_levels = ml_mesh->GetNumberOfLevels();
  std::ostringstream sp_out_base; sp_out_base << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "sp_";
  system._LinSolver[n_levels - 1]->sparsity_pattern_print_nonzeros(sp_out_base.str(), "on");
  system._LinSolver[n_levels - 1]->sparsity_pattern_print_nonzeros(sp_out_base.str(), "off");

  
  const Mesh* msh = ml_mesh->GetLevel(n_levels - 1);
  unsigned nprocs = msh->n_processors();
  unsigned iproc = msh->processor_id();

  
    for(int ivar = 0; ivar < unknowns.size(); ivar++) {
  
        if (unknowns[ivar]._is_sparse == false ) { 
        
  unsigned variable_index = system.GetSolPdeIndex(unknowns[ivar]._name.c_str());
  
  
  
  unsigned n_dofs_var_all_procs = 0;
  for (int ip = 0; ip < nprocs; ip++) {
     n_dofs_var_all_procs += system._LinSolver[n_levels - 1]->KKoffset[variable_index + 1][ip] - system._LinSolver[n_levels - 1]->KKoffset[variable_index][ip];
  // how does this depend on the number of levels and the number of processors? 
  // For the processors I summed over them and it seems to work fine
  // For the levels... should I pick the coarsest level instead of the finest one, or is it the same?
   } 

  unsigned n_vars = system.GetSolPdeIndex().size();   //assume all variables are dense: we should sum 3 sparse + 1 dense... or better n_components_ctrl * dense and the rest sparse
  //check that the dofs are picked correctly, it doesn't seem so 
  
  system.SetSparsityPatternMinimumSize (n_dofs_var_all_procs * n_vars, unknowns[ivar]._name);  ///@todo this is like AddSolution: it increases a vector
        }
    }
    
    
    
  return; 
  
  }


//*******************************************************************************************
//*********************** Domain and Mesh Independent - END *****************************************
//*******************************************************************************************


}
 
 
#endif

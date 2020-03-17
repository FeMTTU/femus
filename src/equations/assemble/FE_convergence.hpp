#ifndef __femus_utils_FE_convergence_hpp__
#define __femus_utils_FE_convergence_hpp__


#include "MultiLevelSolution.hpp"
#include "Math.hpp"
#include "Assemble_unknown.hpp"



namespace femus {

class MultiLevelProblem;
class MultiLevelSolution;
class MultiLevelMesh;
class Files;



class Main_single_level {
    
public: 

// What this function does is to yield a MultiLevelSolution object at each level, to construct a level hierarchy to be used for error analysis 
virtual const MultiLevelSolution  run_on_single_level(const Files & files,
                                                      MultiLevelProblem & ml_prob,
                                                      const std::vector< Unknown > & unknowns,
                                                      const MultiLevelSolution::BoundaryFuncMLProb  SetBoundaryCondition,
                                                      const MultiLevelSolution::InitFuncMLProb SetInitialCondition,
                                                      MultiLevelMesh & ml_mesh,
                                                      const unsigned lev) const = 0;
                                                     
};

template < class type = double >
class FE_convergence {
 
    
    
public: 
 

  void  convergence_study(const Files & files,
                          MultiLevelProblem & ml_prob,
                          const std::vector< Unknown > & unknowns,
                          const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition,
                          const MultiLevelSolution::InitFuncMLProb SetInitialCondition,
                          MultiLevelMesh & ml_mesh,
                          MultiLevelMesh & ml_mesh_all_levels,
                          const unsigned max_number_of_meshes,
                          const unsigned norm_flag,
                          const unsigned conv_order_flag,
                          const Main_single_level & main_in,
                          const Math::Function< double > * exact_sol = NULL);

    
    
static   std::vector < std::vector < std::vector < type > > >  initialize_vector_of_norms(const unsigned unknowns_size, 
                                                                                             const unsigned max_number_of_meshes, 
                                                                                             const unsigned norm_flag);

    
   
static   const MultiLevelSolution  initialize_convergence_study(const MultiLevelProblem & ml_prob,
                                                                const std::vector< Unknown > &  unknowns,  
                                                                MultiLevelMesh & ml_mesh_all_levels, 
                                                                const unsigned max_number_of_meshes, 
                                                                const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition,
                                                                const MultiLevelSolution::InitFuncMLProb SetInitialCondition_in
                                                               );  


//   print the error and the order of convergence between different levels
static  void output_convergence_order(const std::vector < std::vector < std::vector < type > > > &  norm,
                                    const unsigned int u,
                                    const unsigned int i,
                                    const unsigned int n);

static  void output_convergence_order_all(const std::vector< Unknown > &  unknowns,
                                        const std::vector < std::vector < std::vector < type > > > &  norms, 
                                        const unsigned norm_flag, 
                                        const unsigned max_number_of_meshes);

 
 

static  std::vector< type > compute_error_norms(const std::vector < std::vector < /*const*/ elem_type_templ_base<type, double> *  > > & elem_all,
                                                const std::vector<Gauss> & quad_rules,
                                                const MultiLevelSolution* ml_sol, 
                                                const MultiLevelSolution* ml_sol_all_levels,
                                                const std::string & unknown,
                                                const unsigned current_level,
                                                const unsigned norm_flag,
                                                const unsigned conv_order_flag,
                                                const Math::Function< type > * ex_sol_in = NULL
                                             );
 
 

     
static  void compute_error_norms_per_unknown_per_level(const std::vector < std::vector < /*const*/ elem_type_templ_base<type, double> *  > > & elem_all,
                                                       const std::vector<Gauss> & quad_rules,
                                                       const MultiLevelSolution* ml_sol_single_level,
                                                       MultiLevelSolution* ml_sol_all_levels,
                                                       const std::vector< Unknown > &  unknowns,
                                                       const unsigned i,
                                                       const unsigned norm_flag,
                                                       std::vector < std::vector < std::vector < type > > > &  norms,
                                                       const unsigned conv_order_flag,
                                                       const Math::Function< type > * ex_sol_in = NULL
                                         );


 
    
};






    

template < class type>
  void  FE_convergence< type >::convergence_study(const Files & files,
                                                  MultiLevelProblem & ml_prob,
                                                  const std::vector< Unknown > & unknowns,
                                                  const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition,
                                                  const MultiLevelSolution::InitFuncMLProb SetInitialCondition,
                                                  MultiLevelMesh & ml_mesh,
                                                  MultiLevelMesh & ml_mesh_all_levels,
                                                  const unsigned max_number_of_meshes,
                                                  const unsigned norm_flag,
                                                  const unsigned conv_order_flag,
                                                  const Main_single_level & main_in,
                                                  const Math::Function< double > * exact_sol) {


  
    
    vector < vector < vector < double > > > norms = FE_convergence::initialize_vector_of_norms ( unknowns.size(), 
                                                                                                 max_number_of_meshes, 
                                                                                                 norm_flag);
    
     MultiLevelSolution         ml_sol_all_levels = FE_convergence::initialize_convergence_study(ml_prob,
                                                                                                 unknowns,
                                                                                                 ml_mesh_all_levels, 
                                                                                                 max_number_of_meshes, 
                                                                                                 SetBoundaryCondition, 
                                                                                                 SetInitialCondition);
    
  //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<type, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
            
       for (int lev = 0; lev < max_number_of_meshes; lev++) {
                  
            const MultiLevelSolution ml_sol_single_level = main_in.run_on_single_level(files,
                                                                                       ml_prob,
                                                                                       unknowns,
                                                                                       SetBoundaryCondition,
                                                                                       SetInitialCondition, 
                                                                                       ml_mesh,
                                                                                       lev);
            

            FE_convergence::compute_error_norms_per_unknown_per_level ( elem_all,
                                                                        ml_prob.GetQuadratureRuleAllGeomElems(),
                                                                        & ml_sol_single_level,
                                                                        & ml_sol_all_levels,
                                                                        unknowns,
                                                                        lev,
                                                                        norm_flag,
                                                                        norms,
                                                                        conv_order_flag,
                                                                        exact_sol);
        
      }
   
       FE_convergence::output_convergence_order_all(unknowns, norms, norm_flag, max_number_of_meshes);
   
}

    
    
template < class type>
/*static*/   std::vector < std::vector < std::vector < type > > >  FE_convergence< type >::initialize_vector_of_norms(const unsigned unknowns_size, 
                                                                                             const unsigned max_number_of_meshes, 
                                                                                             const unsigned norm_flag) {
   
       //how many Unknowns, how many mesh levels, how many norms
       
   std::vector < std::vector < std::vector < type > > > norms( unknowns_size );
  
     for (unsigned int u = 0; u < unknowns_size; u++) {
              norms[u].resize( max_number_of_meshes );
       for (int i = 0; i < norms[u].size(); i++) {   // loop on the mesh level
               norms[u][i].resize(norm_flag + 1);
           }   
       }
    
    return norms; 
       
}

    
   
template < class type>
/*static*/   const MultiLevelSolution  FE_convergence< type >::initialize_convergence_study(const MultiLevelProblem & ml_prob, 
                                                                                            const std::vector< Unknown > &  unknowns,
                                                                                            MultiLevelMesh & ml_mesh_all_levels,
                                                                                            const unsigned max_number_of_meshes,
                                                                                            const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition_in,                                                                const MultiLevelSolution::InitFuncMLProb SetInitialCondition_in)  {

   //Mesh: construct all levels  ==================
        unsigned numberOfUniformLevels_finest = max_number_of_meshes;
        ml_mesh_all_levels.RefineMesh(numberOfUniformLevels_finest, numberOfUniformLevels_finest, NULL);
//      ml_mesh_all_levels.EraseCoarseLevels(numberOfUniformLevels - 2);  // need to keep at least two levels to send u_(i-1) projected(prolongated) into next refinement

 
 //Solution ==================
//         std::vector < MultiLevelSolution * >   ml_sol_all_levels(unknowns.size());
//                ml_sol_all_levels[u] = new MultiLevelSolution (& ml_mesh_all_levels);  //with the declaration outside and a "new" inside it persists outside the loop scopes
               MultiLevelSolution ml_sol_all_levels(& ml_mesh_all_levels);

               
               ml_sol_all_levels.AttachSetBoundaryConditionFunction(SetBoundaryCondition_in);
               for (unsigned int u = 0; u < unknowns.size(); u++) {
               ml_sol_all_levels.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order);  //We have to do so to avoid buildup of AddSolution with different FE families
               ml_sol_all_levels.Initialize(unknowns[u]._name.c_str(), SetInitialCondition_in, & ml_prob);
               ml_sol_all_levels.GenerateBdc(unknowns[u]._name.c_str(), "Steady", & ml_prob);
               }
 //Solution  ==================

 return ml_sol_all_levels;
} 
    


//   print the error and the order of convergence between different levels
template < class type>
/*static*/  void FE_convergence< type >::output_convergence_order(const std::vector < std::vector < std::vector < type > > > &  norm,
                                    const unsigned int u,
                                    const unsigned int i,
                                    const unsigned int n) {

   if(i < norm[u].size() - 2)  {
//   std::cout << norm_name << " ERROR and ORDER OF CONVERGENCE: " << fam << " " << ord << "\n\n";

    std::cout << i + 1 << "\t\t";
    std::cout.precision(14);

    std::cout << norm[u][i][n] << "\t";

    std::cout << std::endl;
  
      std::cout.precision(3);
      std::cout << "\t\t";

        std::cout << log( norm[u][i][n] / norm[u][i + 1][n] ) / log(2.) << "  \t\t\t\t";

      std::cout << std::endl;
    }
    
  
  
}


template < class type>
/*static */ void FE_convergence< type >::output_convergence_order_all(const std::vector< Unknown > &  unknowns,
                                        const std::vector < std::vector < std::vector < type > > > &  norms, 
                                        const unsigned norm_flag, 
                                        const unsigned max_number_of_meshes) {
    
    assert( unknowns.size() == norms.size() );
    
    const std::vector< std::string > norm_names = {"L2-NORM","H1-SEMINORM"};
  
     for (unsigned int u = 0; u < unknowns.size(); u++) {
       for (int n = 0; n < norm_flag + 1; n++) {
            std::cout << unknowns[u]._name << " : " << norm_names[n] << " ERROR and ORDER OF CONVERGENCE"  << std::endl;
         for (int i = 0; i < max_number_of_meshes; i++) {
                output_convergence_order(norms,u,i,n);
            }
            std::cout << std::endl;
         }
      }
      
 }
     


 
 

template < class type>
/*static*/  std::vector< type > FE_convergence< type >::compute_error_norms(const std::vector < std::vector < /*const*/ elem_type_templ_base<type, double> *  > > & elem_all,
                                                                            const std::vector<Gauss> & quad_rules,
                                                                            const MultiLevelSolution* ml_sol,
                                                                            const MultiLevelSolution* ml_sol_all_levels,
                                                                            const std::string & unknown,
                                                                            const unsigned current_level,
                                                                            const unsigned norm_flag,
                                                                            const unsigned conv_order_flag,
                                                                            const Math::Function< type > * ex_sol_in
                                             ) {
     
  // (//0 = only L2: //1 = L2 + H1)
  
// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv 
//i.e. ||prol_(u_(i-1)) - u_(i)|| = err(i) => err(i-1)/err(i) = 2^alpha ,implemented as log(err(i)/err(i+1))/log2

   if (conv_order_flag == 1 && ex_sol_in == NULL) { std::cout << "Please provide analytical solution" << std::endl; abort(); }
   
    
  const unsigned num_norms = norm_flag + 1;
  //norms that we are computing here //first L2, then H1 ============
  std::vector< type > norms(num_norms);                  std::fill(norms.begin(), norms.end(), 0.);   
  std::vector< type > norms_exact_dofs(num_norms);       std::fill(norms_exact_dofs.begin(), norms_exact_dofs.end(), 0.);
  std::vector< type > norms_inexact_dofs(num_norms);     std::fill(norms_inexact_dofs.begin(), norms_inexact_dofs.end(), 0.);
  //norms that we are computing here //first L2, then H1 ============
  
  
  unsigned level = ml_sol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = ml_sol->_mlMesh->GetLevel(level);
  const Solution* sol = ml_sol->GetSolutionLevel(level);

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

 unsigned iproc = msh->processor_id();

  //solution variable
  unsigned soluIndex = ml_sol->GetIndex(unknown.c_str()); // ml_sol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  
 //***************************************************  
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = 3;

 //***************************************************  
   std::vector < std::vector < double > >  JacI_qp(space_dim);
   std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < dim; d++) { Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < space_dim; d++) { JacI_qp[d].resize(dim); }
    
    double detJac_qp;
 //***************************************************  
   
  

  // ======================================
  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  
  
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  vector < vector < type > > x(dim_offset_grad);    // local coordinates
  for (unsigned i = 0; i < x.size(); i++)   x[i].reserve(maxSize);

//-----------------  
  vector < type > phi_coords;
  vector < type > phi_coords_x;
  vector < type > phi_coords_xx;

  phi_coords.reserve(maxSize);
  phi_coords_x.reserve(maxSize * dim_offset_grad);
  phi_coords_xx.reserve(maxSize * dim2);

  vector < type > phi;
  vector < type > phi_x;
  vector < type > phi_xx;
  
  type weight; // gauss point weight


  vector < type >  solu;                               solu.reserve(maxSize);
  vector < type >  solu_exact_at_dofs;   solu_exact_at_dofs.reserve(maxSize);
  vector < type >  solu_coarser_prol;     solu_coarser_prol.reserve(maxSize);


  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim_offset_grad);
  phi_xx.reserve(maxSize * dim2);

  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);

    // resize local arrays
    solu.resize(nDofu);
    solu_exact_at_dofs.resize(nDofu);
    solu_coarser_prol.resize(nDofu);

    for (int i = 0; i < x.size(); i++) {
      x[i].resize(nDofx);
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim_offset_grad/*dim*/; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);  // global extraction and local storage for the element coordinates
      }
    }
    
    

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
        std::vector< type > x_at_node(dim_offset_grad,0.);
        for (unsigned jdim = 0; jdim < x_at_node.size(); jdim++) x_at_node[jdim] = x[jdim][i];
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);
                   solu[i]  =                                        (*sol->_Sol[soluIndex])(solDof);
      solu_coarser_prol[i]  = (*ml_sol_all_levels->GetSolutionLevel(current_level)->_Sol[soluIndex])(solDof);
      if (ex_sol_in != NULL) solu_exact_at_dofs[i] = ex_sol_in->value(x_at_node);
    }


    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < quad_rules[ielGeom].GetGaussPointsNumber(); ig++) {

      // *** get gauss point weight, test function and test function partial derivatives ***
	elem_all[ielGeom][xType]->JacJacInv(x, ig, Jac_qp, JacI_qp, detJac_qp, space_dim);
    elem_all[ielGeom][soluType]->shape_funcs_current_elem(ig, JacI_qp, phi, phi_x, phi_xx, space_dim);
    weight = detJac_qp * quad_rules[ielGeom].GetGaussWeightsPointer()[ig];

    
      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      type solu_gss = 0.;
      type exactSol_from_dofs_gss = 0.;
      type solu_coarser_prol_gss = 0.;
      vector < type > gradSolu_gss(dim_offset_grad, 0.);
      vector < type > gradSolu_exact_at_dofs_gss(dim_offset_grad, 0.);
      vector < type > gradSolu_coarser_prol_gss(dim_offset_grad, 0.);
      vector < type > x_gss(dim_offset_grad, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss                += phi[i] * solu[i];
        exactSol_from_dofs_gss  += phi[i] * solu_exact_at_dofs[i];
        solu_coarser_prol_gss   += phi[i] * solu_coarser_prol[i];

        for (unsigned jdim = 0; jdim < dim_offset_grad; jdim++) {
          gradSolu_gss[jdim]                += phi_x[i * dim_offset_grad + jdim] * solu[i];
          gradSolu_exact_at_dofs_gss[jdim]  += phi_x[i * dim_offset_grad + jdim] * solu_exact_at_dofs[i];
          gradSolu_coarser_prol_gss[jdim]   += phi_x[i * dim_offset_grad + jdim] * solu_coarser_prol[i];
        }
        
         for (unsigned jdim = 0; jdim < dim_offset_grad; jdim++) {
           x_gss[jdim] += x[jdim][i] * phi_coords[i];
          }

      }

// H^0 ==============      
//     if (norm_flag == 0) {
      type exactSol = 0.; if (ex_sol_in != NULL) exactSol = ex_sol_in->value(x_gss);
      norms[0]               += (solu_gss - exactSol)                * (solu_gss - exactSol)       * weight;
      norms_exact_dofs[0]    += (solu_gss - exactSol_from_dofs_gss)  * (solu_gss - exactSol_from_dofs_gss) * weight;
      norms_inexact_dofs[0]  += (solu_gss - solu_coarser_prol_gss)   * (solu_gss - solu_coarser_prol_gss)  * weight;
//     }
    
// H^1 ==============      
    /*else*/ if (norm_flag == 1) {
      vector < type > exactGradSol(dim_offset_grad,0.);    if (ex_sol_in != NULL) exactGradSol = ex_sol_in->gradient(x_gss);

      for (unsigned j = 0; j < dim_offset_grad ; j++) {
        norms[1]               += ((gradSolu_gss[j] - exactGradSol[j])               * (gradSolu_gss[j]  - exactGradSol[j])) * weight;
        norms_exact_dofs[1]    += ((gradSolu_gss[j] - gradSolu_exact_at_dofs_gss[j]) * (gradSolu_gss[j] - gradSolu_exact_at_dofs_gss[j])) * weight;
        norms_inexact_dofs[1]  += ((gradSolu_gss[j] - gradSolu_coarser_prol_gss[j])  * (gradSolu_gss[j] - gradSolu_coarser_prol_gss[j]))  * weight;
      }
   }
      
   } // end gauss point loop
   
  } //end element loop for each process

  // add the norms of all processes
  NumericVector* norm_vec;
                 norm_vec = NumericVector::build().release();
                 norm_vec->init(msh->n_processors(), 1 , false, AUTOMATIC);

         /*if (norm_flag == 0) {*/ norm_vec->set(iproc, norms[0]);  norm_vec->close();  norms[0] = norm_vec->l1_norm(); /*}*/
    /*else*/ if (norm_flag == 1) { norm_vec->set(iproc, norms[1]);  norm_vec->close();  norms[1] = norm_vec->l1_norm(); }

          delete norm_vec;

   // add the norms of all processes
  NumericVector* norm_vec_exact_dofs;
                 norm_vec_exact_dofs = NumericVector::build().release();
                 norm_vec_exact_dofs->init(msh->n_processors(), 1 , false, AUTOMATIC);

         /*if (norm_flag == 0) {*/ norm_vec_exact_dofs->set(iproc, norms_exact_dofs[0]);  norm_vec_exact_dofs->close();  norms_exact_dofs[0] = norm_vec_exact_dofs->l1_norm(); /*}*/
    /*else*/ if (norm_flag == 1) { norm_vec_exact_dofs->set(iproc, norms_exact_dofs[1]);  norm_vec_exact_dofs->close();  norms_exact_dofs[1] = norm_vec_exact_dofs->l1_norm(); }

          delete norm_vec_exact_dofs;

  // add the norms of all processes
  NumericVector* norm_vec_inexact;
                 norm_vec_inexact = NumericVector::build().release();
                 norm_vec_inexact->init(msh->n_processors(), 1 , false, AUTOMATIC);

         /*if (norm_flag == 0) {*/ norm_vec_inexact->set(iproc, norms_inexact_dofs[0]);  norm_vec_inexact->close();  norms_inexact_dofs[0] = norm_vec_inexact->l1_norm(); /*}*/
    /*else*/ if (norm_flag == 1) { norm_vec_inexact->set(iproc, norms_inexact_dofs[1]);  norm_vec_inexact->close();  norms_inexact_dofs[1] = norm_vec_inexact->l1_norm(); }

          delete norm_vec_inexact;

          
    for (int n = 0; n < norms.size(); n++) { 
                  norms[n] = sqrt(norms[n]);                                            
       norms_exact_dofs[n] = sqrt(norms_exact_dofs[n]);                                            
     norms_inexact_dofs[n] = sqrt(norms_inexact_dofs[n]);
    }
    
    
if (conv_order_flag == 0)  return norms_inexact_dofs;
if (conv_order_flag == 1)  return norms;
 //   return norms_exact_dofs;

 
} 
 
 
 

     
template < class type>
/*static*/  void FE_convergence< type >::compute_error_norms_per_unknown_per_level(const std::vector < std::vector < /*const*/ elem_type_templ_base<type, double> *  > > & elem_all,
                                                                                   const std::vector<Gauss> & quad_rules,
                                                                                   const MultiLevelSolution* ml_sol_single_level,
                                                                                   MultiLevelSolution* ml_sol_all_levels,
                                                                                   const std::vector< Unknown > &  unknowns,
                                                                                   const unsigned i,
                                                                                   const unsigned norm_flag,
                                                                                   std::vector < std::vector < std::vector < type > > > &  norms,
                                                                                   const unsigned conv_order_flag,
                                                                                   const Math::Function< type > * ex_sol_in
                                         ) {
     
     
        if ( i > 0 ) {

            // ======= prolongate to the current level (i) from the coarser level (i-1) (so that you can compare the two) ========================
            ml_sol_all_levels->RefineSolution(i);
            
            // =======  compute the error norm at the current level (i) ========================
            for (unsigned int u = 0; u < unknowns.size(); u++) {  //this loop could be inside the below function
                
            const std::vector< type > norm_out = FE_convergence::compute_error_norms (elem_all, quad_rules, ml_sol_single_level, ml_sol_all_levels, unknowns[u]._name, i, norm_flag, conv_order_flag, ex_sol_in);

              for (int n = 0; n < norms[u][i-1].size(); n++)      norms[u][i-1][n] = norm_out[n];
                                       
                   }
                   
        }
                
                 
              // ======= store the last computed solution to prepare the next iteration (the current level i is now overwritten) ========================
            ml_sol_all_levels->fill_at_level_from_level(i, ml_sol_single_level->_mlMesh->GetNumberOfLevels() - 1, *ml_sol_single_level);
        
                 
                 
                 
 } 
 







}









#endif

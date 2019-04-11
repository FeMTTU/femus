#ifndef __femus_utils_Assemble_jacobian_hpp__
#define __femus_utils_Assemble_jacobian_hpp__

#include "Assemble_unknown.hpp"
#include "MultiLevelProblem.hpp"

namespace femus {
    


// class assemble_jacobian_base {
//     
//     public:
//         
//    assemble_jacobian_base();
//     
// }


class NumericVector;
class SparseMatrix;


    
template < class real_num, class real_num_mov = double >
class assemble_jacobian/* : public assemble_jacobian_base*/ {
 
    
 public:
    
                                               
 void prepare_before_integration_loop(adept::Stack& stack) const;

 
 void  compute_jacobian_inside_integration_loop(const unsigned i,
                                                const unsigned dim,
                                                const std::vector < unsigned int > Sol_n_el_dofs,
                                                const unsigned int sum_Sol_n_el_dofs,
                                                const std::vector< UnknownLocal < real_num > > & unk_vec,
                                                const std::vector< Phi < real_num_mov > > &  phi,
                                                const double weight,
                                                std::vector< double > & Jac) const;
  
                                               
 void  compute_jacobian_outside_integration_loop(adept::Stack & stack,
                                               const std::vector< std::vector< real_num > > & solu,
                                               const std::vector< real_num > & Res,
                                               std::vector< real_num_mov > & Jac, 
                                               const std::vector< int > & loc_to_glob_map,
                                               NumericVector*           RES,
                                               SparseMatrix*             KK
                                                                   ) const;
                                                                   
    
 void  compute_jacobian_outside_integration_loop(adept::Stack & stack,
                                               const std::vector< UnknownLocal < real_num > > & unk_vec,
                                               const std::vector< real_num > & Res,
                                               std::vector< real_num_mov > & Jac, 
                                               const std::vector< int > & loc_to_glob_map,
                                               NumericVector*           RES,
                                               SparseMatrix*             KK
                                                                   ) const;
                                                                   
    
static void print_element_jacobian(const unsigned int iel,
                            const std::vector < double > & Jac, 
                            const std::vector < unsigned int > & Sol_n_el_dofs, 
                            const unsigned int col_width_visualization,     
                            const unsigned int precision);


static void print_element_residual(const unsigned int iel, 
                            const std::vector < real_num > & Res, 
                            const std::vector < unsigned int > & Sol_n_el_dofs,
                            const unsigned int col_width_visualization,     
                            const unsigned int precision); 
    
static inline unsigned int res_row_index(const std::vector<unsigned int>& _Sol_n_el_dofs, const int my_row_pos, const int i);


static inline unsigned int jac_row_col_index(const std::vector<unsigned int>& _Sol_n_el_dofs, const int nDof_AllVars, const int my_row_pos, const int my_col_pos, const int i, const int j);
    




static void print_global_jacobian(const bool assemble_matrix,
                                         const MultiLevelProblem & ml_prob,
                                         const SparseMatrix* JAC,
                                         const unsigned int nonlin_iteration);
  
  
  
static void print_global_residual(const MultiLevelProblem & ml_prob,
                                  NumericVector* RES,
                                  const unsigned int nonlin_iteration);
    
static void mass_residual (std::vector < real_num > & Res,
                           const std::vector < unsigned int > & Sol_n_el_dofs,
                           const unsigned int sum_Sol_n_el_dofs,
                           const std::vector < unsigned int > & SolPdeIndex,
                           const std::vector < std::vector < double > > & phi_dof_qp,
                           const std::vector < real_num > & SolVAR_qp,
                           const double & weight_hat_qp);


 inline static double laplacian_row(const std::vector < std::vector < std::vector < double > > > & phi_x_fe_qp, 
                             const std::vector < std::vector < double > > & sol_grad_qp,
                             const int my_row_pos, 
                             const int my_col_pos, 
                             const int i,
                             const unsigned dim_in);
  
  
 inline static double laplacian_row_col(const std::vector < std::vector < std::vector < double > > > & phi_x_fe_qp,
                                        const std::vector < std::vector < std::vector < double > > > & phi_x_fe_qp_2,
                                 const int my_row_pos, 
                                 const int my_col_pos, 
                                 const int i, 
                                 const int j,
                                 const unsigned dim_in);
 
 // ===============
 
 
 


};



template < class real_num, class real_num_mov >
 inline /*static*/ double assemble_jacobian<real_num, real_num_mov>::laplacian_row(const std::vector < std::vector < std::vector < double > > > & phi_x_fe_qp,
                                                                                   const std::vector < std::vector < double > > & sol_grad_qp,
                                                                                   const int my_row_pos, 
                                                                                   const int my_col_pos, 
                                                                                   const int i,
                                                                                   const unsigned dim_in) {
  
              double laplace_rhs_i = 0.;

              for (unsigned kdim = 0; kdim < dim_in; kdim++) laplace_rhs_i +=   phi_x_fe_qp[my_row_pos][kdim][i] *
                                                                                sol_grad_qp[my_col_pos][kdim];

     return laplace_rhs_i;
}
  
  
template < class real_num, class real_num_mov >
 inline /*static*/ double assemble_jacobian<real_num, real_num_mov>::laplacian_row_col(const std::vector < std::vector < std::vector < double > > > & phi_x_fe_qp,
                                                                                       const std::vector < std::vector < std::vector < double > > > & phi_x_fe_qp_2,
                                                                                       const int my_row_pos,
                                                                                       const int my_col_pos,
                                                                                       const int i,
                                                                                       const int j,
                                                                                       const unsigned dim_in) {
  
              double laplace_mat_i_j = 0.;

              for (unsigned kdim = 0; kdim < dim_in; kdim++) laplace_mat_i_j  += phi_x_fe_qp  [my_row_pos][kdim][i] *
                                                                                 phi_x_fe_qp_2[my_col_pos][kdim][j];

     return laplace_mat_i_j;
}






template < class real_num, class real_num_mov >
/*static*/ void assemble_jacobian<real_num, real_num_mov>::mass_residual ( 
    std::vector < real_num > & Res,
    const std::vector < unsigned int > & Sol_n_el_dofs,
    const unsigned int sum_Sol_n_el_dofs,
    const std::vector < unsigned int > & SolPdeIndex,
    const std::vector < std::vector < double > > & phi_dof_qp,
    const std::vector < real_num > & SolVAR_qp,
    const double & weight_hat_qp
) {
    
 // The matrix is then filled with AD
 std::cout << "Remember to set all boundary conditions to \"dirichlet = false\"" << std::endl;
    
    
    const double test_value = 5.;
    const unsigned int n_unknowns = Sol_n_el_dofs.size();
    
                                                                                  
	for(unsigned i_unk = 0; i_unk < n_unknowns; i_unk++) { 
	    for(unsigned i_dof = 0; i_dof < Sol_n_el_dofs[i_unk]; i_dof++) {
		Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[i_unk], i_dof) ] += ( SolVAR_qp[i_unk] - test_value ) * phi_dof_qp[i_unk][i_dof] * weight_hat_qp;

// 				  for(unsigned j_unk=dim; j_unk<n_unknowns; j_unk++) {
// 		  	for(unsigned j_dof=0; j_dof < Sol_n_el_dofs[j_unk]; j_dof++) {
// 			  
// // 		              if (i_unk == j_unk )   {
// 				Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, SolPdeIndex[i_unk], SolPdeIndex[j_unk], i_unk, j_unk) ] += 
// 				        ( phi_dof_qp[i_unk][i_dof]*phi_dof_qp[j_unk][j_dof] )*weight_qp;
// // 			      }
// 			  
// 			} //j_dof
// 		  }  //j_unk
	    }  //i_dof
	}  //i_unk

 
}
 
    
template < class real_num, class real_num_mov >
/*static*/ void assemble_jacobian<real_num, real_num_mov>::print_element_jacobian(const unsigned int iel,
                            const std::vector < double > & Jac, 
                            const std::vector < unsigned int > & Sol_n_el_dofs, 
                            const unsigned int col_width_visualization,     
                            const unsigned int precision)  {

    
    const unsigned int n_unknowns = Sol_n_el_dofs.size();
    
    unsigned int nDof_AllVars = 0;
    for(unsigned i_unk = 0; i_unk < n_unknowns; i_unk++) nDof_AllVars += Sol_n_el_dofs[i_unk];
    
 std::cout << "++++++++++ iel " << iel << " ++++++++++" << std::endl;
 
    for(unsigned i_block = 0; i_block < n_unknowns; i_block++) {
      for(unsigned i_dof=0; i_dof < Sol_n_el_dofs[i_block]; i_dof++) {
	  for(unsigned j_block=0; j_block< n_unknowns; j_block++) {
               for(unsigned j_dof=0; j_dof < Sol_n_el_dofs[j_block]; j_dof++) {
 std::cout << std::right  << std::setprecision(precision) << std::setfill(' ') << std::setw(col_width_visualization) /* << std::scientific*/ << 
 Jac[ jac_row_col_index(Sol_n_el_dofs, nDof_AllVars, i_block, j_block, i_dof, j_dof) ] << " " ;
                 }
                  std::cout << " | ";
           } 
     std::cout << std::endl;
        }
        
  std::cout << std::endl;

    }

    
}


template < class real_num, class real_num_mov >
/*static*/ void assemble_jacobian<real_num, real_num_mov>::print_element_residual(const unsigned int iel, 
                            const std::vector < real_num > & Res, 
                            const std::vector < unsigned int > & Sol_n_el_dofs,
                            const unsigned int col_width_visualization,     
                            const unsigned int precision)  {
    
    const unsigned int n_unknowns = Sol_n_el_dofs.size();

 std::cout << "++++++++++ iel " << iel << " ++++++++++" << std::endl;
 
    for(unsigned i_block = 0; i_block < n_unknowns; i_block++) {
      for(unsigned i_dof=0; i_dof < Sol_n_el_dofs[i_block]; i_dof++) {
          
 std::cout << std::right  << std::setw(col_width_visualization) << std::setfill(' ') << std::setprecision(precision)  /*<< std::scientific*/ << 
 Res[  res_row_index(Sol_n_el_dofs, i_block, i_dof) ] << " " ;

          
     std::cout << std::endl;
        }
   std::cout << std::endl;

    }
    

    }
    
    
template < class real_num, class real_num_mov >
/*static*/ inline unsigned int assemble_jacobian<real_num, real_num_mov>::res_row_index(const std::vector<unsigned int>& _Sol_n_el_dofs, const int my_row_pos, const int i) {

    assert(i < _Sol_n_el_dofs[my_row_pos]); 
    
    unsigned int pos_previous = 0;
    for (unsigned k = 0; k < my_row_pos; k++) pos_previous += _Sol_n_el_dofs[k];

    return pos_previous + i;
  }
  

template < class real_num, class real_num_mov >
/*static*/  inline unsigned int assemble_jacobian<real_num, real_num_mov>::jac_row_col_index(const std::vector<unsigned int>& _Sol_n_el_dofs, const int nDof_AllVars, const int my_row_pos, const int my_col_pos, const int i, const int j) {

     assert(i < _Sol_n_el_dofs[my_row_pos]); 
     assert(j < _Sol_n_el_dofs[my_col_pos]); 
     
    unsigned int pos_previous_row = 0;
    unsigned int pos_previous_col = 0;
    for (unsigned k = 0; k < my_row_pos; k++) pos_previous_row += _Sol_n_el_dofs[k];
    for (unsigned k = 0; k < my_col_pos; k++) pos_previous_col += _Sol_n_el_dofs[k];

    return (pos_previous_row + i) * nDof_AllVars + (pos_previous_col + j);
  }
    




template < class real_num, class real_num_mov >
/*static*/ void assemble_jacobian<real_num, real_num_mov>::print_global_jacobian(const bool assemble_matrix,
                                         const MultiLevelProblem & ml_prob,
                                         const SparseMatrix* JAC,
                                         const unsigned int nonlin_iteration) {
    
  if (assemble_matrix) {
        
    JAC->close();
    std::ostringstream mat_out; mat_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "matrix_" << nonlin_iteration  << ".txt";
    JAC->print_matlab(mat_out.str(),"ascii"); //  KK->print();

   }
    
  }
  
  
  
template < class real_num, class real_num_mov >
/*static*/ void  assemble_jacobian<real_num, real_num_mov>::print_global_residual(const MultiLevelProblem & ml_prob,
                                  NumericVector* RES,
                                  const unsigned int nonlin_iteration) {
    
    RES->close();  ///@todo why does this close need to be non-const? 
    std::ostringstream res_out; res_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "res_" << nonlin_iteration  << ".txt";
    std::filebuf res_fb;
    res_fb.open (res_out.str().c_str(),std::ios::out);
    std::ostream  res_file_stream(&res_fb);
    RES->print(res_file_stream);

  }  
    
                                                                   
}  //end namespace



#endif


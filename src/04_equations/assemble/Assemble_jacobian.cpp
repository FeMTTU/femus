#include "Assemble_jacobian.hpp"
#include "System.hpp"
#include "MultiLevelSolution.hpp"



namespace femus {
    


 // template specialization for adept::adouble
template < >
 void assemble_jacobian< adept::adouble, double > ::prepare_before_integration_loop(adept::Stack & stack)  const { 
    
  stack.new_recording();    // start a new recording of all the operations involving adept variables

}   
  


 // template specialization for adept::adouble
template < > 
 void  assemble_jacobian< adept::adouble, double >::compute_jacobian_inside_integration_loop(const unsigned i,
                                                const unsigned dim,
                                                const std::vector < unsigned int > Sol_n_el_dofs,
                                                const unsigned int sum_Sol_n_el_dofs,
                                                const std::vector< UnknownLocal < adept::adouble > > & unk_vec,
                                                const std::vector< Phi < double > > &  phi,
                                                const double weight,
                                                std::vector< double > & Jac) const { }


                                               

    
 // template specialization for adept::adouble
template < >
 void  assemble_jacobian< adept::adouble, double >::compute_jacobian_outside_integration_loop (adept::Stack & stack,
                                                                    const std::vector< std::vector< adept::adouble > > & solu,
                                                                    const std::vector< adept::adouble > & Res,
                                                                    std::vector< double > & Jac,
                                                                    const std::vector< int > & loc_to_glob_map,
                                                                    NumericVector*           RES,
                                                                    SparseMatrix*             KK
                                                                   )  const {
    
    //copy the value of the adept::adoube Res in double Res and store
    //all the calculations were done with adept variables
    
 ///convert to vector of double to send to the global matrix
  std::vector < double > Res_double(Res.size());

  for (int i = 0; i < Res_double.size(); i++) {
      Res_double[i] = - Res[i].value();
    }


    stack.dependent(  & Res[0], Res_double.size());      // define the dependent variables
    for (unsigned  k = 0; k < solu.size(); k++) {   
    stack.independent(&solu[k][0],       solu[k].size());    // define the independent variables
    }
    
    stack.jacobian(&Jac[0], true);    // get the jacobian matrix (ordered by row major)

    stack.clear_independents();
    stack.clear_dependents();    
    
    RES->add_vector_blocked(Res_double, loc_to_glob_map);
    KK->add_matrix_blocked(Jac, loc_to_glob_map, loc_to_glob_map);

}
  


  
 // template specialization for adept::adouble
template < >
 void  assemble_jacobian< adept::adouble, double >::compute_jacobian_outside_integration_loop (adept::Stack & stack,
                                                                    const std::vector< UnknownLocal< adept::adouble > > & unk_loc,
                                                                    const std::vector< adept::adouble > & Res,
                                                                    std::vector< double > & Jac,
                                                                    const std::vector< int > & loc_to_glob_map,
                                                                    NumericVector*           RES,
                                                                    SparseMatrix*             KK
                                                                   )  const {
    
    //copy the value of the adept::adoube Res in double Res and store
    //all the calculations were done with adept variables
    
 ///convert to vector of double to send to the global matrix
  std::vector < double > Res_double(Res.size());

  for (int i = 0; i < Res_double.size(); i++) {
      Res_double[i] = - Res[i].value();
    }


    stack.dependent(  & Res[0], Res_double.size());      // define the dependent variables
    for (unsigned  k = 0; k < unk_loc.size(); k++) {   
    stack.independent(&((unk_loc[k].elem_dofs())[0]),       unk_loc[k].elem_dofs().size());    // define the independent variables
    }
    
    stack.jacobian(&Jac[0], true);    // get the jacobian matrix (ordered by row major)

    stack.clear_independents();
    stack.clear_dependents();    
    
    RES->add_vector_blocked(Res_double, loc_to_glob_map);
    KK->add_matrix_blocked(Jac, loc_to_glob_map, loc_to_glob_map);

}
    
    
//***************************************
//explicit instantiations
//****************************************
template class assemble_jacobian< adept::adouble, double >;




} //end namespace femus



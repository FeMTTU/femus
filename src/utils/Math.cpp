#include "Math.hpp"



namespace femus {


 // template specialization for double
template < > 
 void assemble_jacobian< double >::prepare_before_integration_loop(adept::Stack& stack) const { }


 // template specialization for adept::adouble
template < >
 void assemble_jacobian< adept::adouble  > ::prepare_before_integration_loop(adept::Stack & stack)  const { 
    
  stack.new_recording();    // start a new recording of all the operations involving adept variables

}   
  


//  // template specialization for double
// template < >
// void  assemble_jacobian< double >::compute_jacobian_inside_integration_loop(const unsigned i,
//                                                          const unsigned dim, 
//                                                          const unsigned nDofu, 
//                                                          const std::vector< double > &  phi,
//                                                          const std::vector< double > &  phi_x, 
//                                                          const double weight, 
//                                                          std::vector< double > & Jac) const { 
//     
//     std::cout << "Please implement this template specialization in your application" << std::endl; abort();
// 
// }



 // template specialization for adept::adouble
template < > 
 void  assemble_jacobian< adept::adouble >::compute_jacobian_inside_integration_loop(const unsigned i,
                                               const unsigned dim,
                                               const unsigned nDofu,
                                               const std::vector< adept::adouble > & phi,
                                               const std::vector< adept::adouble > &  phi_x, 
                                               const adept::adouble weight,
                                               std::vector< double > & Jac)  const { }

  
 
                                                   



 // template specialization for double
template < >
 void  assemble_jacobian < double > ::compute_jacobian_outside_integration_loop (adept::Stack & stack,
                                               const std::vector< double > & solu,
                                               const std::vector< double > & Res,
                                               std::vector< double > & Jac,
                                               const std::vector< int > & loc_to_glob_map,
                                               NumericVector*           RES,
                                               SparseMatrix*             KK
                                                                   )  const {
    
    RES->add_vector_blocked(Res, loc_to_glob_map);
    KK->add_matrix_blocked(Jac, loc_to_glob_map, loc_to_glob_map);
    
}

    
 // template specialization for adept::adouble
template < >
 void  assemble_jacobian< adept::adouble >::compute_jacobian_outside_integration_loop (adept::Stack & stack,
                                                                    const std::vector< adept::adouble > & solu,
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
    stack.independent(&solu[0],       solu.size());    // define the independent variables
    
    stack.jacobian(&Jac[0], true);    // get the jacobian matrix (ordered by row major)

    stack.clear_independents();
    stack.clear_dependents();    
    
    RES->add_vector_blocked(Res_double, loc_to_glob_map);
    KK->add_matrix_blocked(Jac, loc_to_glob_map, loc_to_glob_map);

}
  


    
    
 //***************************************
//explicit instantiations for double and adept::adouble
//****************************************
template class assemble_jacobian< double >; 
template class assemble_jacobian< adept::adouble >; 


}  //end namespace

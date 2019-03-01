#include "Math.hpp"



namespace femus {

namespace assemble_jacobian {



// template function: definition 
template < > 
void prepare_before_integration_loop< double >(adept::Stack& stack) { }



// template function: definition 
template < > 
void  compute_jacobian_inside_integration_loop< adept::adouble >(const unsigned i,
                                               const unsigned dim,
                                               const unsigned nDofu,
                                               const std::vector< adept::adouble > & phi,
                                               const std::vector< adept::adouble > &  phi_x, 
                                               const adept::adouble weight,
                                               std::vector< double > & Jac) { }

  
 
                                                   
 // template function: specialization for adept::adouble
template < >
void  compute_jacobian_outside_integration_loop < adept::adouble > (adept::Stack & stack,
                                                                    const std::vector< adept::adouble > & solu,
                                                                    const std::vector< adept::adouble > & Res,
                                                                    std::vector< double > & Jac,
                                                                    const std::vector< int > & loc_to_glob_map,
                                                                    NumericVector*           RES,
                                                                    SparseMatrix*             KK
                                                                   ) {
    
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


 // template function: specialization for double
template < >
void  compute_jacobian_outside_integration_loop < double > (adept::Stack & stack,
                                               const std::vector< double > & solu,
                                               const std::vector< double > & Res,
                                               std::vector< double > & Jac,
                                               const std::vector< int > & loc_to_glob_map,
                                               NumericVector*           RES,
                                               SparseMatrix*             KK
                                                                   ) {
    
    RES->add_vector_blocked(Res, loc_to_glob_map);
    KK->add_matrix_blocked(Jac, loc_to_glob_map, loc_to_glob_map);
    
}

    

 // template function: specialization for adept::adouble
template < >
void prepare_before_integration_loop< adept::adouble  > (adept::Stack & stack) { 
    
  stack.new_recording();    // start a new recording of all the operations involving adept variables

}   
    


 // template function: specialization for double
// template < >
// void  compute_jacobian_inside_integration_loop< double >(const unsigned i,
//                                                          const unsigned dim, 
//                                                          const unsigned nDofu, 
//                                                          const std::vector< double > &  phi,
//                                                          const std::vector< double > &  phi_x, 
//                                                          const double weight, 
//                                                          std::vector< double > & Jac) { 
// 
// // *** phi_j loop ***
//         for (unsigned j = 0; j < nDofu; j++) {
//           /*real_num*/double laplace_jac = 0.;
// 
//           for (unsigned kdim = 0; kdim < dim; kdim++) {
//             laplace_jac += (phi_x[i * dim + kdim] * phi_x[j * dim + kdim]);
//           }
// 
//           Jac[i * nDofu + j] += (laplace_jac + phi[i] * phi[j]) * weight;
//         } // end phi_j loop
// 
//         
// }

    
    
 //***************************************
//explicit instantiations for double and adept::adouble
//****************************************
template void prepare_before_integration_loop< double  >        (adept::Stack & stack);

template void prepare_before_integration_loop< adept::adouble  >(adept::Stack & stack);

// template void compute_jacobian_inside_integration_loop< double >(const unsigned i,
//                                                          const unsigned dim, 
//                                                          const unsigned nDofu, 
//                                                          const std::vector< double > &  phi,
//                                                          const std::vector< double > &  phi_x, 
//                                                          const double weight, 
//                                                          std::vector< double > & Jac);

template void compute_jacobian_inside_integration_loop  < adept::adouble >(const unsigned i,
                                                         const unsigned dim, 
                                                         const unsigned nDofu, 
                                                         const std::vector< adept::adouble > &  phi,
                                                         const std::vector< adept::adouble > &  phi_x, 
                                                         const adept::adouble weight, 
                                                         std::vector< double > & Jac);

template void compute_jacobian_outside_integration_loop< double >(adept::Stack & stack,
                                                                  const std::vector< double > & solu,
                                                                  const std::vector< double > & Res,
                                                                  std::vector< double > & Jac,
                                                                   const std::vector< int > & loc_to_glob_map,
                                                                   NumericVector*           RES,
                                                                   SparseMatrix*             KK
                                                                   );

template void compute_jacobian_outside_integration_loop < adept::adouble >(adept::Stack & stack,
                                               const std::vector< adept::adouble > & solu,
                                               const std::vector< adept::adouble > & Res,
                                               std::vector< double > & Jac,
                                               const std::vector< int > & loc_to_glob_map,
                                               NumericVector*           RES,
                                               SparseMatrix*             KK
                                                                   );
  

   } //end namespace assemble_jacobian


}  //end namespace

/*=========================================================================

 Program: FEMUS
 Module: FemusInit
 Authors: Eugenio Aulisa, Giorgio Bornia
 
 Copyright (c) FEMTTU
 All rights reserved. 

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_quadrature_GaussPoints_hpp__
#define __femus_quadrature_GaussPoints_hpp__



#include <string>
#include <iostream>
#include <iomanip>




namespace femus {

  
    
  class Gauss {
     
  public:

    Gauss(const char *geom_elem, const char *order_gauss);
    
  inline const double *  GetGaussWeightsPointer() const {
    return GaussWeight;
  };
  
  inline const double *  GetGaussCoordinatePointer (const unsigned &k) const {
    return GaussWeight + (k+1) * GaussPoints;
  };
  
  
  inline const double  GetGaussWeight(const unsigned ig) const {
    return GaussWeight[ig];
  };
  
  inline const unsigned GetGaussPointsNumber() const {
      return GaussPoints;
  };     

  inline const std::string  GetGaussOrderString() const {
    return _order;
  };
  
  inline int  GetGaussOrderIdx() const {
    return gauss_order;
  };
  
  
  
#define FEMUS_QUAD_DIM  3
  static void print_tensor_product_1D ( const double  weights_and_nodes_in[][6], const unsigned n_points_1D){


             std::cout << std::setprecision(17);

             const bool print_weights =false;
             const bool print_nodes_x =false;
             const bool print_nodes_y =false;
             const bool print_nodes_z =true;


       for ( unsigned int p =0; p< n_points_1D; p++){

         for ( unsigned int q =0; q< n_points_1D; q++){

#if FEMUS_QUAD_DIM == 3
          for ( unsigned int r =0; r< n_points_1D; r++){
#endif

//              std::cout << p << " " << q << std::endl ;
           if(print_weights) {  std::cout << weights_and_nodes_in[0][p] * weights_and_nodes_in[0][q]
 #if FEMUS_QUAD_DIM == 3
       * weights_and_nodes_in[0][r]
#endif
               ; }
           if(print_nodes_x) {  std::cout << weights_and_nodes_in[1][p] ; }
           if(print_nodes_y) {  std::cout << weights_and_nodes_in[1][q] ; }
 #if FEMUS_QUAD_DIM == 3
           if(print_nodes_z) {  std::cout << weights_and_nodes_in[1][r] ; }
#endif

//            std::cout << endl ;
           std::cout << ", " ;

#if FEMUS_QUAD_DIM == 3
          }
#endif
        }

    }



}


  protected:
    
    int gauss_order;
    
    std::string _order;
    
    unsigned GaussPoints;
    
    const double * GaussWeight;
    
    
  public:
    
    static constexpr unsigned _NUM_QUAD_RULE_HEX_QUAD_LINE = 6;
   
  };
     

  






  

     
} //end namespace femus     


#endif

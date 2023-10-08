#ifndef __femus_solution_functions_over_domains_or_mesh_files_hpp__ 
#define __femus_solution_functions_over_domains_or_mesh_files_hpp__







template < class type = double >
class Zero : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return 0.;
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  =  0.;
        solGrad[1]  =  0.;

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return  0.; 
    }

      
};




namespace  Domain_square_01by01  {
    


template < class type = double >
class Function_NonZero_on_boundary_1 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return x[0] * x[0] * (1. - x[0]) + sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  =  - x[0] * x[0] +  (1. - x[0]) * 2. * x[0]  + pi * cos( pi * (x[0]) ) * sin( pi * (x[1]) );
        solGrad[1]  =                                              pi * sin( pi * (x[0]) ) * cos( pi * (x[1]) );

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return   - 2. * x[0] +  2. * (1. - x[0])  -  2. * x[0]  - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) ) - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }



  private: 
    
   static constexpr double pi = acos(-1.);
      
};


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return + sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  = pi * cos( pi * (x[0]) ) * sin( pi * (x[1]) );
        solGrad[1]  = pi * sin( pi * (x[0]) ) * cos( pi * (x[1]) );

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return -pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) ) - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }



  private: 
    
   static constexpr double pi = acos(-1.);
      
};


//this solution shows SUPERCONVERGENCE for SERENDIPITY FE, and it is like SUPER PERFECT for BIQUADRATIC FE... it is because of the MESH!
template < class type = double >
class Function_Zero_on_boundary_2 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return  x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  = (1. - 2. * x[0]) *  x[1] * (1. - x[1]);
        solGrad[1]  = (1. - 2. * x[1]) *  x[0] * (1. - x[0]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return -2. * ( x[0] * (1. - x[0])  + x[1] * (1. - x[1]) );
    }



};


//this solution does not have SUPERCONVERGENCE even with the straight mesh
template < class type = double >
class Function_Zero_on_boundary_3 : public Math::Function< type > {

public:

// manufactured Laplacian =============
    type value(const std::vector < type >& x) const {
        
        return    x[0] *  x[0] * (1. - x[0]) * x[1] * (1. - x[1]) ;
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  =  x[0] * (1. - 2. * x[0]) *  x[1] * (1. - x[1]) +  x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
        solGrad[1]  =  x[0] * (1. - 2. * x[1]) *  x[0] * (1. - x[0]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return     x[0] *  (  -2. *  x[1] * (1. - x[1])   ) + 2. *  (1. - 2. * x[0]) *  x[1] * (1. - x[1])   +  x[0] * ( -2. *  x[0] * (1. - x[0])) ;
    }



};








} //end namespace




namespace  Domain_square_01by01_Mesh_Distorted  {


template < class type = double >
class Function_Zero_on_boundary_Continuous0_NotC1_1 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
       type eps = 1.e-5;
        type val;
        type discontinuing_factor = 1.;//fabs(x[1] - 10. * x[0] + 4.);
        
        if (x[1] - 10. * x[0] + 4. > - eps) val = discontinuing_factor * ( x[1] - 10. * x[0] + 4. ) *  (1. - x[1]) * x[1] * x[0];
         else                               val = -1. * discontinuing_factor * ( x[1] - 10. * x[0] + 4. ) *  (1. - x[1]) * x[1] * (1. - x[0]) ;
//         if (x[0] < 0.5 + eps) val = x[0]        *  x[1] * (1. - x[1]);
//          else                 val = ( 1. - x[0] ) *  x[1] * (1. - x[1]);
       
        return val;
    }


    vector < type >  gradient(const std::vector < type >& x) const {
//   std::cout << "Redo calc"; @todo

        vector < type > solGrad(x.size());

        solGrad[0]  = pi * cos( pi * (x[0]) ) * sin( pi * (x[1]) );
        solGrad[1]  = pi * sin( pi * (x[0]) ) * cos( pi * (x[1]) );

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
//   std::cout << "Redo calc"; @todo
        
        return -pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) ) - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }



  private: 
    
   static constexpr double pi = acos(-1.);
      
};



} //end namespace





namespace  Domain_square_01by01_Mesh_Straight  {
  
    
//this solution does not have SUPERCONVERGENCE even with the straight mesh
template < class type = double >
class Function_NonZero_on_boundary_Continuous0_1 : public Math::Function< type > {

public:

// manufactured Laplacian =============
    type value(const std::vector < type >& x) const {
        
        type eps = 1.e-5;
        type val;
        if (x[0] < 0.5 + eps) val = x[0] * (1. - x[0]);
         else                 val = - 0.25 + x[0];
       
        return  val;
    }


    vector < type >  gradient(const std::vector < type >& x) const {

  std::cout << "Redo calc";
  
        vector < type > solGrad(x.size());
        solGrad[0]  =  x[0] * (1. - 2. * x[0]) *  x[1] * (1. - x[1]) +  x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
        solGrad[1]  =  x[0] * (1. - 2. * x[1]) *  x[0] * (1. - x[0]);
   
        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
  std::cout << "Redo calc";
        
        return     x[0] *  (  -2. *  x[1] * (1. - x[1])   ) + 2. *  (1. - 2. * x[0]) *  x[1] * (1. - x[1])   +  x[0] * ( -2. *  x[0] * (1. - x[0])) ;
    }



};

    
template < class type = double >
class Function_Zero_on_boundary_Continuous1_1 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
       type eps = 1.e-5;
        type val;
        if (x[0] < 0.5 + eps) val = (0.5 * x[0]  *  x[0]  + 0.25  )  *  x[1] * (1. - x[1]);
         else                 val = ( x[0] - 0.5 * x[0] * x[0]  ) *  x[1] * (1. - x[1]);
       
        return val;
    }


    vector < type >  gradient(const std::vector < type >& x) const {
//   std::cout << "Redo calc";  @todo

        vector < type > solGrad(x.size());

        solGrad[0]  = pi * cos( pi * (x[0]) ) * sin( pi * (x[1]) );
        solGrad[1]  = pi * sin( pi * (x[0]) ) * cos( pi * (x[1]) );

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
//   std::cout << "Redo calc";  @todo
        
        return -pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) ) - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }



  private: 
    
   static constexpr double pi = acos(-1.);
      
};
    
    
} //end namespace





namespace  Domain_square_m05p05  {

 template < class type = double >
class Function_Zero_on_boundary_4 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
  return cos(pi * x[0]) * cos(pi * x[1]);
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  = -pi * sin(pi * x[0]) * cos(pi * x[1]);
        solGrad[1]  = -pi * cos(pi * x[0]) * sin(pi * x[1]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return -pi * pi * cos(pi * x[0]) * cos(pi * x[1]) - pi * pi * cos(pi * x[0]) * cos(pi * x[1]);
    }



  private: 
    
   static constexpr double pi = acos(-1.);
      
};

   
}    




namespace  Domain_L_shaped  {
    
  template < class type = double >
class Function_NonZero_on_boundary_2 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return x[0] * x[0] * (1. - x[0]) + sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }


    vector < type >  gradient(const std::vector < type >& x) const {

        vector < type > solGrad(x.size());

        solGrad[0]  =  - x[0] * x[0] +  (1. - x[0]) * 2. * x[0]  + pi * cos( pi * (x[0]) ) * sin( pi * (x[1]) );
        solGrad[1]  =                                              pi * sin( pi * (x[0]) ) * cos( pi * (x[1]) );

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return   - 2. * x[0] +  2. * (1. - x[0])  -  2. * x[0]  - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) ) - pi * pi * sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }



  private: 
    
   static constexpr double pi = acos(-1.);
      
};
  
    
}









#endif

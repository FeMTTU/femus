#ifndef __femus_solution_functions_over_domains_or_mesh_files_hpp__ 
#define __femus_solution_functions_over_domains_or_mesh_files_hpp__


#include "Function.hpp"

using namespace femus;

// Functions for every domain - BEGIN ===============================

template < class type = double >
class Zero : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return 0.;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        solGrad[0]  =  0.;
        solGrad[1]  =  0.;

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return  0.; 
    }

      
};

// Functions for every domain - END ===============================


// Functions for 1D domains - BEGIN ===============================



namespace segment_0x1 {



template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
    return  x[0] * (1. - x[0]);
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        solGrad[0]  =  1. - 2. * x[0];

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
    return  -2.;
    }

    
    
};


}



// Functions for 1D domains - END ===============================



// Functions for 2D domains - BEGIN ===============================

namespace  Domain_square_01by01  {
    


template < class type = double >
class Function_NonZero_on_boundary_1 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return x[0] * x[0] * (1. - x[0]) + sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

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
        
    return x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        solGrad[0]  = (1. - 2. * x[0]) *  x[1] * (1. - x[1]);
        solGrad[1]  = (1. - 2. * x[1]) *  x[0] * (1. - x[0]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
    return -2. * ( x[0] * (1. - x[0])  + x[1] * (1. - x[1]) );
    }

    
    
};



//this solution shows SUPERCONVERGENCE for SERENDIPITY FE, and it is like SUPER PERFECT for BIQUADRATIC FE... it is because of the MESH!
template < class type = double >
class Function_Zero_on_boundary_2 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return  x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

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


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        solGrad[0]  =  x[0] * (1. - 2. * x[0]) *  x[1] * (1. - x[1]) +  x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
        solGrad[1]  =  x[0] * (1. - 2. * x[1]) *  x[0] * (1. - x[0]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return     x[0] *  (  -2. *  x[1] * (1. - x[1])   ) + 2. *  (1. - 2. * x[0]) *  x[1] * (1. - x[1])   +  x[0] * ( -2. *  x[0] * (1. - x[0])) ;
    }



};


template < class type = double >
class Function_Zero_on_boundary_4 : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {
        
        return + sin( pi * (x[0]) ) * sin( pi * (x[1]) );
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

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


    std::vector < type >  gradient(const std::vector < type >& x) const {
//   std::cout << "Redo calc"; @todo

        std::vector < type > solGrad(x.size(), 0.);

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


    std::vector < type >  gradient(const std::vector < type >& x) const {

  std::cout << "Redo calc";
  
        std::vector < type > solGrad(x.size(), 0.);
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


    std::vector < type >  gradient(const std::vector < type >& x) const {
//   std::cout << "Redo calc";  @todo

        std::vector < type > solGrad(x.size(), 0.);

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


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

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


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

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


///@todo
namespace circle {


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
        // return    ;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        // solGrad[0]  = ;
        // solGrad[1]  = ;

        // return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        // return ;
    }

    
    
};



}


namespace semicircle_centered_at_0_by_0 {


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
    double xxx = x[0];
    double yyy = x[1];
    const double r2 = xxx * xxx + yyy * yyy;
    double r = (1. - r2) * yyy;
    // double r = (1. - x[0] * x[0] - x[1] * x[1]) * x[1];
    
    return r;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        solGrad[0]  = -2. * x[0] * x[1];
        solGrad[1]  = - 2. * x[1] * x[1] + 1. - x[0] * x[0] - x[1] * x[1];

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
   return  - 8. * x[1]; 
    }

    
    
};


}


namespace quarter_circle_centered_at_0_by_0 {


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
    // for a quarter-circle in Quadrant 1
    
    double xx = x[0];
    double yy = x[1];
    
    return  xx * yy * (1.0 - (xx*xx + yy*yy) ); // forced to be zero on the x and y axis, and the circle edge
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

    double xx = x[0];
    double yy = x[1];
    
        solGrad[0]  = yy * (1.0 - (xx*xx + yy*yy) ) + xx * yy * (-2. * xx);
        solGrad[1]  = xx * (1.0 - (xx*xx + yy*yy) ) + yy * xx * (-2. * yy);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
    double xx = x[0];
    double yy = x[1];
    
    return  -12. * xx * yy;
    }

    
    
};

 
 
}


namespace annulus_centered_at_0_by_0 {


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
     double r2 = x[0] * x[0] + x[1] * x[1];
     double res = (1. - r2) * ( r2 - 0.25 );
     return res;
        
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        double r2 = x[0] * x[0] + x[1] * x[1];
  
        solGrad[0]  = ( - 2. * x[0] ) * ( r2 - 0.25 ) +   (1. - r2) * ( 2. * x[0] );
        solGrad[1]  = ( - 2. * x[1] ) * ( r2 - 0.25 ) +   (1. - r2) * ( 2. * x[1] );

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
  double r2 = x[0] * x[0] + x[1] * x[1];
  double res = -8. * r2 + 4. * (1. - r2) - 4. * (r2 - 0.25);
  return res;
    }

    
    
};


 
 
}



namespace semiannulus_centered_at_0_by_0_cut_along_y {


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
     double r2 = x[0] * x[0] + x[1] * x[1];

     return   x[0] * (1. - r2 ) * ( r2 - 0.25);
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);
        
     double r2 = x[0] * x[0] + x[1] * x[1];

        solGrad[0]  =           (1. - r2 ) * ( r2 - 0.25) + (x[0]) * ( - 2. * x[0] ) * ( r2 - 0.25) + (x[0]) * (1. - r2 ) * ( 2. * x[0] );
        solGrad[1]  =                                     + (x[0]) * ( - 2. * x[1] ) * ( r2 - 0.25) + (x[0]) * (1. - r2 ) * ( 2. * x[1] );

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
         return -2. * x[0] * ( -5. + 12. * x[0] * x[0] + 12. * x[1] * x[1] );
    }

    
    
};


 
 
 
}



// Functions for 2D domains - END ===============================



// Functions for 3D domains - BEGIN ===============================

namespace cube_01_by_01_by_01 {


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
  return x[0] * (1. - x[0]) * x[1] * (1. - x[1]) * x[2] * (1. - x[2]);
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        solGrad[0]  = (1. - 2. * x[0]) * x[1] * (1. - x[1]) * x[2] * (1. - x[2]);
        solGrad[1]  = (1. - 2. * x[1]) * x[2] * (1. - x[2]) * x[0] * (1. - x[0]);
        solGrad[2]  = (1. - 2. * x[2]) * x[0] * (1. - x[0]) * x[1] * (1. - x[1]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        return -2. * ( x[0] * (1. - x[0])  + x[1] * (1. - x[1]) +  x[2] * (1. - x[2]) );
    }

    
    
};




}



namespace cylinder_along_z_with_base_centered_at_1_by_1 {
 

template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
      double xx = x[0];
      double yy = x[1];
      double zz = x[2];
      double r = zz * (2. - zz) * ( 1. - (xx - 1.) * (xx - 1.) - (yy - 1.) * (yy - 1.) ) ;
      return r;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

      double xx = x[0];
      double yy = x[1];
      double zz = x[2];
        
        solGrad[0]  = zz * (2. - zz) * (  - (2. * xx - 2. ) ) ;
        solGrad[1]  = zz * (2. - zz) * (  - (2. * yy - 2. ) ) ;
        solGrad[2]  = (2. - 2. * zz) * ( 1. - (xx - 1.) * (xx - 1.) - (yy - 1.) * (yy - 1.) ) ;

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
   return  4. * x[2] * (x[2] - 2.)  +  2. * ( (x[0] - 1.) * (x[0] - 1.)  +  (x[1] - 1.) * (x[1] - 1.) - 1.);
    }

    
};


}


///@todo
namespace semicylinder {


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
        // return    ;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        // solGrad[0]  = ;
        // solGrad[1]  = ;
        // solGrad[2]  = ;

        // return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        // return ;
    }

    
    
};

 
 
 
}



namespace quarter_cylinder_along_z_with_base_centered_at_0_by_0 {


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
    double xx = x[0];
    double yy = x[1];
    double zz = x[2];
    
     return  xx * yy * zz * (2. - zz) * (-xx * xx - yy * yy + 1.);
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

    double xx = x[0];
    double yy = x[1];
    double zz = x[2];
    
        solGrad[0]  = yy * zz * (2. - zz) * (-xx * xx - yy * yy + 1.) + xx * yy * zz * (2. - zz) * (- 2. * xx );
        solGrad[1]  = xx * zz * (2. - zz) * (-xx * xx - yy * yy + 1.) + xx * yy * zz * (2. - zz) * (- 2. * yy );
        solGrad[2]  = xx * yy * (2. - 2. * zz) * (-xx * xx - yy * yy + 1.);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
    double xx = x[0];
    double yy = x[1];
    double zz = x[2];
    
    return  2. * xx * yy * (-1. + xx * xx + yy * yy - 12. * zz + 6. * zz * zz);
    }

    
    
};


 
 }


namespace prism_annular_base_along_z_with_base_centered_at_0_by_0 {
 
 

template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
     return  x[2] * (1. - x[2] ) * (1. - x[0] * x[0] - x[1] * x[1]) * (-0.25 + x[0] * x[0] + x[1] * x[1]);
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        solGrad[0]  = x[2] * (1. - x[2] ) * ( ( - 2. * x[0] ) * (-0.25 + x[0] * x[0] + x[1] * x[1]) + (1. - x[0] * x[0] - x[1] * x[1]) * ( 2. * x[0] ) );
        solGrad[1]  = x[2] * (1. - x[2] ) * ( ( - 2. * x[1] ) * (-0.25 + x[0] * x[0] + x[1] * x[1]) + (1. - x[0] * x[0] - x[1] * x[1]) * ( 2. * x[1] ) );
        solGrad[2]  =   (1. - 2. * x[2] ) * (1. - x[0] * x[0] - x[1] * x[1]) * (-0.25 + x[0] * x[0] + x[1] * x[1]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
    return  2. * ( 0.25 + x[0] * x[0] * x[0] * x[0] + x[1] * x[1] * x[1] * x[1]  + 2.5 * x[2] - 2.5 * x[2] * x[2] +  
      x[1] * x[1] * (-1.25 - 8. * x[2] + 8. * x[2] * x[2])
    + x[0] * x[0] * (-1.25 + 2. * x[1] * x[1]  - 8. * x[2] + 8. * x[2] * x[2]) );
    // return  2. * (0.25 + x^4 + y^4 + 2.5 z - 2.5 z^2 + y^2 (-1.25 - 8 z + 8 z^2) + x^2 (-1.25 + 2 y^2 - 8 z + 8 z^2));
    }

    
    
};


 
 
 
}

// Functions for 3D domains - END ===============================



#endif

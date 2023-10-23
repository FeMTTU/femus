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

// Functions for every domain - END ===============================


// 1D - BEGIN ===============================



namespace segment_0x1 {

  
  namespace function_0 {
    

double value(const std::vector<double> & x) {
    
    return  x[0] * (1. - x[0]);
}


// user-made equation - accepts only coordinates
double laplacian(const std::vector<double> & x){
    
    return  -2.;
     }

  }


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
    return  x[0] * (1. - x[0]);
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size());

        solGrad[0]  =  1. - 2. * x[0];

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
    return  -2.;
    }

    
    
};


}



// 1D - END ===============================



// 2D - BEGIN ===============================
namespace square_01_by_01 {
  
  
  namespace function_0 {

double value(const std::vector < double >& x) {
    
  return x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
    
}

double laplacian(const std::vector < double >& x) {
    
  return -2. * ( x[0] * (1. - x[0])  + x[1] * (1. - x[1]) );
  
}


}





template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
    return x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size());

        solGrad[0]  = (1. - 2. * x[0]) *  x[1] * (1. - x[1]);
        solGrad[1]  = (1. - 2. * x[1]) *  x[0] * (1. - x[0]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
    return -2. * ( x[0] * (1. - x[0])  + x[1] * (1. - x[1]) );
    }

    
    
};



}



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


    std::vector < type >  gradient(const std::vector < type >& x) const {

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


namespace circle {


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
        // return    ;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size());

        // solGrad[0]  = ;
        // solGrad[1]  = ;

        // return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        // return ;
    }

    
    
};



}


namespace semicircle {

namespace function_0 {


double value(const std::vector < double >& x) {
    
    double xxx = x[0];
    double yyy = x[1];
    const double r2 = xxx * xxx + yyy * yyy;
    double r = (1. - r2) * yyy;
    
    return r;
   
}

double laplacian(const std::vector < double >& x) {

   return  - 8. * x[1]; 
    
}


}




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

        std::vector < type > solGrad(x.size());

        solGrad[0]  = -2. * x[0] * x[1];
        solGrad[1]  = - 2. * x[1] * x[1] + 1. - x[0] * x[0] - x[1] * x[1];

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
   return  - 8. * x[1]; 
    }

    
    
};


}


namespace quarter_circle {

 namespace function_0 {
   
double value(const std::vector<double> & x_qp){
    
    // for a quarter-circle in Quadrant 1
    
    double x = x_qp[0];
    double y = x_qp[1];
    
    return  x * y * (1.0 - (x*x + y*y)); // forced to be zero on the x and y axis, and the circle edge
}

 

double laplacian(const std::vector<double> & x_qp){
    
    // for a quarter-circle in Quadrant 1
    
    double x = x_qp[0];
    double y = x_qp[1];
    
    return  -12. * x * y;
}


 } 




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

        std::vector < type > solGrad(x.size());

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


namespace annulus {
  
 namespace function_0 {

double value(const std::vector < double >& x) {
    
  double r2 = x[0] * x[0] + x[1] * x[1];
  double res = (1. - r2) * ( r2 - 0.25 );
//   double yprime = (1. - r2)' * ( r2 - 0.25 ) +   (1. - r2) * ( r2 - 0.25 )'; 
  return res;

    
}


double laplacian(const std::vector < double >& x) {
    
  double r2 = x[0] * x[0] + x[1] * x[1];
  double res = -8. * r2 + 4. * (1. - r2) - 4. * (r2 - 0.25);
//   double y = 16. * (0.3125 - r2);
  return res;
  
}


 }

 


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
        // return    ;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size());

        // solGrad[0]  = ;
        // solGrad[1]  = ;

        // return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        // return ;
    }

    
    
};


 
 
}



namespace semiannulus {
  
 namespace function_0 {
   

double value(const std::vector < double >& x) {
    
     double r2 = x[0] * x[0] + x[1] * x[1];

     return   x[0] * (1. - sqrt(r2) ) * ( sqrt(r2) - .5);
//      return   x[0] * (r2 - 1. ) * (.25 - r2);
}


double laplacian(const std::vector < double >& x) {
  
    double r2 = x[0] * x[0] + x[1] * x[1];
    double temp = -x[0] * (8. - 4.5 / (sqrt(r2)));
    //double temp = (4. - 1.5 / (sqrt(r2)));
  return temp;

//   return - 20. * x[0] * (-0.45 + x[0] * x[0] - 0.6 * x[1] * x[1] );
    
}


 }
 

 

template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
        // return    ;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size());

        // solGrad[0]  = ;
        // solGrad[1]  = ;

        // return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        // return ;
    }

    
    
};


 
 
 
}



// 2D - END ===============================



// 3D - BEGIN ===============================

namespace cube_01_by_01_by_01 {

namespace function_0 {


double value(const std::vector < double >& x) {
    
  return x[0] * (1. - x[0]) * x[1] * (1. - x[1]) * x[2] * (1. - x[2]);
    
}


double laplacian(const std::vector < double >& x) {
    
  return -2. * ( x[0] * (1. - x[0])  + x[1] * (1. - x[1]) +  x[2] * (1. - x[2]) );
  
}


}





template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
        // return    ;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size());

        // solGrad[0]  = ;
        // solGrad[1]  = ;

        // return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        // return ;
    }

    
    
};




}



namespace cylinder {

 namespace function_0 {
   
   
double value(const std::vector<double> & xxx){
  
      double x = xxx[0];
      double y = xxx[1];
      double z = xxx[2];
     double r = z * (2. - z) * ( 1 - (x-1) * (x-1) - (y-1) * (y-1) ) ;
  return r;

     
}
 
double  laplacian(const std::vector < double >& x_qp) {
    
   double r = 4*x_qp[2]*(x_qp[2] - 2.)  +  2*( (x_qp[0] - 1.)*(x_qp[0] - 1.)  +  (x_qp[1] - 1.)*(x_qp[1] - 1.) - 1.);
  return r;
}

 }
 
 

template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
        // return    ;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size());

        // solGrad[0]  = ;
        // solGrad[1]  = ;

        // return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        // return ;
    }

    
    
};


}


namespace semicylinder {

 namespace function_0 {
  

double value(const std::vector < double >& x) {
}


double laplacian(const std::vector < double >& x) {
}



 }

 


template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
        // return    ;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size());

        // solGrad[0]  = ;
        // solGrad[1]  = ;

        // return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        // return ;
    }

    
    
};

 
 
 
}



namespace quarter_cylinder {

 namespace function_0 {
  


double value(const std::vector < double >& x_qp) {

    double x = x_qp[0];
    double y = x_qp[1];
    double z = x_qp[2];
    
     return  x*y*z * (2.0 - z)*(-x*x - y*y + 1.0);
}


double laplacian(const std::vector < double >& x_qp) {

      
    // Quarter cylinder of radius 1 and length 2
    
    double x = x_qp[0];
    double y = x_qp[1];
    double z = x_qp[2];
    
    // Function = x*y*z*(z-2.0)*(x*x + y*y - 1.0)
    
    // Return -Delta U0
    return ( - 12.0 * x * y * z * (2.0 - z) + 2 * x * y * ( x * x + y * y ) );
  
    
    
}

 }
 
 
 

template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
        // return    ;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size());

        // solGrad[0]  = ;
        // solGrad[1]  = ;

        // return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        // return ;
    }

    
    
};


 
 }


namespace prism_annular_base {

 namespace function_0 {



 double value(const std::vector<double> & x) {
     
     return  x[2] * (1. - x[2] ) * (1. - x[0]*x[0] - x[1]*x[1]) * (-0.25 + x[0]*x[0] + x[1]*x[1]);
 }
 

//calculator: f = z(z - 1)(1 - x^2 - y^2)(1/4 - x^2 - y^2);
double laplacian(const std::vector < double > & x) {
    
  double r2 = 0.5 - 2.5*pow(x[0],2) + 2*pow(x[0],4) - 2.5*pow(x[1],2) + 4*pow(x[0],2)*pow(x[1],2) + 2*pow(x[1],4) + 5.*x[2] - 16*pow(x[0],2)*x[2] - 16*pow(x[1],2)*x[2] - 5.*pow(x[2],2) + 16*pow(x[0],2)*pow(x[2],2) + 16*pow(x[1],2)*pow(x[2],2);
  
   return  r2;
   
  }
 
 }
 
 
 

template < class type = double >
class Function_Zero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {
        
        // return    ;
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size());

        // solGrad[0]  = ;
        // solGrad[1]  = ;

        // return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {
        
        // return ;
    }

    
    
};


 
 
 
}

// 3D - END ===============================



#endif

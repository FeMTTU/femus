#ifndef __femus_solution_functions_over_domains_or_mesh_files1_hpp__ 
#define __femus_solution_functions_over_domains_or_mesh_files1_hpp__


#include "Function.hpp"

using namespace femus;


namespace Domains {


// Functions for 2D domains - BEGIN ===============================


namespace  square_m05p05  {

    
template < class type = double >
class Function_Zero_on_boundary_4_Laplacian : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {

        return  -2.* pi * pi * cos(pi * x[0]) * cos(pi * x[1]);
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        solGrad[0]  = 2. * pi * pi * pi * sin(pi * x[0]) * cos(pi * x[1]);
        solGrad[1]  = 2. * pi * pi * pi * cos(pi * x[0]) * sin(pi * x[1]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {

        return  4. * pi * pi * pi * pi * cos(pi * x[0]) * cos(pi * x[1]);
    }



  private:

   static constexpr double pi = acos(-1.);

};



template < class type = double >
class Function_NonZero_on_boundary_4 : public Math::Function< type > {

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
class Function_NonZero_on_boundary_4_Laplacian : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {

        return  2. - 6. * x[0] - 2. *  pi * pi * sin(pi * x[0]) * sin(pi * x[1]);
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        solGrad[0]  = -6. - 2. * pi * pi * pi * cos(pi * x[0]) * sin(pi * x[1]);
        solGrad[1]  = -2. * pi * pi * pi * sin(pi * x[0]) * cos(pi * x[1]);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {

        return  4. * pi * pi * pi * pi * sin(pi * x[0]) * sin(pi * x[1]);
    }



  private:

   static constexpr double pi = acos(-1.);

};
    
    

   
}    



} //end Domains

#endif

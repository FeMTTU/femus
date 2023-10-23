#ifndef __femus_utils_Math_Function_hpp__
#define __femus_utils_Math_Function_hpp__


#include <vector>


namespace femus {


namespace Math {
    
    

template < class type = double >
  class Function {
 
  public:
      
 virtual type value(const std::vector < type >& x) const {}/*= 0*/;

 virtual std::vector < type >  gradient(const std::vector < type >& x) const  {}/*= 0*/;

 virtual type laplacian(const std::vector < type >& x) const {}/*= 0*/;
 
  type helmholtz(const std::vector < type >& x) const { return ( - laplacian(x) + value(x) ); };

};




}

}



#endif

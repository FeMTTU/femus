#ifndef __femus_app_specifics_hpp__
#define __femus_app_specifics_hpp__


#include <vector>
#include <string>


namespace femus {
    
    
class MultiLevelProblem;


class app_specifics {
 
  public:
      
      app_specifics() {  _mesh_files.resize(2); }
  
   std::vector< std::string >   _mesh_files;   
 
   //func pointer of Boundary Conditions: typedef
    typedef bool (*BoundaryFunction) (const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double &value, const int FaceName, const double time);
    
   //func pointer of Boundary Conditions
    BoundaryFunction   _bdry_func;

   //typedef for func pointer of Boundary Conditions

    
};


}



#endif

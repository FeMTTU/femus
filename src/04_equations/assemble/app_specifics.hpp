#ifndef __femus_app_specifics_hpp__
#define __femus_app_specifics_hpp__


#include <vector>
#include <string>

#include "CurrentElem.hpp"
#include "ElemType_template.hpp"

namespace femus {
    
    
class MultiLevelProblem;
class Mesh;
class MultiLevelSolution;


   /// @todo this should be templated as a class, or some of its functions should be templated
class app_specifics {
 
  public:
      
      app_specifics() {  }
  
   std::vector< std::string >   _mesh_files;  //same domain, only potentially multiple mesh discretizations 
 
   std::string  _system_name;  //for now we only accept 1 System in a certain App. This name is needed to retrieve the equation from the Problem
   
   //func pointer of EQUATION
   /// @todo this should be templated
   typedef void (* AssembleFunctionType) (MultiLevelProblem &ml_prob);

   AssembleFunctionType  _assemble_function;


   //func pointer of RHS
   typedef  double    (* AssembleFunctionRHS )  (const std::vector<double> & x_qp);
   
   AssembleFunctionRHS  _assemble_function_rhs;

    //func pointer of true solution
   typedef  double    (* TrueSolution )  (const std::vector<double> & x_qp);
   
   TrueSolution  _true_solution;
  
   //func pointer of Boundary Conditions
    typedef bool (*BoundaryFunction) (const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double &value, const int FaceName, const double time);
    
    BoundaryFunction   _boundary_conditions_types_and_values;

    
    
   
};


}



#endif

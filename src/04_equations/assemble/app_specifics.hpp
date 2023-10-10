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
      
      app_specifics() {
         
         _assemble_function = NULL;
         _assemble_function_rhs = NULL;
         _boundary_conditions_types_and_values = NULL;
         
         _true_solution = NULL;
         
      }
  
   std::vector< std::string >   _mesh_files;  //same domain, only potentially multiple mesh discretizations 
 
   std::string  _system_name;  //for now we only accept 1 System in a certain App. This name is needed to retrieve the equation from the Problem
   
   //func pointer of EQUATION - BEGIN
   /// @todo this should be templated
   typedef void (* AssembleFunctionType) (MultiLevelProblem &ml_prob);

   AssembleFunctionType  _assemble_function;
   //func pointer of EQUATION - END


   //func pointer of RHS - BEGIN
   typedef  double    (* AssembleFunctionRHS )  (const std::vector<double> & x_qp);
   
   AssembleFunctionRHS  _assemble_function_rhs;
   //func pointer of RHS - END

   //func pointer of Boundary Conditions - BEGIN
    typedef bool (*BoundaryFunction) (const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double &value, const int FaceName, const double time);
    
    BoundaryFunction   _boundary_conditions_types_and_values;
   //func pointer of Boundary Conditions - END


    //func pointer of true solution - BEGIN
   typedef  double    (* TrueSolution )  (const std::vector<double> & x_qp);
   
   TrueSolution  _true_solution;
    //func pointer of true solution - END
  
    
    
   
};


}



#endif

#ifndef __femus_app_specifics_hpp__
#define __femus_app_specifics_hpp__


#include "CurrentElem.hpp"
#include "ElemType_template.hpp"
#include "Function.hpp"

#include <vector>
#include <string>



namespace femus {
    
    
class MultiLevelProblem;
class Mesh;
class MultiLevelSolution;


   /// @todo this should be templated as a class, or some of its functions should be templated
class app_specifics {
 
  public:
      
      app_specifics() {
         
         _assemble_function = NULL;
         _boundary_conditions_types_and_values = NULL;
         
         _assemble_function_rhs = NULL;
         
         _assemble_function_for_rhs = NULL;
         _true_solution_function = NULL;         
         
      }
  
   //mesh files - BEGIN
   std::vector< std::string >   _mesh_files;  //same domain, only potentially multiple mesh discretizations 

   std::vector< std::string >   _mesh_files_path_relative_to_executable;  //same domain, only potentially multiple mesh discretizations 
   //mesh files - END
 
   //System - BEGIN
   std::string  _system_name;  //for now we only accept 1 System in a certain App. This name is needed to retrieve the equation from the Problem
   //System - END
   
   //func pointer of EQUATION - BEGIN
   /// @todo this should be templated
   typedef void (* AssembleFunctionType) (MultiLevelProblem &ml_prob);

   AssembleFunctionType  _assemble_function;
   //func pointer of EQUATION - END


   //func pointer of Boundary Conditions - BEGIN
    typedef bool (*BoundaryFunction) (const MultiLevelProblem * ml_prob, 
                                      const std::vector < double >& x,
                                      const char name[], 
                                      double &value, 
                                      const int FaceName,
                                      const double time);
    
    
    BoundaryFunction   _boundary_conditions_types_and_values;
   //func pointer of Boundary Conditions - END
    

   //func pointer of RHS - BEGIN
   // [[deprecated]]
   typedef  double    (* AssembleFunctionRHS )  (const std::vector<double> & x_qp);
   
   AssembleFunctionRHS  _assemble_function_rhs;
   //func pointer of RHS - END
   
  
   Math::Function< double > *  _assemble_function_for_rhs;
   
   Math::Function< double > *  _true_solution_function;
    
};


}



#endif

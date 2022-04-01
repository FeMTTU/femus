#ifndef __femus_app_specifics_hpp__
#define __femus_app_specifics_hpp__


#include <vector>
#include <string>

#include "CurrentElem.hpp"

namespace femus {
    
    
class MultiLevelProblem;
class Mesh;
class MultiLevelSolution;


class app_specifics {
 
  public:
      
      app_specifics() {  _mesh_files.resize(2); }
  
   std::vector< std::string >   _mesh_files;  //same domain, only potentially multiple mesh discretizations 
 
 
   //typedef for func pointer of EQUATION
   typedef void (* AssembleFunctionType) (MultiLevelProblem &ml_prob);

   AssembleFunctionType  _assemble_function;


   //typedef for func pointer of RHS
   typedef  double    (* AssembleFunctionRHS )  (const std::vector<double> & x_qp);
   
   AssembleFunctionRHS  _assemble_function_rhs;
   
    //typedef for natural boundary integral loop - 1d
   typedef  void   (* AssembleFunctionNaturalBoundaryLoop1d )   
                     (const MultiLevelProblem *    ml_prob, 
                     const Mesh *                    msh,
                     const MultiLevelSolution *    ml_sol, 
                     const unsigned iel,
                     CurrentElem < double > & geom_element,
                     const unsigned xType,
                     const std::string solname_u,
                     const unsigned solFEType_u,
                     std::vector< double > & Res
                    );
   
   
   AssembleFunctionNaturalBoundaryLoop1d   _assemble_function_natural_boundary_loop_1d;
   
   
   //func pointer of Boundary Conditions: typedef
    typedef bool (*BoundaryFunction) (const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double &value, const int FaceName, const double time);
    
   //func pointer of Boundary Conditions
    BoundaryFunction   _bdry_func;

    
    
   
};


}



#endif

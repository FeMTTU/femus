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
 
 
   //typedef for func pointer of EQUATION
   /// @todo this should be templated
   typedef void (* AssembleFunctionType) (MultiLevelProblem &ml_prob);

   AssembleFunctionType  _assemble_function;


   //typedef for func pointer of RHS
   typedef  double    (* AssembleFunctionRHS )  (const std::vector<double> & x_qp);
   
   AssembleFunctionRHS  _assemble_function_rhs;

    //typedef for func pointer of RHS
   typedef  double    (* NormTrueSolution )  (const std::vector<double> & x_qp);
   
   NormTrueSolution  _norm_true_solution;
  
    //typedef for natural boundary integral loop - 1d
   /// @todo this should be templated
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
   
   
    //typedef for natural boundary integral loop - 2d3d
   /// @todo this should be templated
   typedef  void   (* AssembleFunctionNaturalBoundaryLoop2d3d )   
                      (const MultiLevelProblem *    ml_prob, 
                       const Mesh *                    msh,
                       const MultiLevelSolution *    ml_sol, 
                       const unsigned iel,
                       CurrentElem < double > & geom_element,
                       const unsigned solType_coords,
                       const std::string solname_u,
                       const unsigned solFEType_u,
                       std::vector< double > & Res,
                       //-----------
                       std::vector < std::vector < /*const*/ elem_type_templ_base<double, double/*real_num, real_num_mov*/> *  > >  elem_all,
                       const unsigned dim,
                       const unsigned space_dim,
                       const unsigned max_size
                    );
   
   AssembleFunctionNaturalBoundaryLoop2d3d   _assemble_function_natural_boundary_loop_2d3d;
   
   //func pointer of Boundary Conditions: typedef
    typedef bool (*BoundaryFunction) (const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double &value, const int FaceName, const double time);
    
   //func pointer of Boundary Conditions
    BoundaryFunction   _bdry_func;

    
    
   
};


}



#endif

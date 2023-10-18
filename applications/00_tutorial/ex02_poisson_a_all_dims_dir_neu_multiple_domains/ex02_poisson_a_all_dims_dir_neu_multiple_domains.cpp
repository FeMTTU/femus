/**
 * Solve 
 *     - \Delta u = 1
*/


#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NonLinearImplicitSystem.hpp"

 


#include "00_poisson_eqn_with_dirichlet_or_neumann_bc.hpp"
#include "00_norm_of_errors_of_unknowns.hpp"

#include "app_specifics.hpp"

using namespace femus;

// ======= Domain-related and Mesh-related stuff - BEGIN =========================

// SEGMENT - BEGIN

namespace segment {

  
  namespace function_0 {
    

double value(const std::vector<double> & x) {
    
    // for a 1d segment
    
    return  x[0] * (1. - x[0]);
}


// user-made equation - accepts only coordinates
double laplacian(const std::vector<double> & x){
    
    // for a 1d segment
    
    return  -2.;
     }

  }

  
  
 
// This depends on: 
// mesh file (for the face flags)
// underlying exact solution (if provided)
// equation (name of the unknowns; also, if it is Laplace, biharmonic, Stokes, etc)
// 
// 
bool bc_all_dirichlet(const MultiLevelProblem * ml_prob,
                      const std::vector < double >& x,
                      const char name[], 
                      double& value,
                      const int face_name,
                      const double time) {

  bool dirichlet = false;
  value = 0.;
  
  
 if (ml_prob->GetMLMesh()->GetDimension() != 1 )  abort();
  
  if (face_name == 1) {
      dirichlet = true;
        value = 0.;
    }
  else if (face_name == 2) {
      dirichlet = true;
        value = 0.;
    }

 
  return dirichlet;
  
 }
 
  }
  
  

// SEGMENT - END



// SQUARE - BEGIN
namespace square {
  
  
namespace function_0 {

double value(const std::vector < double >& x) {
    
  return x[0] * (1. - x[0]) * x[1] * (1. - x[1]);
    
}

double laplacian(const std::vector < double >& x) {
    
  return -2. * ( x[0] * (1. - x[0])  + x[1] * (1. - x[1]) );
  
}


}



bool bc_all_dirichlet(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

  if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();
  
   bool dirichlet = false;
  value = 0.;
    
     
  if (face_name == 1) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 2) {
      dirichlet = true;
        value = 0.;
  }

 else  if (face_name == 3) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 4) {
      dirichlet = true;
        value = 0.;
  }
  
  

   return dirichlet;
   
}



}
// SQUARE - END



// CUBE - BEGIN
namespace cube {

namespace function_0 {


double value(const std::vector < double >& x) {
    
  return x[0] * (1. - x[0]) * x[1] * (1. - x[1]) * x[2] * (1. - x[2]);
    
}


double laplacian(const std::vector < double >& x) {
    
  return -2. * ( x[0] * (1. - x[0])  + x[1] * (1. - x[1]) +  x[2] * (1. - x[2]) );
  
}


}



bool bc_all_dirichlet(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

  if (ml_prob->GetMLMesh()->GetDimension() != 3 )  abort();
  
   bool dirichlet = false;
  value = 0.;
    
     
  if (face_name == 1) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 2) {
      dirichlet = true;
        value = 0.;
  }

 else  if (face_name == 3) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 4) {
      dirichlet = true;
        value = 0.;
  }
 else  if (face_name == 5) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 6) {
      dirichlet = true;
        value = 0.;
  }
  
  

   return dirichlet;
   
}





}
// CUBE - END





// CIRCLE - BEGIN
namespace circle {
  
  
namespace function_0 {
  
  
double value(const std::vector < double >& x) {
}

double laplacian(const std::vector < double >& x) {
}


}


bool bc_all_dirichlet(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

 if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();
  


//  return dirichlet;
    
}






}
// CIRCLE - END



// SEMICIRCLE - BEGIN
namespace semicircle {

namespace function_0 {

double laplacian(const std::vector < double >& x) {

   return  - 8. * x[1]; 
    
}

double value(const std::vector < double >& x) {
    
    double xxx = x[0];
    double yyy = x[1];
    const double r2 = xxx * xxx + yyy * yyy;
    double r = (1. - r2) * yyy;
    
    return r;
   
}



}



bool bc_all_dirichlet(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {
    
   if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();

  bool dirichlet = false;
  value = 0.;
   
   
     if (face_name == 1) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 2) {
      dirichlet = true;
        value = 0.;
  }   
  
       return dirichlet;

    
}


}
// SEMICIRCLE - END



// QUARTER CIRCLE - BEGIN
namespace quarter_circle {

 namespace function_0 {
   
double value(const std::vector<double> & x_qp){
    
    // for a quarter-circle in Quadrant 1
    
    double x = x_qp[0];
    double y = x_qp[1];
    
    return  x * y * (1.0 - (x*x + y*y)); // forced to be zero on the x and y axis, and the circle edge
}

 

// this is specifically the laplacian of the function given above
// flynn, user-made equation - accepts only coordinates
double laplacian(const std::vector<double> & x_qp){
    
    // for a quarter-circle in Quadrant 1
    
    double x = x_qp[0];
    double y = x_qp[1];
    
    return  -12. * x * y;
}


 } 
 
 

bool bc_all_dirichlet(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {
    
 if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();
 
  bool dirichlet = false;
  value = 0.;
  
  const double tolerance = 1.e-5;
  
  
    if (face_name == 1) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 2) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 3) {
      dirichlet = true;
        value = 0.;
  }

  
  

  return dirichlet;
  
 }
 
 

 }
// QUARTER CIRCLE - END




// ANNULUS - BEGIN
namespace annulus {
  
 namespace function_0 {

 
double laplacian(const std::vector < double >& x) {
    
  double r2 = x[0] * x[0] + x[1] * x[1];
  double y = -8. * r2 + 4. * (1. - r2) - 4. * (r2 - 0.25);
//   double y = 16. * (0.3125 - r2);
  return y;
  
}

double value(const std::vector < double >& x) {
    
  double r2 = x[0] * x[0] + x[1] * x[1];
  double y = (1. - r2) * ( r2 - 0.25 );
//   double yprime = (1. - r2)' * ( r2 - 0.25 ) +   (1. - r2) * ( r2 - 0.25 )'; 
  return y;

    
}

 }
 
 
 
bool bc_all_dirichlet(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

     if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();
     
  bool dirichlet = false;
  value = 0.;
    

    if (face_name == 1) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 2) {
      dirichlet = true;
        value = 0.;
  }
   
   return dirichlet;
    
}

 
}
// ANNULUS - END



// SEMIANNULUS - BEGIN
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
 
 
 bool bc_all_dirichlet(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

 if (ml_prob->GetMLMesh()->GetDimension() != 2 )  abort();
 
  bool dirichlet = false;
  value = 0.;
  
     
    if (face_name == 1) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 2) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 3) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 4) {
      dirichlet = true;
        value = 0.;
  }
   
 
 
  return dirichlet;
  
 }

 
 

}
// SEMIANNULUS - END





// CYLINDER - BEGIN
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
 
 
bool bc_all_dirichlet(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

 if (ml_prob->GetMLMesh()->GetDimension() != 3 )  abort();

 bool dirichlet = false;
  value = 0.;
  
  
    if (face_name == 1) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 2) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 3) {
      dirichlet = true;
        value = 0.;
  }

 
  return dirichlet;
  
 }
 

 
}
// CYLINDER - END


// SEMICYLINDER - BEGIN
namespace semicylinder {

 namespace function_0 {
  

double laplacian(const std::vector < double >& x) {
}

double value(const std::vector < double >& x) {
}


 }
 
 
 bool bc_all_dirichlet(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

    if (ml_prob->GetMLMesh()->GetDimension() != 3 )  abort();

    
    
//       return dirichlet; 
    
}

 
 
}
// SEMICYLINDER - END


// QUARTER CYLINDER - BEGIN
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
 
 
   
bool bc_all_dirichlet(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

     if (ml_prob->GetMLMesh()->GetDimension() != 3 )  abort();

  bool dirichlet = false;
  value = 0.;
    

     
    if (face_name == 1) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 2) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 3) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 4) {
      dirichlet = true;
        value = 0.;
  }
  else if (face_name == 5) {
      dirichlet = true;
        value = 0.;
  }
 
 
  return dirichlet;
  
}

 
 
}
// QUARTER CYLINDER - END





// PRISM WITH ANNULAR BASE - BEGIN
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
 
 bool bc_all_dirichlet(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

 if (ml_prob->GetMLMesh()->GetDimension() != 3 )  abort();   
 
 bool dirichlet = false;
    value = 0.;


        if (face_name == 1) {
            dirichlet = true;
            value = 0.;
        }
        else if (face_name == 2) {
            dirichlet = true;
            value = 0.;
        }
        else if (face_name == 3) {
            dirichlet = true;
            value = 0.;
        }
        else if (face_name == 4) {
            dirichlet = true;
            value = 0.;
        }

    
    return dirichlet;

}



}
// PRISM WITH ANNULAR BASE - END



// ======= Domain-related and Mesh-related stuff - END =========================




 //Unknown initial condition - BEGIN  ==================
double InitialValueU(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {
    
  return 0.;
  
}
 //Unknown initial condition - END  ==================




 //Unknown definition - BEGIN  ==================
 const std::vector< Unknown >  provide_list_of_unknowns() {
     
     
  std::vector< FEFamily > feFamily;
  std::vector< FEOrder >   feOrder;

                        feFamily.push_back(LAGRANGE);
 
                        feOrder.push_back(/*FIRST*/SECOND);
 

  assert( feFamily.size() == feOrder.size() );
 
 std::vector< Unknown >  unknowns(feFamily.size());

   unknowns[0]._name      = "u";

   unknowns[0]._is_sparse = true;
   
     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
              unknowns[u]._fe_family  = feFamily[u];
              unknowns[u]._fe_order   = feOrder[u];
              unknowns[u]._time_order = 0;
              unknowns[u]._is_pde_unknown = true;
              
     }
 
 
   return unknowns;
     
}
 //Unknown definition - END  ==================





int main(int argc, char** args) {

  // ======= Init ========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Problem ========================
  MultiLevelProblem ml_prob;
  
  // ======= Files - BEGIN  ========================
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = true;
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);

  // ======= Problem, Files ========================
  ml_prob.SetFilesHandler(&files);
  // ======= Files - END  ========================
  
  // ======= Problem, Quad Rule - BEGIN ========================
  std::string fe_quad_rule("seventh");
  
  ml_prob.SetQuadratureRuleAllGeomElems(fe_quad_rule);
  ml_prob.set_all_abstract_fe_multiple();
  // ======= Problem, Quad Rule - END  ========================
  
    // ======= App Specifics - BEGIN  ==================
  std::string system_common_name = "Laplace";
  std::vector< app_specifics >   my_specifics;
  
  app_specifics  app_segment;   //me
  app_specifics  app_square;   //me
//   app_specifics  app_circle;   //Gayani
  app_specifics  app_semicircle;   //Himali
  app_specifics  app_quarter_circle;      //Max
  app_specifics  app_annulus;       //Jon
  app_specifics  app_semiannulus;        //Fahad
  app_specifics  app_cube;   //me
  app_specifics  app_cylinder;            //Aman
//   app_specifics  app_semicylinder;   //Rifat
  app_specifics  app_quarter_cylinder; //Armando   //coarser mesh to be done
  app_specifics  app_prism_annular_base;  //Abu
 
  const std::string relative_path_to_build_directory =  "../../../";
  
  //segment - BEGIN
  app_segment._system_name = system_common_name;
  app_segment._assemble_function = poisson_equation::equation_with_dirichlet_or_neumann_bc<double, double>;
  
  app_segment._mesh_files.push_back("segment_16_dir_neu.med");
  app_segment._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/01_1d/segment/0-1/");
  
  app_segment._boundary_conditions_types_and_values             = segment::bc_all_dirichlet;

  app_segment._true_solution    = segment::function_0::value;
  app_segment._assemble_function_rhs = segment::function_0::laplacian;
  //segment - END
  
  
 //assignment_square - BEGIN
  app_square._system_name = system_common_name;
  app_square._assemble_function                            = poisson_equation::equation_with_dirichlet_or_neumann_bc<double, double>;
  
  app_square._mesh_files.push_back("assignment_square_quadrilateral.med");
  app_square._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/02_2d/square/0-1x0-1/");
  
  app_square._boundary_conditions_types_and_values             = square::bc_all_dirichlet;
  app_square._assemble_function_rhs = square::function_0::laplacian;
  app_square._true_solution    = square::function_0::value;
 //assignment_square - END


  //assignment_semicircle - BEGIN
  app_semicircle._system_name = system_common_name;
  app_semicircle._assemble_function = poisson_equation::equation_with_dirichlet_or_neumann_bc<double, double>;
  
//   app_semicircle._mesh_files.push_back("assignment_semicircle_triangular.med");
  app_semicircle._mesh_files.push_back("assignment_semicircle_quadrilateral.med");
  app_semicircle._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/02_2d/circle_semi/");
  
  app_semicircle._boundary_conditions_types_and_values             = semicircle::bc_all_dirichlet;
  app_semicircle._assemble_function_rhs = semicircle::function_0::laplacian;
  app_semicircle._true_solution    = semicircle::function_0::value;
  //assignment_semicircle - END
   
  
  //assignment_quarter_circle - BEGIN
  app_quarter_circle._system_name = system_common_name;
  app_quarter_circle._assemble_function = poisson_equation::equation_with_dirichlet_or_neumann_bc<double, double>;
  
  app_quarter_circle._mesh_files.push_back("assignment_quarter_circle_triangular.med");
  app_quarter_circle._mesh_files.push_back("assignment_quarter_circle_quadrilateral.med");
  app_quarter_circle._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/02_2d/circle_quarter/");
  app_quarter_circle._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/02_2d/circle_quarter/");
  
  app_quarter_circle._boundary_conditions_types_and_values             = quarter_circle::bc_all_dirichlet;
  app_quarter_circle._assemble_function_rhs = quarter_circle::function_0::laplacian;
  app_quarter_circle._true_solution    = quarter_circle::function_0::value;
  //assignment_quarter_circle - END
 

  //assignment_annulus - BEGIN
  app_annulus._system_name = system_common_name;
  app_annulus._assemble_function = poisson_equation::equation_with_dirichlet_or_neumann_bc<double, double>;
  
//   app_annulus._mesh_files.push_back("assignment_annulus_triangular.med");
  app_annulus._mesh_files.push_back("assignment_annulus_quadrilateral.med");
  app_annulus._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/02_2d/annulus/");
  
  app_annulus._boundary_conditions_types_and_values             = annulus::bc_all_dirichlet;

  app_annulus._assemble_function_rhs = annulus::function_0::laplacian;
  app_annulus._true_solution    = annulus::function_0::value;
  //assignment_annulus - END

  
  //assignment_semiannulus - BEGIN
  app_semiannulus._system_name = system_common_name;
  app_semiannulus._assemble_function = poisson_equation::equation_with_dirichlet_or_neumann_bc<double, double>;
  
  app_semiannulus._mesh_files.push_back("assignment_semiannulus_triangular.med");
  app_semiannulus._mesh_files.push_back("assignment_semiannulus_quadrilateral.med");
  app_semiannulus._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/02_2d/annulus_semi/");
  app_semiannulus._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/02_2d/annulus_semi/");
  
  app_semiannulus._boundary_conditions_types_and_values             = semiannulus::bc_all_dirichlet;

  app_semiannulus._assemble_function_rhs = semiannulus::function_0::laplacian;
  app_semiannulus._true_solution    = semiannulus::function_0::value;
  //assignment_semiannulus - END
 

 //assignment_cube - BEGIN
  app_cube._system_name = system_common_name;
  app_cube._assemble_function                            = poisson_equation::equation_with_dirichlet_or_neumann_bc<double, double>;
  
  app_cube._mesh_files.push_back("assignment_cube_hexahedral.med");
  app_cube._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/03_3d/cube/0-1x0-1x0-1/");
  
  app_cube._boundary_conditions_types_and_values             = cube::bc_all_dirichlet;
  app_cube._assemble_function_rhs = cube::function_0::laplacian;
  app_cube._true_solution    = cube::function_0::value;
 //assignment_cube - END
  

  
  //assignment_cylinder - BEGIN
  app_cylinder._system_name = system_common_name;
  app_cylinder._assemble_function = poisson_equation::equation_with_dirichlet_or_neumann_bc<double, double>;
  
  app_cylinder._mesh_files.push_back("assignment_cylinder_tetrahedral.med");
  app_cylinder._mesh_files.push_back("assignment_cylinder_hexahedral.med");
  app_cylinder._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/03_3d/cylinder/");
  app_cylinder._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/03_3d/cylinder/");
  
  app_cylinder._boundary_conditions_types_and_values             = cylinder::bc_all_dirichlet;

  app_cylinder._assemble_function_rhs = cylinder::function_0::laplacian;
  app_cylinder._true_solution    = cylinder::function_0::value;
  //assignment_cylinder - END


  //assignment_quarter_cylinder - BEGIN
  app_quarter_cylinder._system_name = system_common_name;
  app_quarter_cylinder._assemble_function = poisson_equation::equation_with_dirichlet_or_neumann_bc<double, double>;

  app_quarter_cylinder._mesh_files.push_back("assignment_quarter_cylinder_tetrahedral.med");
  app_quarter_cylinder._mesh_files.push_back("assignment_quarter_cylinder_hexahedral.med");
//   app_quarter_cylinder._mesh_files.push_back("assignment_quarter_cylinder_hexahedral_0.med");
//   app_quarter_cylinder._mesh_files.push_back("assignment_quarter_cylinder_hexahedral_1.med");
  app_quarter_cylinder._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/03_3d/cylinder_quarter/");
  app_quarter_cylinder._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/03_3d/cylinder_quarter/");
  
  app_quarter_cylinder._boundary_conditions_types_and_values             = quarter_cylinder::bc_all_dirichlet;
  app_quarter_cylinder._assemble_function_rhs = quarter_cylinder::function_0::laplacian;
  app_quarter_cylinder._true_solution    = quarter_cylinder::function_0::value;
  //assignment_quarter_cylinder - END


  
  
  //assignment_tetra_prism_annular_base - BEGIN
  app_prism_annular_base._system_name = system_common_name;
  app_prism_annular_base._assemble_function = poisson_equation::equation_with_dirichlet_or_neumann_bc<double, double>;
  
  app_prism_annular_base._mesh_files.push_back("assignment_prism_annular_base_tetrahedral.med");
  app_prism_annular_base._mesh_files.push_back("assignment_prism_annular_base_hexahedral.med");
  app_prism_annular_base._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/03_3d/prism_annular_base/");
  app_prism_annular_base._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/03_3d/prism_annular_base/");
  
  app_prism_annular_base._boundary_conditions_types_and_values             = prism_annular_base::bc_all_dirichlet;

  app_prism_annular_base._assemble_function_rhs = prism_annular_base::function_0::laplacian;
  app_prism_annular_base._true_solution    = prism_annular_base::function_0::value;
  //assignment_tetra_prism_annular_base - END
  
 
  


  
  
  my_specifics.push_back(app_segment);
  my_specifics.push_back(app_square);
  my_specifics.push_back(app_semicircle);
  my_specifics.push_back(app_quarter_circle);
  my_specifics.push_back(app_annulus);
  my_specifics.push_back(app_semiannulus);
  my_specifics.push_back(app_cube);
  my_specifics.push_back(app_cylinder);
  my_specifics.push_back(app_quarter_cylinder);
  my_specifics.push_back(app_prism_annular_base);

    // ======= App Specifics - END  ==================
  
  
 //begin app loop  
  for (unsigned int app = 0; app < my_specifics.size(); app++)  {
        
  // ======= Problem, App (every Problem has 1 App for now) ========================
  ml_prob.set_app_specs_pointer(&my_specifics[app]);
  
  
  // ======= Mesh - BEGIN  ==================
  MultiLevelMesh ml_mesh;
  // ======= Mesh - END  ==================
 

 for (unsigned int m = 0; m < my_specifics[app]._mesh_files.size(); m++)  {
   
  // ======= Mesh, Coarse reading - BEGIN ==================
  double Lref = 1.;
  
  const bool read_groups = true; //with this being false, we don't read any group at all. Therefore, we cannot even read the boundary groups that specify what are the boundary faces, for the boundary conditions
  const bool read_boundary_groups = true;
  

  const std::string mesh_file = my_specifics[app]._mesh_files_path_relative_to_executable[m] + my_specifics[app]._mesh_files[m];
  
  ml_mesh.ReadCoarseMeshFileReadingBeforePartitioning(mesh_file.c_str(), Lref, read_groups, read_boundary_groups);
    
  ml_mesh.GetLevelZero(0)->build_dofmap_all_fe_families_and_elem_and_node_structures();
 

  ml_mesh.BuildFETypesBasedOnExistingCoarseMeshGeomElements();
  
  ml_mesh.PrepareNewLevelsForRefinement();
  // ======= Mesh, Coarse reading - END ==================



 for (unsigned int r = 1; r < 3; r++)  {

  // ======= Mesh: Refinement - BEGIN  ==================
  unsigned numberOfUniformLevels = r;
  unsigned numberOfSelectiveLevels = 0;
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  // ======= Mesh: Refinement - END  ==================
  
  // ======= Mesh: Coarse erasing - BEGIN  ========================
  ml_mesh.EraseCoarseLevels(numberOfUniformLevels + numberOfSelectiveLevels - 1);
  ml_mesh.PrintInfo();
  // ======= Mesh: Coarse erasing - END  ========================
  
  // ======= Solution - BEGIN ==================
  MultiLevelSolution ml_sol(&ml_mesh);

  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);
  
  // ======= Problem, Mesh and Solution  ==================
  ml_prob.SetMultiLevelMeshAndSolution(& ml_sol);
  // ======= Solution - END ==================
  
  // ======= Solutions that are Unknowns - BEGIN ==================

 // ======= Solution: Add ==================
  std::vector< Unknown > unknowns = provide_list_of_unknowns();   ///@todo probably this should go in the Problem, to be accessed from other places
  
   for (unsigned int u = 0; u < unknowns.size(); u++)  { 
      ml_sol.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);
  }
 // ======= Solution: Initial Conditions ==================
   for (unsigned int u = 0; u < unknowns.size(); u++)  { 
      ml_sol.Initialize(unknowns[u]._name.c_str(), InitialValueU, & ml_prob);
  }
  // ======= Solution: Boundary Conditions ==================
  ml_sol.AttachSetBoundaryConditionFunction(my_specifics[app]._boundary_conditions_types_and_values);
   for (unsigned int u = 0; u < unknowns.size(); u++)  { 
      ml_sol.GenerateBdc(unknowns[u]._name.c_str(), (unknowns[u]._time_order == 0) ? "Steady" : "Time_dependent", & ml_prob);
  }

  // ======= Solutions that are Unknowns - END ==================

  

//   std::vector < std::vector < const elem_type_templ_base<double, double> *  > > elem_all = ml_prob.evaluate_all_fe<double, double>();
  
    // ======= Problem, System - BEGIN ========================
  ml_prob.clear_systems();
  
  NonLinearImplicitSystem & system = ml_prob.add_system < NonLinearImplicitSystem > (my_specifics[app]._system_name);
  
  system.SetDebugNonlinear(true);
 
       // ======= System Unknowns ========================
  for (unsigned int u = 0; u < unknowns.size(); u++)  system.AddSolutionToSystemPDE(unknowns[u]._name.c_str());  
 
  // attach the assembling function to system
  system.SetAssembleFunction( my_specifics[app]._assemble_function );

//   system.SetMaxNumberOfLinearIterations(2);
  // initialize and solve the system
  system.SetMgType(V_CYCLE/*F_CYCLE*//*M_CYCLE*/); //it doesn't matter if I use only 1 level


  system.SetOuterSolver(GMRES);
 
  system.init();
  
  system.MGsolve();
  
  compute_L2_norm_of_errors_of_unknowns_with_analytical_sol<double, double>(ml_prob);
    // ======= Problem, System - END ========================
  
  
    // ======= Print - BEGIN ========================
  const std::string print_order = "biquadratic"; //"linear", "quadratic", "biquadratic"
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
 
  ml_sol.GetWriter()->Write(my_specifics[app]._system_name + "_" + my_specifics[app]._mesh_files/*mesh_files*/[m], files.GetOutputPath(), print_order.c_str(), variablesToBePrinted);
    // ======= Print - END ========================
  
    }  //end refinement loop
  
  }  //end mesh file loop
 
 
} //end app loop
 
  return 0;
}











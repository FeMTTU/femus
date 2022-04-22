/**
 * Solve 
 *     - \Delta u = 1
*/


#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "NumericVector.hpp"

#include "CurrentElem.hpp"
#include "ElemType_template.hpp"

#include "Assemble_jacobian.hpp"
#include "Assemble_unknown_jacres.hpp"

#include "app_specifics.hpp"


using namespace femus;
 

bool semiannulus__laplacian__bc(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

  bool dirichlet = false;
  value = 0.;
  
  const double tolerance = 1.e-5;
  
 
 if (ml_prob->GetMLMesh()->GetDimension() == 2 )  {
     
     
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
        value = 0. * ( x[0] * x[0]); //Neumann function, here we specify the WHOLE normal derivative, which is a scalar, not each Cartesian component
  }
   
 
 }
 
 
 
  return dirichlet;
  
 }


double semiannulus__laplacian__rhs(const std::vector < double >& x) {
  
    double r2 = x[0] * x[0] + x[1] * x[1];
    double temp = x[0] * (8. - 4.5 / (sqrt(r2)));
    //double temp = (4. - 1.5 / (sqrt(r2)));
  return temp;

//   return - 20. * x[0] * (-0.45 + x[0] * x[0] - 0.6 * x[1] * x[1] );
    
}

double semiannulus__laplacian__true_solution(const std::vector < double >& x) {
    
     double r2 = x[0] * x[0] + x[1] * x[1];

     return   x[0] * (sqrt(r2) - 1. ) * (.5 - sqrt(r2) );
//      return   x[0] * (r2 - 1. ) * (.25 - r2);
}


// (z)(z-2)((x-1)^2 + (y-1)^2 - 1) 
 double cylinder__laplacian__true_solution(const std::vector<double> & x_qp){
  
     double r = 4*x_qp[2]*(x_qp[2] - 2)  +  2*((x_qp[0] - 1)*(x_qp[0] - 1)  +  (x_qp[1] - 1)*(x_qp[1] - 1) - 1);
  return -r;

     
}
 
double  cylinder__laplacian__rhs(const std::vector < double >& x_qp) {
   double r = 4*x_qp[2]*(x_qp[2] - 2)  +  2*((x_qp[0] - 1)*(x_qp[0] - 1)  +  (x_qp[1] - 1)*(x_qp[1] - 1) - 1);
  return -r;
}


 
bool cylinder__laplacian__bc(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

  bool dirichlet = false;
  value = 0.;
  
  const double tolerance = 1.e-5;
  
  
  if (ml_prob->GetMLMesh()->GetDimension() == 3 )  {
     
     
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
   
 
 }
 
 
  return dirichlet;
  
 }

 double quarter_circle__laplacian__true_solution(const std::vector<double> & x_qp){
    
    // for a quarter-circle in Quadrant 1
    
    double x = x_qp[0];
    double y = x_qp[1];
    
    return  x * y * (1.0 - (x*x + y*y)); // forced to be zero on the x and y axis, and the circle edge
}

 

// this is specifically the laplacian of the function given above
// flynn, user-made equation - accepts only coordinates
double quarter_circle__laplacian__rhs(const std::vector<double> & x_qp){
    
    // for a quarter-circle in Quadrant 1
    
    double x = x_qp[0];
    double y = x_qp[1];
    
    return -1.0 * -12. * x * y;
}



bool quarter_circle__laplacian__bc(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

  bool dirichlet = false;
  value = 0.;
  
  const double tolerance = 1.e-5;
  
  
 if (ml_prob->GetMLMesh()->GetDimension() == 2 )  {
     
     
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


 }

 
  return dirichlet;
  
 }




bool prism_annular_base__laplacian__bc(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

    bool dirichlet = false;
    value = 0.;

    const double tolerance = 1.e-5;


    if (ml_prob->GetMLMesh()->GetDimension() == 3 )  {


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

    }
    return dirichlet;

}


 double prism_annular_base__laplacian__true_solution(const std::vector<double> & x) {
     
     return  x[2] * (x[2] - 1.) * (1. - x[0]*x[0] - x[1]*x[1]) * (0.25 - x[0]*x[0] - x[1]*x[1]);
 }
 

//calculator: f = z(z - 1)(1 - x^2 - y^2)(1/4 - x^2 - y^2);
double prism_annular_base__laplacian__rhs(const std::vector < double > & x) {
    
  double r2 = 0.5 - 2.5*pow(x[0],2) + 2*pow(x[0],4) - 2.5*pow(x[1],2) + 4*pow(x[0],2)*pow(x[1],2) + 2*pow(x[1],4) + 5.*x[2] - 16*pow(x[0],2)*x[2] - 16*pow(x[1],2)*x[2] - 5.*pow(x[2],2) + 16*pow(x[0],2)*pow(x[2],2) + 16*pow(x[1],2)*pow(x[2],2);
  
   return - r2;
   
  }




double segment_dir_neu_fine__laplacian__true_solution(const std::vector<double> & x) {
    
    // for a 1d segment
    
    return  x[0] * (1. - x[0]);
}


// user-made equation - accepts only coordinates
double segment_dir_neu_fine__laplacian__rhs(const std::vector<double> & x_qp){
    
    // for a 1d segment
    
    return  2.;
}



 
 
bool segment_dir_neu_fine__laplacian__bc(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

  bool dirichlet = false;
  value = 0.;
  
  const double tolerance = 1.e-5;
  
 if (ml_prob->GetMLMesh()->GetDimension() == 1 )  {
  
  if (face_name == 1) {
      dirichlet = true;
        value = 0.; //Dirichlet value
    }
  else if (face_name == 2) {
      dirichlet = true;
        value = 0.; //Dirichlet value
    }

    
 }
 
 
 
  return dirichlet;
  
 }

 
double InitialValueU(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {
    
  return 0.;
  
}

 
 //====== NEUMANN LOOP 1D =============================================   
void laplacian_natural_loop_1d(const MultiLevelProblem *    ml_prob, 
                     const Mesh *                    msh,
                     const MultiLevelSolution *    ml_sol, 
                     const unsigned iel,
                     CurrentElem < double > & geom_element,
                     const unsigned xType,
                     const std::string solname_u,
                     const unsigned solFEType_u,
                     std::vector< double > & Res
                    ) {
    
     double grad_u_dot_n;
    
    for (unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
        
       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, xType);
       
       geom_element.set_elem_center_bdry_3d();
       
       std::vector  <  double > xx_face_elem_center(3, 0.); 
          xx_face_elem_center = geom_element.get_elem_center_bdry_3d();
        
       const int boundary_index = msh->el->GetFaceElementIndex(iel, jface);
       
       if ( boundary_index < 0) { //I am on the boundary
                  
         unsigned int face = - (boundary_index + 1);
    
         bool is_dirichlet =  ml_sol->GetBdcFunctionMLProb()(ml_prob, xx_face_elem_center, solname_u.c_str(), grad_u_dot_n, face, 0.);                     
         //we have to be careful here, because in GenerateBdc those coordinates are passed as NODE coordinates, 
         //while here we pass the FACE ELEMENT CENTER coordinates. 
         // So, if we use this for enforcing space-dependent Dirichlet or Neumann values, we need to be careful!
         
             if ( !(is_dirichlet)  &&  (grad_u_dot_n != 0.) ) {  //dirichlet == false and nonhomogeneous Neumann
                 
                 
                 
                   unsigned n_dofs_face = msh->GetElementFaceDofNumber(iel, jface, solFEType_u);

                  for (unsigned i = 0; i < n_dofs_face; i++) {
                      
                 unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);

                 Res[i_vol] +=  grad_u_dot_n /* * phi[node] = 1. */;
                 
                         }
                         
                         
                         
        
                    }
                  
              }

    }
    
}




template < class real_num, class real_num_mov >
void laplacian_natural_loop_2d3d(const MultiLevelProblem *    ml_prob, 
                       const Mesh *                    msh,
                       const MultiLevelSolution *    ml_sol, 
                       const unsigned iel,
                       CurrentElem < double > & geom_element,
                       const unsigned solType_coords,
                       const std::string solname_u,
                       const unsigned solFEType_u,
                       std::vector< double > & Res,
                       //-----------
                       std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > >  elem_all,
                       const unsigned dim,
                       const unsigned space_dim,
                       const unsigned max_size
                    ) {
    
    
    /// @todo - should put these outside the iel loop --
    std::vector < std::vector < double > >  JacI_iqp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_iqp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_iqp_bdry.size(); d++) {   Jac_iqp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_iqp_bdry.size(); d++) { JacI_iqp_bdry[d].resize(dim-1); }
    
    double detJac_iqp_bdry;
  double weight_iqp_bdry = 0.;
// ---
  //boundary state shape functions
  vector <double> phi_u_bdry;  
  vector <double> phi_u_x_bdry; 

  phi_u_bdry.reserve(max_size);
  phi_u_x_bdry.reserve(max_size * space_dim);
// ---
  

    
    
     double grad_u_dot_n;
    
    for (unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
        
       geom_element.set_coords_at_dofs_bdry_3d(iel, jface, solType_coords);
       
       geom_element.set_elem_center_bdry_3d();

       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, jface);    
       

       std::vector <  double > xx_face_elem_center(3, 0.); 
       xx_face_elem_center = geom_element.get_elem_center_bdry_3d();
        
       const int boundary_index = msh->el->GetFaceElementIndex(iel, jface);
       
       if ( boundary_index < 0) { //I am on the boundary
                  
         unsigned int face = - (boundary_index + 1);
    
         bool is_dirichlet =  ml_sol->GetBdcFunctionMLProb()(ml_prob, xx_face_elem_center, solname_u.c_str(), grad_u_dot_n, face, 0.);                     
         //we have to be careful here, because in GenerateBdc those coordinates are passed as NODE coordinates, 
         //while here we pass the FACE ELEMENT CENTER coordinates. 
         // So, if we use this for enforcing space-dependent Dirichlet or Neumann values, we need to be careful!
         
             if ( !(is_dirichlet) /* &&  (grad_u_dot_n != 0.)*/ ) {  //dirichlet == false and nonhomogeneous Neumann
                 
                 
                 
                        const unsigned n_gauss_bdry = ml_prob->GetQuadratureRule(ielGeom_bdry).GetGaussPointsNumber();
        
    
		for(unsigned ig_bdry = 0; ig_bdry < n_gauss_bdry; ig_bdry++) {
  
     elem_all[ielGeom_bdry][solType_coords]->JacJacInv(geom_element.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_iqp_bdry, JacI_iqp_bdry, detJac_iqp_bdry, space_dim);
//      elem_all[ielGeom_bdry][solType_coords]->compute_normal(Jac_iqp_bdry, normal);
    
    weight_iqp_bdry = detJac_iqp_bdry * ml_prob->GetQuadratureRule(ielGeom_bdry).GetGaussWeightsPointer()[ig_bdry];
                
    elem_all[ielGeom_bdry][solFEType_u ]->shape_funcs_current_elem(ig_bdry, JacI_iqp_bdry, phi_u_bdry, phi_u_x_bdry,  boost::none, space_dim);
                 
                   unsigned n_dofs_face = msh->GetElementFaceDofNumber(iel, jface, solFEType_u);

                  for (unsigned i_bdry = 0; i_bdry < n_dofs_face; i_bdry++) {
                      
                 unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

                 Res[i_vol] +=  weight_iqp_bdry * grad_u_dot_n * phi_u_bdry[i_bdry];
                 
                           }
                         
                         
                         
                        }
        
                         
        
                    }
                  
              }
    }
    
}






 

template < class real_num, class real_num_mov >
void laplacian_dir_neu_eqn(MultiLevelProblem& ml_prob);

template < class real_num, class real_num_mov >
void compute_norm(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Files ========================
  const bool use_output_time_folder = false;
  const bool redirect_cout_to_file = false;
  Files files; 
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);

  // ======= Quad Rule ========================
  std::string fe_quad_rule("seventh");

    // ======= App Specifics  ==================
  std::vector< app_specifics >   my_specifics;
  
  app_specifics  app_segment;
  app_specifics  app_prism_annular_base;
  app_specifics  app_quarter_circle;
  app_specifics  app_cylinder;
  app_specifics  app_semiannulus;
  
  
  //segment_dir_neu_fine
  app_segment._mesh_files.push_back("assignment_segment_dir_neu_fine.med");
  app_segment._mesh_files.push_back("assignment_segment_dir_neu_fine.med");
  
  app_segment._assemble_function = laplacian_dir_neu_eqn<double, double>;
  app_segment._assemble_function_natural_boundary_loop_1d = laplacian_natural_loop_1d;
  app_segment._assemble_function_natural_boundary_loop_2d3d = laplacian_natural_loop_2d3d;
  app_segment._assemble_function_rhs = segment_dir_neu_fine__laplacian__rhs;
  app_segment._bdry_func = segment_dir_neu_fine__laplacian__bc;
  
  app_segment._norm_true_solution = segment_dir_neu_fine__laplacian__true_solution;
  
  //assignment_tetra_prism_annular_base
  app_prism_annular_base._mesh_files.push_back("assignment_prism_annular_base_tetrahedral.med");
  app_prism_annular_base._mesh_files.push_back("assignment_prism_annular_base_hexahedral.med");
  
  app_prism_annular_base._assemble_function = laplacian_dir_neu_eqn<double, double>;
  app_prism_annular_base._assemble_function_natural_boundary_loop_1d = laplacian_natural_loop_1d;
  app_prism_annular_base._assemble_function_natural_boundary_loop_2d3d = laplacian_natural_loop_2d3d;
  app_prism_annular_base._assemble_function_rhs = prism_annular_base__laplacian__rhs;
  app_prism_annular_base._bdry_func = prism_annular_base__laplacian__bc;

  app_prism_annular_base._norm_true_solution = prism_annular_base__laplacian__true_solution;
  
  //assignment_quarter_circle
  app_quarter_circle._mesh_files.push_back("assignment_quarter_circle_triangular.med");
  app_quarter_circle._mesh_files.push_back("assignment_quarter_circle_quadrangular.med");
  
  app_quarter_circle._assemble_function = laplacian_dir_neu_eqn<double, double>;
  app_quarter_circle._assemble_function_natural_boundary_loop_1d = laplacian_natural_loop_1d;
  app_quarter_circle._assemble_function_natural_boundary_loop_2d3d = laplacian_natural_loop_2d3d;
  app_quarter_circle._assemble_function_rhs = quarter_circle__laplacian__rhs;
  app_quarter_circle._bdry_func = quarter_circle__laplacian__bc;

  app_quarter_circle._norm_true_solution = quarter_circle__laplacian__true_solution;
  
  //assignment_cylinder
  app_cylinder._mesh_files.push_back("assignment_cylinder_tetrahedral.med");
  app_cylinder._mesh_files.push_back("assignment_cylinder_hexahedral.med");
  
  app_cylinder._assemble_function = laplacian_dir_neu_eqn<double, double>;
  app_cylinder._assemble_function_natural_boundary_loop_1d = laplacian_natural_loop_1d;
  app_cylinder._assemble_function_natural_boundary_loop_2d3d = laplacian_natural_loop_2d3d;
  app_cylinder._assemble_function_rhs = cylinder__laplacian__rhs;
  app_cylinder._bdry_func = cylinder__laplacian__bc;

  app_cylinder._norm_true_solution = cylinder__laplacian__true_solution;
 
  //assignment_semiannulus
  app_semiannulus._mesh_files.push_back("assignment_semiannulus_triangular.med");
  app_semiannulus._mesh_files.push_back("assignment_semiannulus_quadrangular.med");
  
  app_semiannulus._assemble_function = laplacian_dir_neu_eqn<double, double>;
  app_semiannulus._assemble_function_natural_boundary_loop_1d = laplacian_natural_loop_1d;
  app_semiannulus._assemble_function_natural_boundary_loop_2d3d = laplacian_natural_loop_2d3d;
  
  app_semiannulus._assemble_function_rhs = semiannulus__laplacian__rhs;
  app_semiannulus._bdry_func = semiannulus__laplacian__bc;
  app_semiannulus._norm_true_solution = semiannulus__laplacian__true_solution;
 

//   //assignment_annulus - jon
//   my_specifics[2]._mesh_files.push_back("assignment_quarter_circle_triangular.med");
//   my_specifics[2]._mesh_files.push_back("assignment_quarter_circle_quadrangular.med");
//   
//   my_specifics[2]._assemble_function = laplacian_dir_neu_eqn<double, double>;
//   my_specifics[2]._assemble_function_natural_boundary_loop_1d = laplacian_natural_loop_1d;
//   my_specifics[2]._assemble_function_natural_boundary_loop_2d3d = laplacian_natural_loop_2d3d;
//   my_specifics[2]._assemble_function_rhs = quarter_circle__laplacian__rhs;
//   my_specifics[2]._bdry_func = quarter_circle__laplacian__bc;

  
//   my_specifics.push_back(app_segment);
//   my_specifics.push_back(app_quarter_circle);
//   my_specifics.push_back(app_prism_annular_base);
//   my_specifics.push_back(app_cylinder);
  my_specifics.push_back(app_semiannulus);
  
  
  
  for (unsigned int app = 0; app < my_specifics.size(); app++)  { //begin app loop
      
   // ======= Mesh  ==================
//    std::vector<std::string> mesh_files;
//    
//    mesh_files.push_back(my_specifics[app]._mesh_files[0]);
//    mesh_files.push_back(my_specifics[app]._mesh_files[1]);
//    mesh_files.push_back("Mesh_2_xy_boundaries_groups_4x4.med");
//    mesh_files.push_back("Mesh_1_x_all_dir.med");
//    mesh_files.push_back("Mesh_1_y_all_dir.med");
//    mesh_files.push_back("Mesh_1_z_all_dir.med");
//    mesh_files.push_back("Mesh_2_xz_all_dir.med");
//    mesh_files.push_back("Mesh_2_yz_all_dir.med");
//    mesh_files.push_back("Mesh_3_xyz_all_dir.med");
//    mesh_files.push_back("dome_tri.med");
//    mesh_files.push_back("dome_quad.med");
//    mesh_files.push_back("disk_quad.med");
//    mesh_files.push_back("disk_quad_45x.med");
//    mesh_files.push_back("disk_quad_90x.med");
//    mesh_files.push_back("disk_tri.med");
//    mesh_files.push_back("disk_tri_45x.med");
//    mesh_files.push_back("disk_tri_90x.med");
   


 for (unsigned int m = 0; m < my_specifics[app]._mesh_files/*mesh_files*/.size(); m++)  {
   
  // ======= Mesh  ==================
  // define multilevel mesh
  MultiLevelMesh ml_mesh;
  double scalingFactor = 1.;
 
  const bool read_groups = true; //with this being false, we don't read any group at all. Therefore, we cannot even read the boundary groups that specify what are the boundary faces, for the boundary conditions
  const bool read_boundary_groups = true;
  
  std::string mesh_file_tot = "./input/" + my_specifics[app]._mesh_files/*mesh_files*/[m];
  
  ml_mesh.ReadCoarseMesh(mesh_file_tot.c_str(), fe_quad_rule.c_str(), scalingFactor, read_groups, read_boundary_groups);
//     ml_mesh.GenerateCoarseBoxMesh(2,0,0,0.,1.,0.,0.,0.,0.,EDGE3,fe_quad_rule.c_str());
//     ml_mesh.GenerateCoarseBoxMesh(0,2,0,0.,0.,0.,1.,0.,0.,EDGE3,fe_quad_rule.c_str());
 
  unsigned numberOfUniformLevels = 6;
  unsigned numberOfSelectiveLevels = 0;
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
  ml_mesh.EraseCoarseLevels(numberOfUniformLevels + numberOfSelectiveLevels - 1);
  ml_mesh.PrintInfo();

  
  
  // ======= Solution  ==================
  MultiLevelSolution ml_sol(&ml_mesh);

  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);
  
  // ======= Problem ========================
  // define the multilevel problem attach the ml_sol object to it
  MultiLevelProblem ml_prob(&ml_sol);
  
  // add variables to ml_sol
  ml_sol.AddSolution("u", LAGRANGE, SECOND/*DISCONTINUOUS_POLYNOMIAL, ZERO*/);
  
  // ======= Solution: Initial Conditions ==================
  ml_sol.Initialize("All");    // initialize all variables to zero
  ml_sol.Initialize("u", InitialValueU, & ml_prob);

  // ======= Solution: Boundary Conditions ==================
  ml_sol.AttachSetBoundaryConditionFunction(my_specifics[app]._bdry_func);
  ml_sol.GenerateBdc("u", "Steady",  & ml_prob);

  

  // ======= Problem, II ========================
  ml_prob.SetFilesHandler(&files);
  ml_prob.set_app_specs_pointer(&my_specifics[app]);
  ml_prob.SetQuadratureRuleAllGeomElems(fe_quad_rule);
  ml_prob.set_all_abstract_fe_multiple();
  
//   std::vector < std::vector < const elem_type_templ_base<double, double> *  > > elem_all = ml_prob.evaluate_all_fe<double, double>();
  
    // ======= System ========================
 // add system  in ml_prob as a Linear Implicit System
  NonLinearImplicitSystem& system = ml_prob.add_system < NonLinearImplicitSystem > ("Laplace");
  
  system.SetDebugNonlinear(true);
 
  system.AddSolutionToSystemPDE("u");
 
  // attach the assembling function to system
  system.SetAssembleFunction( my_specifics[app]._assemble_function );

//   system.SetMaxNumberOfLinearIterations(2);
  // initialize and solve the system
  system.SetMgType(V_CYCLE/*F_CYCLE*//*M_CYCLE*/); //it doesn't matter if I use only 1 level


  system.SetOuterSolver(GMRES);
 
  system.init();
  
  system.MGsolve();
  
  compute_norm<double, double>(ml_prob);
  
    // ======= Print ========================
  // print solutions
  const std::string print_order = "biquadratic"; //"linear", "quadratic", "biquadratic"
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("all");
 
  ml_sol.GetWriter()->Write(my_specifics[app]._mesh_files/*mesh_files*/[m], files.GetOutputPath(), print_order.c_str(), variablesToBePrinted);
  
  }
 
 
} //end app loop
 
  return 0;
}





template < class real_num, class real_num_mov >
void laplacian_dir_neu_eqn(MultiLevelProblem& ml_prob) {

  NonLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystem> ("Laplace");  
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
  SparseMatrix*             JAC = pdeSys->_KK;
  NumericVector*           RES = pdeSys->_RES;

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id();

  //=============== Geometry ========================================
  unsigned xType = BIQUADR_FE; // the FE for the domain need not be biquadratic
  
  CurrentElem < double > geom_element(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
//***************************************************  


 //******************** quadrature *******************************  
  double jacXweight_qp; 

 //********************* unknowns *********************** 
 //***************************************************  
  const int n_vars = mlPdeSys->GetSolPdeIndex().size();
  std::cout << "************" << n_vars << "************";
  std::vector <double> phi_u;
  std::vector <double> phi_u_x; 
  std::vector <double> phi_u_xx;

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * space_dim);
  phi_u_xx.reserve(maxSize * dim2);
  
  const std::string solname_u = "u";
  unsigned solIndex_u;
  solIndex_u = ml_sol->GetIndex(solname_u.c_str()); 
  unsigned solFEType_u = ml_sol->GetSolutionType(solIndex_u); 

  unsigned solPdeIndex_u;
  solPdeIndex_u = mlPdeSys->GetSolPdeIndex(solname_u.c_str());

  std::vector < double >  sol_u;     sol_u.reserve(maxSize);
  std::vector< int > l2GMap_u;    l2GMap_u.reserve(maxSize);
 //***************************************************  
 //***************************************************  

  
 //***************************************************  
 //********* WHOLE SET OF VARIABLES ****************** 

  std::vector< int > l2GMap_AllVars; l2GMap_AllVars.reserve(n_vars*maxSize); // local to global mapping
  std::vector< double >         Res;            Res.reserve(n_vars*maxSize);  // local redidual vector
  std::vector < double >        Jac;            Jac.reserve(n_vars*maxSize * n_vars*maxSize);
 //***************************************************  

  RES->zero();
  if (assembleMatrix)  JAC->zero();

  
 //***************************************************  
     std::vector < std::vector < real_num_mov > >  JacI_qp(space_dim);
     std::vector < std::vector < real_num_mov > >  Jac_qp(dim);
    for (unsigned d = 0; d < dim; d++) { Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < space_dim; d++) { JacI_qp[d].resize(dim); }
    
    real_num_mov detJac_qp;
    

  //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
 //***************************************************  
  
  

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
      

    geom_element.set_coords_at_dofs_and_geom_type(iel, xType);
        
    const short unsigned ielGeom = geom_element.geom_type();

 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solFEType_u);
    sol_u    .resize(nDof_u);
    l2GMap_u.resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solFEType_u);
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
      l2GMap_u[i] = pdeSys->GetSystemDof(solIndex_u, solPdeIndex_u, i, iel);
    }
 //***************************************************  
 
 //******************** ALL VARS ********************* 
    unsigned nDof_AllVars = nDof_u; 
    int nDof_max    =  nDof_u;   // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
    Res.resize(nDof_AllVars);                  std::fill(Res.begin(), Res.end(), 0.);
    Jac.resize(nDof_AllVars * nDof_AllVars);   std::fill(Jac.begin(), Jac.end(), 0.);
    l2GMap_AllVars.resize(0);                  l2GMap_AllVars.insert(l2GMap_AllVars.end(),l2GMap_u.begin(),l2GMap_u.end());
 //*************************************************** 
    

 //========= BOUNDARY ==================   
    if (dim == 1)   ml_prob.get_app_specs_pointer()->_assemble_function_natural_boundary_loop_1d(& ml_prob, msh, ml_sol,
                      iel, geom_element, xType,
                      solname_u, solFEType_u,
                      Res
                     );

    if (dim == 2 || dim == 3)   ml_prob.get_app_specs_pointer()->_assemble_function_natural_boundary_loop_2d3d(& ml_prob, msh, ml_sol,
                      iel, geom_element, xType,
                      solname_u, solFEType_u,
                      Res,
                      elem_all,
                      dim,
                      space_dim,
                      maxSize
                     );
 
 //========= VOLUME ==================   
   
 //========= gauss value quantities ==================   
	std::vector<double> sol_u_x_gss(space_dim);     std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
 //===================================================   
    
    //--------------    from giorgiobornia.cpp
 /// @assignment You need to evaluate your manufactured right hand side at the quadrature point qp.
 /// Hence, you need to compute the coordinates of the quadrature point. Let us call them x_qp.
 /// These are obtained just like every quantity at a quadrature point, i.e., by interpolating the values of the quantity at the element nodes.
 /// The interpolation is performed by using the shape functions.
 /// In other words, 
 ///         (x_qp) = summation of (x_nodes) * (shape function of that node, evaluated at qp)
 /// 
 ///   (x_nodes) are obtained from   geom_element.get_coords_at_dofs_3d()  (this is a  vector< vector >,  where the outer index is the dimension and the inner index ranges over the nodes) 
 ///   (shape function of that node, evaluated at qp)  is obtained from phi_u  (this is a vector, whose index ranges over the nodes)
 
 
//--------------   
    
    
    
      // *** Quadrature point loop ***
      for (unsigned i_qp = 0; i_qp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); i_qp++) {
          
        // *** get gauss point weight, test function and test function partial derivatives ***
// 	msh->_finiteElement[ielGeom][solFEType_u]->Jacobian(geom_element.get_coords_at_dofs_3d(),    i_qp, weight,    phi_u,    phi_u_x,    boost::none /*phi_u_xx*/);
          
	elem_all[ielGeom][xType]->JacJacInv(geom_element.get_coords_at_dofs_3d(), i_qp, Jac_qp, JacI_qp, detJac_qp, space_dim);
    jacXweight_qp = detJac_qp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[i_qp];
    elem_all[ielGeom][solFEType_u]->shape_funcs_current_elem(i_qp, JacI_qp, phi_u, phi_u_x, boost::none /*phi_u_xx*/, space_dim);


//--------------    
	std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < nDof_u; i++) {
// 	                                                sol_u_gss      += sol_u[i] * phi_u[i];
                   for (unsigned d = 0; d < sol_u_x_gss.size(); d++)   sol_u_x_gss[d] += sol_u[i] * phi_u_x[i * space_dim + d];
          }
//--------------    

//--------------    
 /// @assignment You need to evaluate your manufactured right hand side at the quadrature point qp.
 /// Hence, you need to compute the coordinates of the quadrature point. Let us call them x_qp.
 /// These are obtained just like every quantity at a quadrature point, i.e., by interpolating the values of the quantity at the element nodes.
 /// The interpolation is performed by using the shape functions.
 /// In other words, 
 ///         (x_qp) = summation of (x_nodes) * (shape function of that node, evaluated at qp)
 /// 
 ///   (x_nodes) are obtained from   geom_element.get_coords_at_dofs_3d()  (this is a  vector< vector >,  where the outer index is the dimension and the inner index ranges over the nodes) 
 ///   (shape function of that node, evaluated at qp)  is obtained from phi_u  (this is a vector, whose index ranges over the nodes)
 
 std::vector<double> x_qp(dim, 0.);
          
        for (unsigned i = 0; i < nDof_u; i++) {
          	for (unsigned d = 0; d < dim; d++) {
	                                                x_qp[d]    += geom_element.get_coords_at_dofs_3d()[d][i] * phi_u[i]; // fetch of coordinate points
            }
        }
 

          
//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
        for (unsigned i = 0; i < nDof_max; i++) {
	  
//--------------    
	      double laplace_res_du_u_i = 0.;
              if ( i < nDof_u ) {
                  for (unsigned kdim = 0; kdim < space_dim; kdim++) {
                       laplace_res_du_u_i             +=  phi_u_x   [i * space_dim + kdim] * sol_u_x_gss[kdim];
	             }
              }
//--------------    
              
	      
//======================Residuals=======================
          // FIRST ROW
 /// @assignment for your manufactured right-hand side, implement a function that receives the coordinate of the quadrature point
 /// Put it after the includes, in the top part of this file
 if (i < nDof_u)                      Res[0      + i] +=  jacXweight_qp * ( phi_u[i] * ( ml_prob.get_app_specs_pointer()->_assemble_function_rhs(x_qp)  ) - laplace_res_du_u_i);
//           if (i < nDof_u)                      Res[0      + i] += jacXweight_qp * ( phi_u[i] * (  1. ) - laplace_beltrami_res_du_u_i);
//======================Residuals=======================
	      
          if (assembleMatrix) {
	    
            // *** phi_j loop ***
            for (unsigned j = 0; j < nDof_max; j++) {


//--------------    
              double laplace_mat_du_u_i_j = 0.;

              if ( i < nDof_u && j < nDof_u ) {
                  for (unsigned kdim = 0; kdim < space_dim; kdim++) {
                         laplace_mat_du_u_i_j           += phi_u_x   [i * space_dim + kdim] *
                                                           phi_u_x   [j * space_dim + kdim];
	              }
             }
//--------------    



              //============ delta_state row ============================
              //DIAG BLOCK delta_state - state
		  if ( i < nDof_u && j < nDof_u )       Jac[ (0 + i) * nDof_AllVars   + 	(0 + j) ]  += jacXweight_qp * laplace_mat_du_u_i_j;
// 		  if ( i < nDof_u && j < nDof_u )       Jac[ (0 + i) * nDof_AllVars   + 	(0 + j) ]  += jacXweight_qp * laplace_beltrami_mat_du_u_i_j; ///@todo On a flat domain, this must coincide with the standard Laplacian, so we can do a double check with this
            } // end phi_j loop
          } // endif assemble_matrix

        } // end phi_i loop
        
      } // end gauss point loop


    RES->add_vector_blocked(Res, l2GMap_AllVars);

    if (assembleMatrix) {
      JAC->add_matrix_blocked(Jac, l2GMap_AllVars, l2GMap_AllVars);
    }
   
   
  } //end element loop for each process

  RES->close();

  if (assembleMatrix) JAC->close();

     //print JAC and RES to files
    const unsigned nonlin_iter = 0/*mlPdeSys->GetNonlinearIt()*/;
    assemble_jacobian< double, double >::print_global_jacobian(assembleMatrix, ml_prob, JAC, nonlin_iter);
    assemble_jacobian< double, double >::print_global_residual(ml_prob, RES, nonlin_iter);
  

  // ***************** END ASSEMBLY *******************

  return;
}



template < class real_num, class real_num_mov >
void compute_norm(MultiLevelProblem& ml_prob) {

  NonLinearImplicitSystem* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystem> ("Laplace");  
  const unsigned level = mlPdeSys->GetLevelToAssemble();
  const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->el;

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id();

  //=============== Geometry ========================================
  unsigned xType = BIQUADR_FE; // the FE for the domain need not be biquadratic
  
  CurrentElem < double > geom_element(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
//***************************************************  


 //******************** quadrature *******************************  
  double jacXweight_qp; 

 //********************* unknowns *********************** 
 //***************************************************  
  const int n_vars = mlPdeSys->GetSolPdeIndex().size();
  std::cout << "************" << n_vars << "************";
  std::vector <double> phi_u;
  std::vector <double> phi_u_x; 
  std::vector <double> phi_u_xx;

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * space_dim);
  phi_u_xx.reserve(maxSize * dim2);
  
  const std::string solname_u = "u";
  unsigned solIndex_u;
  solIndex_u = ml_sol->GetIndex(solname_u.c_str()); 
  unsigned solFEType_u = ml_sol->GetSolutionType(solIndex_u); 

  unsigned solPdeIndex_u;
  solPdeIndex_u = mlPdeSys->GetSolPdeIndex(solname_u.c_str());

  std::vector < double >  sol_u;     sol_u.reserve(maxSize);
  std::vector < double >  sol_u_true;     sol_u_true.reserve(maxSize);
 //***************************************************  
 //***************************************************  


  
 //***************************************************  
     std::vector < std::vector < real_num_mov > >  JacI_qp(space_dim);
     std::vector < std::vector < real_num_mov > >  Jac_qp(dim);
    for (unsigned d = 0; d < dim; d++) { Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < space_dim; d++) { JacI_qp[d].resize(dim); }
    
    real_num_mov detJac_qp;
    

  //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, real_num_mov> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
 //***************************************************  
  
  double  norm = 0.;
  

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
      

    geom_element.set_coords_at_dofs_and_geom_type(iel, xType);
        
    const short unsigned ielGeom = geom_element.geom_type();

 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solFEType_u);
    sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solFEType_u);
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
    }
    
 //**************** true sol **************************** 
    sol_u_true    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
        std::vector< double > xyz_i(dim);
             	for (unsigned d = 0; d < xyz_i.size(); d++) {
                   xyz_i[d] = geom_element.get_coords_at_dofs_3d()[d][i];
                }
               sol_u_true[i]  =   ml_prob.get_app_specs_pointer()->_norm_true_solution( xyz_i);  
    }

    //***************************************************  
 
 //******************** ALL VARS ********************* 
    unsigned nDof_AllVars = nDof_u; 
    int nDof_max    =  nDof_u;   // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
    
 //*************************************************** 
    

// // //  //========= BOUNDARY ==================   
// // //     if (dim == 1)   ml_prob.get_app_specs_pointer()->_assemble_function_natural_boundary_loop_1d(& ml_prob, msh, ml_sol,
// // //                       iel, geom_element, xType,
// // //                       solname_u, solFEType_u,
// // //                       Res
// // //                      );
// // // 
// // //     if (dim == 2 || dim == 3)   ml_prob.get_app_specs_pointer()->_assemble_function_natural_boundary_loop_2d3d(& ml_prob, msh, ml_sol,
// // //                       iel, geom_element, xType,
// // //                       solname_u, solFEType_u,
// // //                       Res,
// // //                       elem_all,
// // //                       dim,
// // //                       space_dim,
// // //                       maxSize
// // //                      );
 
 //========= VOLUME ==================   
   
 //========= gauss value quantities ==================   
	std::vector<double> sol_u_x_gss(space_dim);     std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
 //===================================================   
    
    
    
      // *** Quadrature point loop ***
      for (unsigned i_qp = 0; i_qp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); i_qp++) {
          
        // *** get gauss point weight, test function and test function partial derivatives ***
// 	msh->_finiteElement[ielGeom][solFEType_u]->Jacobian(geom_element.get_coords_at_dofs_3d(),    i_qp, weight,    phi_u,    phi_u_x,    boost::none /*phi_u_xx*/);
          
	elem_all[ielGeom][xType]->JacJacInv(geom_element.get_coords_at_dofs_3d(), i_qp, Jac_qp, JacI_qp, detJac_qp, space_dim);
    jacXweight_qp = detJac_qp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[i_qp];
    elem_all[ielGeom][solFEType_u]->shape_funcs_current_elem(i_qp, JacI_qp, phi_u, phi_u_x, boost::none /*phi_u_xx*/, space_dim);


//--------------    
	std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < nDof_u; i++) {
// 	                                                sol_u_gss      += sol_u[i] * phi_u[i];
                   for (unsigned d = 0; d < sol_u_x_gss.size(); d++)   sol_u_x_gss[d] += sol_u[i] * phi_u_x[i * space_dim + d];
          }
//--------------    

//--------------    
 /// @assignment You need to evaluate your manufactured right hand side at the quadrature point qp.
 /// Hence, you need to compute the coordinates of the quadrature point. Let us call them x_qp.
 /// These are obtained just like every quantity at a quadrature point, i.e., by interpolating the values of the quantity at the element nodes.
 /// The interpolation is performed by using the shape functions.
 /// In other words, 
 ///         (x_qp) = summation of (x_nodes) * (shape function of that node, evaluated at qp)
 /// 
 ///   (x_nodes) are obtained from   geom_element.get_coords_at_dofs_3d()  (this is a  vector< vector >,  where the outer index is the dimension and the inner index ranges over the nodes) 
 ///   (shape function of that node, evaluated at qp)  is obtained from phi_u  (this is a vector, whose index ranges over the nodes)
 
 std::vector<double> x_qp(dim, 0.);
          
        for (unsigned i = 0; i < nDof_u; i++) {
          	for (unsigned d = 0; d < dim; d++) {
	                                                x_qp[d]    += geom_element.get_coords_at_dofs_3d()[d][i] * phi_u[i]; // fetch of coordinate points
            }
        }
 
//--------------    

          
//==========FILLING WITH THE EQUATIONS ===========
	// *** phi_i loop ***
	double u_qp = 0.;
        for (unsigned i = 0; i < nDof_max; i++) {
            u_qp +=  sol_u[i] * phi_u[i];
        }

 	double u_true_qp = 0.;
        for (unsigned i = 0; i < nDof_max; i++) {
            u_true_qp +=  sol_u_true[i] * phi_u[i];
        }
       
        // //--------------    
// 	      double laplace_res_du_u_i = 0.;
//               if ( i < nDof_u ) {
//                   for (unsigned kdim = 0; kdim < space_dim; kdim++) {
//                        laplace_res_du_u_i             +=  phi_u_x   [i * space_dim + kdim] * sol_u_x_gss[kdim];
// 	             }
//               }
// //--------------    
              
	      
//======================Residuals=======================
          // FIRST ROW
 /// @assignment for your manufactured right-hand side, implement a function that receives the coordinate of the quadrature point
 /// Put it after the includes, in the top part of this file
/* if (i < nDof_u) */                     /*Res[0      + i]*/ norm +=  jacXweight_qp * (  u_qp - u_true_qp) * (  u_qp - u_true_qp) ;
//======================Residuals=======================
	      
// // //           if (assembleMatrix) {
// // // 	    
// // //             // *** phi_j loop ***
// // //             for (unsigned j = 0; j < nDof_max; j++) {
// // // 
// // // 
// // // //--------------    
// // //               double laplace_mat_du_u_i_j = 0.;
// // // 
// // //               if ( i < nDof_u && j < nDof_u ) {
// // //                   for (unsigned kdim = 0; kdim < space_dim; kdim++) {
// // //                          laplace_mat_du_u_i_j           += phi_u_x   [i * space_dim + kdim] *
// // //                                                            phi_u_x   [j * space_dim + kdim];
// // // 	              }
// // //              }
// // // //--------------    
// // // 
// // // 
// // // 
// // //               //============ delta_state row ============================
// // //               //DIAG BLOCK delta_state - state
// // // 		  if ( i < nDof_u && j < nDof_u )       Jac[ (0 + i) * nDof_AllVars   + 	(0 + j) ]  += jacXweight_qp * laplace_mat_du_u_i_j;
// // // // 		  if ( i < nDof_u && j < nDof_u )       Jac[ (0 + i) * nDof_AllVars   + 	(0 + j) ]  += jacXweight_qp * laplace_beltrami_mat_du_u_i_j; ///@todo On a flat domain, this must coincide with the standard Laplacian, so we can do a double check with this
// // //             } // end phi_j loop
// // //           } // endif assemble_matrix

        
      } // end gauss point loop


   
  } //end element loop for each process

  
  std::cout << std::scientific << std::setw(20) << std::setprecision(15) << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&Norm: " << norm << std::endl;
  



  return;
}


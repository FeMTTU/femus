/** tutorial/Ex3
 * This example shows how to set and solve the weak form of the nonlinear problem
 *                     -\Delta^2 u = f(x) \text{ on }\Omega,
 *            u=0 \text{ on } \Gamma,
 *      \Delta u=0 \text{ on } \Gamma,
 * on a box domain $\Omega$ with boundary $\Gamma$,
 * by using a system of second order partial differential equation.
 * all the coarse-level meshes are removed;
 * a multilevel problem and an equation system are initialized;
 * a direct solver is used to solve the problem.
 **/



#include "FemusInit.hpp"
#include "Files.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "LinearEquationSolver.hpp"
#include "VTKWriter.hpp"
#include "NumericVector.hpp"

#include "FE_convergence.hpp"

#include "Solution_functions_over_domains_or_mesh_files.hpp"

#include "adept.h"



#define LIBRARY_OR_USER   1 //0: library; 1: user

#if LIBRARY_OR_USER == 0
   #include "01_biharmonic_coupled.hpp"
   #define NAMESPACE_FOR_BIHARMONIC   femus
#elif LIBRARY_OR_USER == 1
   #include "biharmonic_coupled.hpp"
   #define NAMESPACE_FOR_BIHARMONIC_COUPLED   karthik
#endif



using namespace femus;


namespace Domains {

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



}

namespace  quarter_circle_centered_at_0_by_0  {
template < class type = double >
class Function_Zero_on_boundary_1_Laplacian : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {

        return  -12. *x[0] * x[1];
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        solGrad[0]  = -12. * x[1];
        solGrad[1]  = -12. * x[0];

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {

        return  0.;
    }


};

template < class type = double >
class Function_NonZero_on_boundary_1 : public Math::Function< type > {



public:

    type value(const std::vector < type >& x) const {

    // for a quarter-circle in Quadrant 1

    double xx = x[0];
    double yy = x[1];

    return  xx * yy * (1.0 + (xx*xx + yy*yy) ); // forced to be zero on the x and y axis, and the circle edge
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

    double xx = x[0];
    double yy = x[1];

        solGrad[0]  = yy * (1.0 + (xx*xx + yy*yy) ) + xx * yy * (2. * xx);
        solGrad[1]  = xx * (1.0 + (xx*xx + yy*yy) ) + yy * xx * (2. * yy);

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {

    double xx = x[0];
    double yy = x[1];

    return  12. * xx * yy;
    }



};


template < class type = double >
class Function_NonZero_on_boundary_1_Laplacian : public Math::Function< type > {

public:

    type value(const std::vector < type >& x) const {

        return  12. *x[0] * x[1];
    }


    std::vector < type >  gradient(const std::vector < type >& x) const {

        std::vector < type > solGrad(x.size(), 0.);

        solGrad[0]  = 12. * x[1];
        solGrad[1]  = 12. * x[0];

        return solGrad;
    }


    type laplacian(const std::vector < type >& x) const {

        return  0.;
    }


};



template <class type = double>
class Function_NonZero_on_boundary_2 : public Math::Function<type> {

public:
    type value(const std::vector<type>& x) const {
// // //         double xx = x[0];
// // //         double yy = x[1];

        return x[0] * x[1] * cos(x[0] * x[0] + x[1] * x[1]);
    }

    std::vector<type> gradient(const std::vector<type>& x) const {
        std::vector<type> solGrad(x.size(), 0.);

// // //         double xx = x[0];
// // //         double yy = x[1];

        solGrad[0] = x[1] * cos(x[0] * x[0] + x[1] * x[1]) - 2. * x[0] * x[0] * x[1] * sin(x[0] * x[0] + x[1] * x[1]);
        solGrad[1] = x[0] * cos(x[0] * x[0] + x[1] * x[1]) - 2. * x[0] * x[1] * x[1] * sin(x[0] * x[0] + x[1] * x[1]);

        return solGrad;
    }

    type laplacian(const std::vector<type>& x) const {

        return -8. * x[0] * x[1] * x[1] * x[1] * cos( x[0] * x[0]+ x[1] * x[1]) - 12. * x[0]*x[1]*sin(x[0] * x[0] + x[1] * x[1]);
    }
};

template <class type = double>
class Function_NonZero_on_boundary_2_Laplacian : public Math::Function<type> {

public:
    type value(const std::vector<type>& x) const {
        return -8. * x[0] * x[1] * x[1] * x[1] * cos( x[0] * x[0]+ x[1] * x[1]) - 12. * x[0]*x[1]*sin(x[0] * x[0] + x[1] * x[1]);
    }

    std::vector<type> gradient(const std::vector<type>& x) const {
        std::vector<type> solGrad(x.size(), 0.);

        solGrad[0] = 16. * x[0] * x[0] * x[1] * x[1] * x[1] * sin(x[0] * x[0] + x[1] * x[1]) - 8. * x[1] * x[1] * x[1] * cos(x[0] * x[0] + x[1] * x[1]) - 24. * x[0] * x[0] * x[1] * cos(x[0] * x[0] + x[1] * x[1])- 12. * x[1] * sin(x[0] * x[0] + x[1] * x[1]);
        solGrad[1] = 16. * x[0] * x[1] * x[1] * x[1] * x[1] * sin(x[0] * x[0] + x[1] * x[1]) - 48. * x[0] * x[1] * x[1] * cos(x[0] * x[0] + x[1] * x[1])- 12. * x[0] * sin(x[0] * x[0] + x[1] * x[1]);

        return solGrad;
    }

    type laplacian(const std::vector<type>& x) const {
        return 64. * x[0] * x[1] * x[1] * x[1] * x[1] * x[1] * sin(x[0] * x[0] + x[1] * x[1]) + 320. * x[0] * x[1] * x[1] * x[1] * sin(x[0] * x[0] + x[1] * x[1]) - 240. * x[0] * x[1] * cos(x[0] * x[0] + x[1] * x[1]);
    }
};




}




}



// namespace boundary_conditions{

//====Set boundary condition-BEGIN==============================
bool SetBoundaryCondition_bc_all_dirichlet_homogeneous(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& Value, const int facename, const double time) {
  bool dirichlet = true; //dirichlet

  if (!strcmp(SolName, "u")) {
      Math::Function <double> * u = ml_prob -> get_ml_solution() -> get_analytical_function(SolName);
      // strcmp compares two string in lexiographic sense.
    Value = u -> value(x);
  }
  else if (!strcmp(SolName, "v")) {
      Math::Function <double> * v = ml_prob -> get_ml_solution() -> get_analytical_function(SolName);
    Value = v -> value(x);
  }
  return dirichlet;
}
// }

//====Set boundary condition-END==============================


// // // //====Set boundary condition Dirichlet-Neumann-BEGIN==============================
// // //
// // // bool SetBoundaryCondition_bc_all_neumann_dirichlet(const MultiLevelProblem* ml_prob, const std::vector<double>& x, const char SolName[], double& value, const int facename, const double time) {
// // //     bool is_dirichlet = true;
// // //
// // //     if (!strcmp(SolName, "u")) {
// // //         // Assuming "u" corresponds to your first variable
// // //         Math::Function<double>* u = ml_prob->get_ml_solution()->get_analytical_function(SolName);
// // //            is_dirichlet = true;
// // //            value = u->value(x);
// // //     } else if (!strcmp(SolName, "v")) {
// // //         // Assuming "v" corresponds to your second variable
// // //            is_dirichlet = true;
// // //         Math::Function<double>* v = ml_prob->get_ml_solution()->get_analytical_function(SolName);
// // //         Math::Function<double>* u = ml_prob->get_ml_solution()->get_analytical_function("u");
// // //         // Set Dirichlet condition for "v" as the Laplacian of "u"
// // //         // value = u->laplacian(x);
// // //         value = v->value(x);
// // //     }
// // //
// // //     return is_dirichlet;
// // // }
// // // //====Set boundary condition Dirichlet-Neumann-END==============================






int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  MultiLevelProblem ml_prob;




  // ======= Files - BEGIN  ========================
  const bool use_output_time_folder = false; // This allows you to run the code multiple times without overwriting. This will generate an output folder each time you run.
  const bool redirect_cout_to_file = true; // puts the output in a log file instead of the term
  Files files;
        files.CheckIODirectories(use_output_time_folder);
        files.RedirectCout(redirect_cout_to_file);

  // ======= Files - END  ========================
  ml_prob.SetFilesHandler(& files);
  std::string fe_quad_rule("seventh");





// // //   std::string system_common_name1 = "Coupled_Biharmonic1";
// // //   std::string system_common_name2 = "Coupled_Biharmonic2";

  std::vector <system_specifics>  my_specifics;

  system_specifics app_square_m05p05_1;

  system_specifics app_square_m05p05_2;

  system_specifics quarter_circle;
  system_specifics quarter_circle_Nzero;



  const std::string relative_path_to_build_directory =  "../../../../";



  // ======= square 1 - BEGIN  ==================

  app_square_m05p05_1._system_name = "Coupled_Biharmonic1";

  app_square_m05p05_1._mesh_files.push_back("square_-0p5-0p5x-0p5-0p5_divisions_2x2.med");
  app_square_m05p05_1._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/square/minus0p5-plus0p5_minus0p5-plus0p5/");

   app_square_m05p05_1._assemble_function = NAMESPACE_FOR_BIHARMONIC_COUPLED :: biharmonic_coupled_equation::AssembleBilaplaceProblem_AD;

   app_square_m05p05_1._boundary_conditions_types_and_values             = SetBoundaryCondition_bc_all_dirichlet_homogeneous;


   Domains::square_m05p05::Function_Zero_on_boundary_4<>   app_square_function_zero_on_boundary_4_1;
   Domains::square_m05p05::Function_Zero_on_boundary_4_Laplacian<>   app_square_function_zero_on_boundary_4_1_laplacian;

   app_square_m05p05_1._assemble_function_for_rhs        = & app_square_function_zero_on_boundary_4_1_laplacian;
   app_square_m05p05_1._true_solution_function           = & app_square_function_zero_on_boundary_4_1;

    // ======= square 1 - END  ==================


    // ======= square 2 - BEGIN  ==================

   app_square_m05p05_2._system_name = "Coupled_Biharmonic2";

   app_square_m05p05_2._mesh_files.push_back("square_-0p5-0p5x-0p5-0p5_divisions_2x2.med");
   app_square_m05p05_2._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/square/minus0p5-plus0p5_minus0p5-plus0p5/");
   app_square_m05p05_2._assemble_function =NAMESPACE_FOR_BIHARMONIC_COUPLED:: biharmonic_coupled_equation::AssembleBilaplaceProblem_AD;

   app_square_m05p05_2._boundary_conditions_types_and_values         = SetBoundaryCondition_bc_all_dirichlet_homogeneous;

   Domains::square_m05p05::Function_NonZero_on_boundary_4<>   app_square_function_zero_on_boundary_4_2;
   Domains::square_m05p05::Function_NonZero_on_boundary_4_Laplacian<>   app_square_function_zero_on_boundary_4_2_laplacian;

   app_square_m05p05_2._assemble_function_for_rhs        = & app_square_function_zero_on_boundary_4_2_laplacian;
   app_square_m05p05_2._true_solution_function           = & app_square_function_zero_on_boundary_4_2;

    // ======= square 2 - END  ==================


    // ======= quarter_circle - BEGIN  ==================

   quarter_circle._system_name = "Coupled_Biharmonic3";

   quarter_circle._mesh_files.push_back("assignment_quarter_circle_quadrilateral.med");
   quarter_circle._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/circle_quarter/");
   quarter_circle._assemble_function =NAMESPACE_FOR_BIHARMONIC_COUPLED:: biharmonic_coupled_equation::AssembleBilaplaceProblem_AD;

   quarter_circle._boundary_conditions_types_and_values         = SetBoundaryCondition_bc_all_dirichlet_homogeneous;

   Domains::quarter_circle_centered_at_0_by_0::Function_Zero_on_boundary_1<>     app_quarter_circle_function_zero_on_boundary_1;
   Domains::quarter_circle_centered_at_0_by_0::Function_Zero_on_boundary_1_Laplacian<>   app_quarter_circle_function_zero_on_boundary_1_laplacian;

   quarter_circle._assemble_function_for_rhs        = & app_quarter_circle_function_zero_on_boundary_1_laplacian;
   quarter_circle._true_solution_function           = & app_quarter_circle_function_zero_on_boundary_1;

    // ======= quarter_circle - END  ==================


    // ======= quarter_circle_Nzero - BEGIN  ==================

   quarter_circle_Nzero._system_name = "Coupled_Biharmonic4";

   quarter_circle_Nzero._mesh_files.push_back("assignment_quarter_circle_quadrilateral.med");
   quarter_circle_Nzero._mesh_files_path_relative_to_executable.push_back(relative_path_to_build_directory + Files::mesh_folder_path() + "00_salome/2d/circle_quarter/");
   quarter_circle_Nzero._assemble_function =NAMESPACE_FOR_BIHARMONIC_COUPLED:: biharmonic_coupled_equation::AssembleBilaplaceProblem_AD;

   quarter_circle_Nzero._boundary_conditions_types_and_values         = SetBoundaryCondition_bc_all_dirichlet_homogeneous;

   Domains::quarter_circle_centered_at_0_by_0::Function_NonZero_on_boundary_2<>     app_quarter_circle_function_nonzero_on_boundary_1;
   Domains::quarter_circle_centered_at_0_by_0::Function_NonZero_on_boundary_2_Laplacian<>   app_quarter_circle_function_nonzero_on_boundary_1_laplacian;

   quarter_circle_Nzero._assemble_function_for_rhs        = & app_quarter_circle_function_nonzero_on_boundary_1_laplacian;
   quarter_circle_Nzero._true_solution_function           = & app_quarter_circle_function_nonzero_on_boundary_1;

    // ======= quarter_circle_Nzero - END  ==================

  my_specifics.push_back(app_square_m05p05_1);
  my_specifics.push_back(app_square_m05p05_2);
  my_specifics.push_back(quarter_circle);
  my_specifics.push_back(quarter_circle_Nzero);


 // ======= App loop - BEGIN  ==================

  for (unsigned int app = 0; app < my_specifics.size(); app++)  {


// //     ml_prob.set_app_specs_pointer(&my_specifics[app]);


    // ======= Mesh - BEGIN  ==================
  MultiLevelMesh ml_mesh;
  // ======= Mesh - END  ==================


   for (unsigned int m = 0; m < my_specifics[app]._mesh_files.size(); m++)  {

    // ======= Mesh, Coarse reading - BEGIN ==================
  double Lref = 1.;

  const bool read_groups = true; //with this being false, we don't read any group at all. Therefore, we cannot even read the boundary groups that specify what are the boundary faces, for the boundary conditions
  const bool read_boundary_groups = true;


  const std::string mesh_file = my_specifics[app]._mesh_files_path_relative_to_executable[m] + my_specifics[app]._mesh_files[m];

  ml_mesh.ReadCoarseMesh(mesh_file, Lref, read_groups, read_boundary_groups);
  // ======= Mesh, Coarse reading - END ==================


    double scalingFactor = 1.;
    ml_mesh.ReadCoarseMesh(mesh_file.c_str(), "seventh", scalingFactor);



  unsigned maxNumberOfMeshes = 5;

  std::vector < std::vector < double > > l2Norm;
  l2Norm.resize(maxNumberOfMeshes);

  std::vector < std::vector < double > > semiNorm;
  semiNorm.resize(maxNumberOfMeshes);

    std::vector<FEOrder> feOrder;
    feOrder.push_back(FIRST);
    feOrder.push_back(SERENDIPITY);
    feOrder.push_back(SECOND);



  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {   // loop on the mesh level

    unsigned numberOfUniformLevels = i + 1;
    unsigned numberOfSelectiveLevels = 0;
    ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);

    // erase all the coarse mesh levels
    ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);

    // print mesh info
    ml_mesh.PrintInfo();



    l2Norm[i].resize( feOrder.size() );
    semiNorm[i].resize( feOrder.size() );

    for (unsigned j = 0; j < feOrder.size(); j++) {   // loop on the FE Order

      // define the multilevel solution and attach the mlMsh object to it
      MultiLevelSolution mlSol(&ml_mesh);


      mlSol.SetWriter(VTK);
      mlSol.GetWriter()->SetDebugOutput(true);

      ml_prob.SetMultiLevelMeshAndSolution(& mlSol);



      mlSol.AddSolution("u", LAGRANGE, feOrder[j]);
      mlSol.set_analytical_function("u",  my_specifics[app]._true_solution_function);


      mlSol.AddSolution("v", LAGRANGE, feOrder[j]);
      mlSol.set_analytical_function("v",  my_specifics[app]._assemble_function_for_rhs);

      mlSol.Initialize("All");



      // define the multilevel problem attach the mlSol object to it
      MultiLevelProblem ml_prob(&mlSol);

      ml_prob.set_app_specs_pointer(& my_specifics[app]);

      mlSol.AttachSetBoundaryConditionFunction(my_specifics[app]._boundary_conditions_types_and_values);


      mlSol.GenerateBdc("u", "Steady", & ml_prob);
      mlSol.GenerateBdc("v", "Steady", & ml_prob);



      ml_prob.clear_systems();

      NonLinearImplicitSystem& system = ml_prob.add_system < NonLinearImplicitSystem > (my_specifics[app]._system_name);


      system.SetDebugNonlinear(true);


      // add solution "u" to system
      system.AddSolutionToSystemPDE("u");
      system.AddSolutionToSystemPDE("v");


      // attach the assembling function to system
      system.SetAssembleFunction(my_specifics[app]._assemble_function);

      // initialize and solve the system
      system.init();

      system.MGsolve();


      std::pair< double , double > norm = GetErrorNorm_L2_H1_with_analytical_sol(&mlSol, "u",my_specifics[app]._true_solution_function );

      l2Norm[i][j]  = norm.first;
      semiNorm[i][j] = norm.second;


      // print solutions
      const std::string print_order = fe_fams_for_files[ FILES_CONTINUOUS_BIQUADRATIC ];
      std::vector < std::string > variablesToBePrinted;
      variablesToBePrinted.push_back("All");

      VTKWriter vtkIO(&mlSol);
      vtkIO.Write(my_specifics[app]._system_name + "_" + my_specifics[app]._mesh_files[m], files.GetOutputPath(), print_order, variablesToBePrinted);


    }
  }


  // FE_convergence::output_convergence_order();


  // ======= L2 - BEGIN  ========================
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "l2 ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < feOrder.size(); j++) {
      std::cout << l2Norm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < feOrder.size(); j++) {
        std::cout << log(l2Norm[i][j] / l2Norm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }
  // ======= L2 - END  ========================



  // ======= H1 - BEGIN  ========================

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "SEMINORM ERROR and ORDER OF CONVERGENCE:\n\n";
  std::cout << "LEVEL\tFIRST\t\t\tSERENDIPITY\t\tSECOND\n";

  for (unsigned i = 0; i < maxNumberOfMeshes; i++) {
    std::cout << i + 1 << "\t";
    std::cout.precision(14);

    for (unsigned j = 0; j < feOrder.size(); j++) {
      std::cout << semiNorm[i][j] << "\t";
    }

    std::cout << std::endl;

    if (i < maxNumberOfMeshes - 1) {
      std::cout.precision(3);
      std::cout << "\t\t";

      for (unsigned j = 0; j < feOrder.size(); j++) {
        std::cout << log(semiNorm[i][j] / semiNorm[i + 1][j]) / log(2.) << "\t\t\t";
      }

      std::cout << std::endl;
    }

  }

  // ======= H1 - END  ========================

    }  //end mesh file loop


} //end app loop

 // ======= App loop - END  ==================


  return 0;
}







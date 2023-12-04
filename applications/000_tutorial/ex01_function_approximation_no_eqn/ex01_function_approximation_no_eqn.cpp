/** tutorial/Ex1
 * This example shows how to:
 * initialize a femus application;
 * define the multilevel-mesh object mlMsh;
 * read from the file ./input/square.neu the coarse-level mesh and associate it to mlMsh;
 * add in mlMsh uniform refined level-meshes;
 * define the multilevel-solution object mlSol associated to mlMsh;
 * add in mlSol different types of finite element solution variables;
 * initialize the solution variables;
 * define vtk and gmv writer objects associated to mlSol;
 * print vtk and gmv binary-format files in ./output directory.
 **/

#include "FemusInit.hpp"
#include "MultiLevelMesh.hpp"
#include "MultiLevelSolution.hpp"
#include "Files.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"

#include "Solution_functions_over_domains_or_mesh_files.hpp"


#include <tuple>
#include <vector>
#include <string>



using namespace femus;

double InitialValueU(const std::vector < double >& x) {
  return x[0] + x[1];
}

double InitialValueP(const std::vector < double >& x) {
  return x[0];
}

double InitialValueT(const std::vector < double >& x) {
  return x[1];
}


// Mesh-dependent functions - BEGIN ===============

bool SetRefinementFlag(const std::vector < double >& x, const int &elemgroupnumber,const int &level) {

  bool refine = false;
  if ( elemgroupnumber == 6 && level < 1 ) refine = true;
  if ( elemgroupnumber == 7 && level < 2 ) refine = true;
  if ( elemgroupnumber == 8 && level < 3 ) refine = true;

  return refine;

}

// Mesh-dependent functions - END ===============

// // // bool SetBoundaryCondition(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
// // //   
// // //   bool dirichlet = true; //dirichlet
// // //   value = 0.;
// // // 
// // //   return dirichlet;
// // // }



int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  

  // Mesh - BEGIN
  
  // define multilevel mesh
  MultiLevelMesh mlMsh;
  double scalingFactor = 1.;
  
  // read coarse level mesh and generate finer level meshes 

  typedef std::pair< std::string, std::string > Mesh_file_name_info;
  

// === BEGIN  
  const std::string relative_path_to_build_directory =  "../../../";
  
  
  const std::string mesh_file_1_path = relative_path_to_build_directory + Files::mesh_folder_path() + 
  "01_gambit/02_2d/square/minus0p5-plus0p5_minus0p5-plus0p5/";

  const std::string mesh_file_1_name = "square_16x16_quad_One_boundary_group.neu";
  
  
  const std::string mesh_file_2_path = relative_path_to_build_directory + Files::mesh_folder_path() +
  "01_gambit/02_2d/square/minus0p5-plus0p5_minus0p5-plus0p5/";
  
  const std::string mesh_file_2_name = "square_4x4_quad_Four_boundary_groups_Four_volume_groups_AMR.neu";

  const std::string mesh_file_3_path = relative_path_to_build_directory + Files::mesh_folder_path() + 
  "01_gambit/03_3d/cube/minus0p5-plus0p5_minus0p5-plus0p5_minus0p5-plus0p5/";

  const std::string mesh_file_3_name = "cube_mixed_Two_boundary_groups_Four_volume_groups_AMR.neu";

  
  std::vector< Mesh_file_name_info >   meshes_path_and_filename;
  
  meshes_path_and_filename.push_back(  Mesh_file_name_info ( mesh_file_1_path, mesh_file_1_name ) );
  meshes_path_and_filename.push_back(  Mesh_file_name_info ( mesh_file_2_path, mesh_file_2_name ) );
  meshes_path_and_filename.push_back(  Mesh_file_name_info ( mesh_file_3_path, mesh_file_3_name ) );
// === END
  
  
// === BEGIN  
  typedef std::tuple < Mesh_file_name_info, MultiLevelMesh::RefinementFunctionBasedOnVolumeGroups, unsigned, unsigned  >   Mesh_file_with_its_refinement_information;

  constexpr unsigned index_for_mesh_file_info = 0;
  constexpr unsigned index_for_function_pointer = 1;
  constexpr unsigned index_for_first_refine_arg = 2;
  constexpr unsigned index_for_second_refine_arg = 3;

  
  std::vector< Mesh_file_with_its_refinement_information >  meshes_and_refinements;
  
  
  meshes_and_refinements.push_back( Mesh_file_with_its_refinement_information (meshes_path_and_filename[0], NULL,              3, 3) );
  meshes_and_refinements.push_back( Mesh_file_with_its_refinement_information (meshes_path_and_filename[1], SetRefinementFlag, 4, 1) );
  meshes_and_refinements.push_back( Mesh_file_with_its_refinement_information (meshes_path_and_filename[2], SetRefinementFlag, 4, 1) );
// === END


    for (unsigned mesh_file_index = 0; mesh_file_index < meshes_and_refinements.size(); mesh_file_index ++) {
  
    
  const std::string  mesh_full_filename =
  std::get< index_for_mesh_file_info >( meshes_and_refinements[ mesh_file_index ] ).first + 
  std::get< index_for_mesh_file_info >( meshes_and_refinements[ mesh_file_index ] ).second;
  
  mlMsh.ReadCoarseMesh( mesh_full_filename.c_str(), "seventh", scalingFactor);
  /* "seventh" is the order of accuracy that is used in the gauss integration scheme
      probably in the future it is not going to be an argument of this function   */

  MultiLevelMesh::RefinementFunctionBasedOnVolumeGroups  refinement_function_pointer = std::get< index_for_function_pointer >( meshes_and_refinements[ mesh_file_index ] );


// non amr
  // unsigned numberOfUniformLevels = 3;
  // unsigned numberOfSelectiveLevels = 0;
// amr
  // unsigned numberOfUniformLevels = 1;
  // unsigned numberOfSelectiveLevels = 3;
  
  mlMsh.RefineMesh( std::get< index_for_first_refine_arg > ( meshes_and_refinements[ mesh_file_index ] ), 
                    std::get< index_for_second_refine_arg >( meshes_and_refinements[ mesh_file_index ] ), 
                    refinement_function_pointer);
  
  mlMsh.PrintInfo();
  
  // Mesh - END

  
  // Solution - BEGIN
  
  // define the multilevel solution and attach the mlMsh object to it
  MultiLevelSolution mlSol(&mlMsh);
  
  // add variables to mlSol
  mlSol.AddSolution("U", LAGRANGE, FIRST);
  mlSol.AddSolution("V", LAGRANGE, SERENDIPITY);
  mlSol.AddSolution("W", LAGRANGE, SECOND);
  mlSol.AddSolution("P", DISCONTINUOUS_POLYNOMIAL, ZERO);
  mlSol.AddSolution("T", DISCONTINUOUS_POLYNOMIAL, FIRST);


  Domains::Function_X<>  function_x;
  Domains::Function_Y<>  function_y;
  Domains::Function_XplusY<>   function_x_plus_y;
  
  mlSol.set_analytical_function("U", & function_x_plus_y);
  mlSol.set_analytical_function("V", & function_x_plus_y);   
  mlSol.set_analytical_function("W", & function_x_plus_y);   
  mlSol.set_analytical_function("P", & function_x);   
  mlSol.set_analytical_function("T", & function_y);   
  
  mlSol.Initialize("All");    // initialize all variables to zero

  mlSol.Initialize("U", InitialValueU);
  mlSol.Initialize("V", InitialValueU);
  mlSol.Initialize("W", InitialValueU);
  mlSol.Initialize("P", InitialValueP);
  mlSol.Initialize("T", InitialValueT);
  
  // // // mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  // // // mlSol.GenerateBdc("All");
  // Solution - END


  // Solution, print - BEGIN
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.SetDebugOutput(true);
  
  // Strip extension from filename (to avoid a Paraview message) - BEGIN
  const size_t lastindex    = (std::get< index_for_mesh_file_info >( meshes_and_refinements[ mesh_file_index ] ).second).find_last_of("."); 
  const std::string name_without_extension = (std::get< index_for_mesh_file_info >( meshes_and_refinements[ mesh_file_index ] ).second).substr(0, lastindex); 
  // Strip extension from filename (to avoid a Paraview message) - END
  
  vtkIO.Write( name_without_extension,
               Files::_application_output_directory, 
               fe_fams_for_files[ FILES_CONTINUOUS_BIQUADRATIC ], 
               variablesToBePrinted);
  vtkIO.Write( name_without_extension,
               Files::_application_output_directory, 
               fe_fams_for_files[ FILES_CONTINUOUS_LINEAR ], 
               variablesToBePrinted);

  // // // GMVWriter gmvIO(&mlSol);
  // // // gmvIO.SetDebugOutput(true);
  // // // gmvIO.Write(std::get< index_for_mesh_file_info >( meshes_and_refinements[ mesh_file_index ] ).second,
  // // //             Files::_application_output_directory, 
  // // //             fe_fams_for_files[ FILES_CONTINUOUS_BIQUADRATIC ],
  // // //             variablesToBePrinted);
  // Solution, print - END

  
    }
    
    
  
  return 0;
}

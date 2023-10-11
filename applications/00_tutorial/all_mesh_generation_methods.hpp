#ifndef __femus_applications_all_mesh_generation_methods_hpp__
#define __femus_applications_all_mesh_generation_methods_hpp__



#include "MultiLevelMesh.hpp"


namespace  Domain_square_m05p05  {


std::string quad9_all_mesh_generation_methods_Structured(const unsigned method_flag, MultiLevelMesh & ml_mesh, const std::string fe_quad_rule) {
    
    
  const std::string relative_path_to_build_directory =  "../../../"; 

    
  double scalingFactor = 1.;

    std::string mesh_name("square");
    
    switch(method_flag) {
        case 0: {
            mesh_name += "_femus";
  // ======= Mesh, coarse gen, I: from function - BEGIN  ========================
  const unsigned int nsub_x = 2;
  const unsigned int nsub_y = 2;
  const unsigned int nsub_z = 0;
  const std::vector<double> xyz_min = {-0.5,-0.5,0.};
  const std::vector<double> xyz_max = { 0.5, 0.5,0.};
  const ElemType geom_elem_type = QUAD9/*TRI6*/;
  ml_mesh.GenerateCoarseBoxMesh(nsub_x, nsub_y, nsub_z, xyz_min[0], xyz_max[0], xyz_min[1], xyz_max[1], xyz_min[2], xyz_max[2], geom_elem_type, fe_quad_rule.c_str() );
  // ======= Mesh, coarse gen, I: from function - END ========================
 
            break;
        }
        
        case 1: {
            mesh_name += "_salome";
//   // ======= Mesh, coarse gen, III: from Salome - BEGIN  ========================
    const std::string mesh_file = relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "00_salome/02_2d/square_-0p5-0p5x-0p5-0p5/square_2x2_centered_at_origin.med";
   ml_mesh.ReadCoarseMesh(mesh_file.c_str(), fe_quad_rule.c_str(), scalingFactor);
//   // ======= Mesh, coarse gen, III: from Salome - END ========================
            break;
        }
        case 2: {
            mesh_name += "_gambit";
//   // ======= Mesh, coarse gen, II: from Gambit - BEGIN  ========================
    const std::string mesh_file = relative_path_to_build_directory + DEFAULT_MESH_FILES_PATH + "01_gambit/square_-0p5-0p5x-0p5-0p5/square_2x2_quad.neu";
    ml_mesh.ReadCoarseMesh(mesh_file.c_str(), fe_quad_rule.c_str(), scalingFactor);
  // ======= Mesh, coarse gen, II: from Gambit - END ========================
            break;
        }
        
        default: { abort(); }
    }
    
    
    return mesh_name;
    
}

    
}


#endif

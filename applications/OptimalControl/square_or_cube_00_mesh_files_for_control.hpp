#ifndef SQUARE_OR_CUBE_MESH_FILES_FOR_CONTROL_HPP
#define SQUARE_OR_CUBE_MESH_FILES_FOR_CONTROL_HPP
 
//******************************** Mesh files oriented to Optimal Control - BEGIN *****************************************

namespace femus {


namespace ctrl {

namespace mesh {

  const std::string mesh_2d_square_1x1 = "parametric_square_1x1.med";
  const std::string mesh_2d_square_1x2 = "parametric_square_1x2.med";
  const std::string mesh_2d_square_2x2 = "parametric_square_2x2.med";
  const std::string mesh_2d_square_4x5 = "parametric_square_4x5.med";

  const std::string mesh_3d_cube_single_face_control_1             = "GroupofNode_Face1.med";
  const std::string mesh_3d_cube_single_face_control_2             = "GroupofNode_Face2.med";
  const std::string mesh_3d_cube_single_face_control_2_old         = "Mesh_3_groups_with_bdry_nodes.med";
  const std::string mesh_3d_cube_single_face_control_2_old_coarser = "Mesh_3_groups_with_bdry_nodes_coarser.med";
  const std::string mesh_3d_cube_single_face_control_3             = "GroupofNode_Face3.med";
  const std::string mesh_3d_cube_single_face_control_4             = "GroupofNode_Face4.med";

  
  
   
   
    
  std::string file_with_prefix(const std::string input_file)  {
     
     
     const std::string mesh_location = "../../../";
     
      std::ostringstream mystream; mystream << mesh_location /*"./"*/ << DEFAULT_INPUTDIR << "/" << input_file;
      const std::string infile = mystream.str();
      
      return infile;

   }
     
     
}
      
}

}
//******************************** Mesh files oriented to Optimal Control - END *****************************************


#endif

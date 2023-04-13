#ifndef SQUARE_OR_CUBE_MESH_FILES_FOR_CONTROL_HPP
#define SQUARE_OR_CUBE_MESH_FILES_FOR_CONTROL_HPP
 
//******************************** Mesh files oriented to Optimal Control - BEGIN *****************************************

namespace femus {


namespace ctrl {

namespace square {

namespace mesh {

  const std::string _2d_square_1x1 = "parametric_square_1x1.med";
  const std::string _2d_square_1x2 = "parametric_square_1x2.med";
  const std::string _2d_square_2x2 = "parametric_square_2x2.med";
  const std::string _2d_square_4x5 = "parametric_square_4x5.med";
}

}


namespace cube {

namespace mesh {
    
  const std::string _3d_cube_single_face_control_1             = "GroupofNode_Face1.med";
  const std::string _3d_cube_single_face_control_2             = "GroupofNode_Face2.med";
  const std::string _3d_cube_single_face_control_2_old         = "Mesh_3_groups_with_bdry_nodes.med";
  const std::string _3d_cube_single_face_control_2_old_coarser = "Mesh_3_groups_with_bdry_nodes_coarser.med";
  const std::string _3d_cube_single_face_control_2_old_coarser_between = "Mesh_3_groups_with_bdry_nodes_coarser_BETWEEN_BETWEEN.med";
  const std::string _3d_cube_single_face_control_3             = "GroupofNode_Face3.med";
  const std::string _3d_cube_single_face_control_4             = "GroupofNode_Face4.med";
  const std::string _3d_cube_all_faces_between_between         = "All_Face_For_Control_Parametric_Group_of_Nodes_MERGED_FACES_BETWEEN_BETWEEN.med";


}
      
}

}

}
//******************************** Mesh files oriented to Optimal Control - END *****************************************


#endif

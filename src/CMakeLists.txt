CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

PROJECT(femus)

#############################################################################################
### Source files - BEGIN
#############################################################################################

SET(femus_src 
00_utils/00_application_initialization/FemusInit.cpp
00_file_handling/application_file_handling/Files.cpp
00_file_handling/runtime_input_parsing/file/InputParser.cpp
00_file_handling/runtime_input_parsing/file/JsonInputParser.cpp
01_materials/Material.cpp
01_materials/Fluid.cpp
01_materials/Solid.cpp
01_materials/Parameter.cpp
02_calculus/function_parser/FunctionBase.cpp
02_calculus/function_parser/ParsedFunction.cpp
02_reference_geom_elements/00_definition/GeomElemBase.cpp
02_reference_geom_elements/00_definition/1d/GeomElemEdge2.cpp
02_reference_geom_elements/00_definition/1d/GeomElemEdge3.cpp
02_reference_geom_elements/00_definition/2d/GeomElemTri3.cpp
02_reference_geom_elements/00_definition/2d/GeomElemTri6.cpp
02_reference_geom_elements/00_definition/2d/GeomElemTri7.cpp
02_reference_geom_elements/00_definition/2d/GeomElemQuad4.cpp
02_reference_geom_elements/00_definition/2d/GeomElemQuad9.cpp
02_reference_geom_elements/00_definition/3d/GeomElemHex8.cpp
02_reference_geom_elements/00_definition/3d/GeomElemHex27.cpp
02_reference_geom_elements/00_definition/3d/GeomElemTet4.cpp
02_reference_geom_elements/00_definition/3d/GeomElemTet10.cpp
02_reference_geom_elements/01_fe/1d/Edge.cpp
02_reference_geom_elements/01_fe/2d/Quadrilateral.cpp
02_reference_geom_elements/01_fe/2d/Triangle.cpp
02_reference_geom_elements/01_fe/3d/Hexahedron.cpp
02_reference_geom_elements/01_fe/3d/Tetrahedron.cpp
02_reference_geom_elements/01_fe/3d/Wedge.cpp
02_reference_geom_elements/02_quadrature/quadrature_interface.cpp
02_reference_geom_elements/02_quadrature/0d/quadrature_Point.cpp
02_reference_geom_elements/02_quadrature/1d/quadrature_Line.cpp
02_reference_geom_elements/02_quadrature/2d/quadrature_Quadrangle.cpp
02_reference_geom_elements/02_quadrature/2d/quadrature_Triangle.cpp
02_reference_geom_elements/02_quadrature/3d/quadrature_Hexahedron.cpp
02_reference_geom_elements/02_quadrature/3d/quadrature_Tetrahedron.cpp
02_reference_geom_elements/02_quadrature/3d/quadrature_Wedge.cpp
02_reference_geom_elements/03_fe_evaluations_at_quadrature/ElemType_template.cpp
02_reference_geom_elements/03_fe_evaluations_at_quadrature/ElemType_template_base.cpp
02_reference_geom_elements/03_fe_evaluations_at_quadrature/ElemType.cpp
03_algebra_dense/00_vectors/DenseVectorBase.cpp
03_algebra_dense/01_matrices/DenseMatrixBase.cpp
03_algebra_dense/01_matrices/DenseMatrix.cpp
03_algebra_dense/01_matrices/DenseSubmatrix.cpp
03_algebra_only_mpi/00_vectors/MyVector.cpp
03_algebra_only_mpi/01_matrices/MyMatrix.cpp
03_algebra/00_vectors/NumericVector.cpp
03_algebra/00_vectors/PetscVector.cpp
03_algebra/01_matrices/PetscMatrix.cpp
03_algebra/01_matrices/SparseMatrix.cpp
03_algebra/01_matrices/Graph.cpp
03_algebra/01_matrices/SlepcSVD.cpp
03_algebra/02_preconditioners/PetscPreconditioner.cpp
03_algebra/02_preconditioners/Preconditioner.cpp
05_problem/MultiLevelProblem.cpp
06_mesh/00_single_level/00_definition/Mesh.cpp
06_mesh/00_single_level/00_definition/Elem.cpp
06_mesh/00_single_level/00_definition/fe_projection_matrices_Lagrange_continuous.cpp
06_mesh/00_single_level/00_definition/fe_prolongation_matrices.cpp
06_mesh/00_single_level/01_input/01_from_external_file/GambitIO.cpp
06_mesh/00_single_level/01_input/01_from_external_file/MED_IO.cpp
06_mesh/00_single_level/01_input/02_from_implemented_code/MeshGeneration.cpp
06_mesh/00_single_level/02_partitioning/MeshPartitioning.cpp
06_mesh/00_single_level/02_partitioning/MeshMetisPartitioning.cpp
06_mesh/00_single_level/02_partitioning/MeshASMPartitioning.cpp
06_mesh/00_single_level/03_refinement/MeshRefinement.cpp
06_mesh/01_multiple_levels/00_definition/MultiLevelMesh.cpp
06_mesh/gencase_deprecated/Box.cpp
06_mesh/gencase_deprecated/Domain.cpp
06_mesh/gencase_deprecated/ElemSto.cpp
06_mesh/gencase_deprecated/MultiLevelMeshTwo.cpp
06_mesh/gencase_deprecated/GenCase.cpp
06_solution/00_single_level/00_definition/Solution.cpp
06_solution/01_multiple_levels/00_definition/MultiLevelSolution.cpp
07_mesh_or_solution/00_single_level/01_output/Writer_one_level.cpp
07_mesh_or_solution/01_multiple_levels/01_output/Writer.cpp
07_mesh_or_solution/01_multiple_levels/01_output/VTKWriter.cpp
07_mesh_or_solution/01_multiple_levels/01_output/GMVWriter.cpp
07_mesh_or_solution/01_multiple_levels/01_output/XDMFWriter.cpp
08_algebra_dependent_on_Mesh_and_Solution_but_independent_of_Systems/03_solvers_with_preconditioner/LinearEquation.cpp
08_algebra_dependent_on_Mesh_and_Solution_but_independent_of_Systems/03_solvers_with_preconditioner/LinearEquationSolver.cpp
08_algebra_dependent_on_Mesh_and_Solution_but_independent_of_Systems/03_solvers_with_preconditioner/LinearEquationSolverPetsc.cpp
08_algebra_dependent_on_Mesh_and_Solution_but_independent_of_Systems/03_solvers_with_preconditioner/petsc_asm/LinearEquationSolverPetscAsm.cpp
08_algebra_dependent_on_Mesh_and_Solution_but_independent_of_Systems/03_solvers_with_preconditioner/petsc_field_split/FieldSplitTree.cpp
08_algebra_dependent_on_Mesh_and_Solution_but_independent_of_Systems/03_solvers_with_preconditioner/petsc_field_split/LinearEquationSolverPetscFieldSplit.cpp
08_equations/00_stationary/System.cpp
08_equations/00_stationary/LinearImplicitSystem.cpp
08_equations/00_stationary/NonLinearImplicitSystem.cpp
08_equations/01_time_dependent/TransientSystem.cpp
08_equations/assemble/Quantity.cpp
08_equations/assemble/DofMap.cpp
08_equations/assemble/Assemble_jacobian.cpp
08_equations/assemble/Assemble_unknown.cpp
08_equations/assemble/BoundaryConditions.cpp
08_equations/assemble/CurrentElem.cpp
08_equations/assemble/CurrentQuantity.cpp
08_equations/assemble/SystemTwo.cpp
08_equations/assemble/TimeLoop.cpp
08_fsi/00_stationary/MonolithicFSINonLinearImplicitSystem.cpp
08_fsi/01_time_dependent/TransientFSI.cpp
09_optimal_control/00_cost_functional.cpp
09_optimal_control/NonLinearImplicitSystemWithPrimalDualActiveSetMethod.cpp
ism/Marker.cpp
ism/PolynomialBases.cpp
ism/Line.cpp
uq/uq.cpp
uq/sparseGrid.cpp
)

# IF(HAVE_HDF5)
#   SET(femus_src ${femus_src} 
#   06_mesh/input/MED_IO.cpp
#   07_solution/01_output/XDMFWriter.cpp
#   06_mesh/gencase/GenCase.cpp
#   08_equations/TimeLoop.cpp
#   )
# ENDIF(HAVE_HDF5) 



#############################################################################################
### Source files - END
#############################################################################################



#############################################################################################
### Library files - BEGIN
#############################################################################################


IF (NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib64 CACHE PATH
       "Single output directory for building all libraries.")
  MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
ENDIF (NOT LIBRARY_OUTPUT_PATH)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${femus_src})



IF(APPLE)

target_link_libraries(${PROJECT_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libhdf5.dylib)
target_link_libraries(${PROJECT_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libmpi.dylib)
target_link_libraries(${PROJECT_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.dylib)
target_link_libraries(${PROJECT_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libmetis.dylib)
target_link_libraries(${PROJECT_NAME} ${SLEPC_DIR}/${PETSC_ARCH}/lib/libslepc.dylib)
#target_link_libraries(${PROJECT_NAME} /Users/anthonygruber/4Femus/slepc/${PETSC_ARCH}/lib/libslepc.dylib)
# target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/libjsoncpp.dylib)
# target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/libadept.a)
# target_link_libraries(${PROJECT_NAME} ${LIBRARY_OUTPUT_PATH}/libb64.a)

#target_link_libraries(${PROJECT_NAME} ${PETSC_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${B64_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${JSONCPP_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${ADEPT_LIBRARIES})

ENDIF(APPLE)

#############################################################################################
### Library files - END
#############################################################################################

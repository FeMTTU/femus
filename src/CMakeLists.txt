CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(femus)

#############################################################################################
### Source files - BEGIN
#############################################################################################

SET(femus_src 
00_utils/00_application_initialization/FemusInit.cpp
00_utils/01_application_file_handling/Files.cpp
00_utils/input_parser/InputParser.cpp
00_utils/input_parser/JsonInputParser.cpp
00_utils/function_parser/FunctionBase.cpp
00_utils/function_parser/ParsedFunction.cpp
00_utils/parallel/MyMatrix.cpp
00_utils/parallel/MyVector.cpp
01_problem/MultiLevelProblem.cpp
02_geom_elements/00_definition/GeomElemBase.cpp
02_geom_elements/00_definition/GeomElemHex8.cpp
02_geom_elements/00_definition/GeomElemHex27.cpp
02_geom_elements/00_definition/GeomElemQuad4.cpp
02_geom_elements/00_definition/GeomElemQuad9.cpp
02_geom_elements/00_definition/GeomElemTet4.cpp
02_geom_elements/00_definition/GeomElemTet10.cpp
02_geom_elements/00_definition/GeomElemTri3.cpp
02_geom_elements/00_definition/GeomElemTri6.cpp
02_geom_elements/00_definition/GeomElemTri7.cpp
02_geom_elements/00_definition/GeomElemEdge2.cpp
02_geom_elements/00_definition/GeomElemEdge3.cpp
02_geom_elements/01_fe/ElemType.cpp
02_geom_elements/01_fe/Hexaedron.cpp
02_geom_elements/01_fe/Line.cpp
02_geom_elements/01_fe/Quadrilateral.cpp
02_geom_elements/01_fe/Tetrahedron.cpp
02_geom_elements/01_fe/Triangle.cpp
02_geom_elements/01_fe/Wedge.cpp
02_geom_elements/02_quadrature/GaussPoints.cpp
02_geom_elements/03_fe_evaluations_at_quadrature/ElemType_template.cpp
02_geom_elements/03_fe_evaluations_at_quadrature/ElemType_template_base.cpp
02_mesh/00_single_level/Elem.cpp
02_mesh/00_single_level/Mesh.cpp
02_mesh/00_single_level/01_input/01_from_external_file/GambitIO.cpp
02_mesh/00_single_level/01_input/01_from_external_file/MED_IO.cpp
02_mesh/00_single_level/01_input/02_from_implemented_code/MeshGeneration.cpp
02_mesh/00_single_level/02_partitioning/MeshPartitioning.cpp
02_mesh/00_single_level/02_partitioning/MeshMetisPartitioning.cpp
02_mesh/00_single_level/02_partitioning/MeshASMPartitioning.cpp
02_mesh/00_single_level/03_refinement/MeshRefinement.cpp
02_mesh/01_multiple_levels/MultiLevelMesh.cpp
02_mesh/gencase/Box.cpp
02_mesh/gencase/Domain.cpp
02_mesh/gencase/ElemSto.cpp
02_mesh/gencase/MultiLevelMeshTwo.cpp
02_mesh/gencase/GenCase.cpp
03_solution/00_single_level/Solution.cpp
03_solution/01_multiple_levels/MultiLevelSolution.cpp
03_solution/01_multiple_levels/01_output/Writer.cpp
03_solution/01_multiple_levels/01_output/VTKWriter.cpp
03_solution/01_multiple_levels/01_output/GMVWriter.cpp
03_solution/01_multiple_levels/01_output/XDMFWriter.cpp
04_equations/00_stationary/LinearImplicitSystem.cpp
04_equations/00_stationary/MonolithicFSINonLinearImplicitSystem.cpp
04_equations/00_stationary/NonLinearImplicitSystem.cpp
04_equations/00_stationary/System.cpp
04_equations/01_time_dependent/TransientSystem.cpp
04_equations/assemble/Quantity.cpp
04_equations/assemble/DofMap.cpp
04_equations/assemble/Assemble_jacobian.cpp
04_equations/assemble/Assemble_unknown.cpp
04_equations/assemble/BoundaryConditions.cpp
04_equations/assemble/CurrentElem.cpp
04_equations/assemble/CurrentQuantity.cpp
04_equations/assemble/SystemTwo.cpp
04_equations/assemble/TimeLoop.cpp
algebra/00_vectors/DenseVectorBase.cpp
algebra/00_vectors/NumericVector.cpp
algebra/00_vectors/PetscVector.cpp
algebra/01_matrices/DenseMatrixBase.cpp
algebra/01_matrices/DenseMatrix.cpp
algebra/01_matrices/DenseSubmatrix.cpp
algebra/01_matrices/PetscMatrix.cpp
algebra/01_matrices/SparseMatrix.cpp
algebra/01_matrices/Graph.cpp
algebra/02_preconditioners/PetscPreconditioner.cpp
algebra/02_preconditioners/Preconditioner.cpp
algebra/03_solvers_with_preconditioner/LinearEquation.cpp
algebra/03_solvers_with_preconditioner/LinearEquationSolver.cpp
algebra/03_solvers_with_preconditioner/LinearEquationSolverPetsc.cpp
algebra/03_solvers_with_preconditioner/petsc_asm/LinearEquationSolverPetscAsm.cpp
algebra/03_solvers_with_preconditioner/petsc_field_split/FieldSplitTree.cpp
algebra/03_solvers_with_preconditioner/petsc_field_split/LinearEquationSolverPetscFieldSplit.cpp
algebra/SlepcSVD.cpp
ism/Marker.cpp
ism/PolynomialBases.cpp
ism/Line.cpp
optimal_control/00_cost_functional.cpp
optimal_control/NonLinearImplicitSystemWithPrimalDualActiveSetMethod.cpp
physics/Fluid.cpp
physics/Material.cpp
physics/Parameter.cpp
physics/Solid.cpp
uq/uq.cpp
uq/sparseGrid.cpp
)

# IF(HAVE_HDF5)
#   SET(femus_src ${femus_src} 
#   02_mesh/input/MED_IO.cpp
#   03_solution/01_output/XDMFWriter.cpp
#   02_mesh/gencase/GenCase.cpp
#   04_equations/TimeLoop.cpp
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

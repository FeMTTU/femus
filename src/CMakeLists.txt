CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(femus)

#############################################################################################
### Source files
#############################################################################################

SET(femus_src 
00_utils/FemusInit.cpp
00_utils/Files.cpp
00_utils/input_parser/InputParser.cpp
00_utils/input_parser/JsonInputParser.cpp
01_mesh/gencase/Box.cpp
01_mesh/gencase/Domain.cpp
01_mesh/gencase/ElemSto.cpp
01_mesh/gencase/MultiLevelMeshTwo.cpp
01_mesh/gencase/GenCase.cpp
01_mesh/00_geom_elements/GeomElemBase.cpp
01_mesh/00_geom_elements/GeomElemHex8.cpp
01_mesh/00_geom_elements/GeomElemHex27.cpp
01_mesh/00_geom_elements/GeomElemQuad4.cpp
01_mesh/00_geom_elements/GeomElemQuad9.cpp
01_mesh/00_geom_elements/GeomElemTet4.cpp
01_mesh/00_geom_elements/GeomElemTet10.cpp
01_mesh/00_geom_elements/GeomElemTri3.cpp
01_mesh/00_geom_elements/GeomElemTri6.cpp
01_mesh/00_geom_elements/GeomElemTri7.cpp
01_mesh/00_geom_elements/GeomElemEdge2.cpp
01_mesh/00_geom_elements/GeomElemEdge3.cpp
01_mesh/01_input/MeshGeneration.cpp
01_mesh/01_input/GambitIO.cpp
01_mesh/01_input/MED_IO.cpp
01_mesh/02_partitioning/MeshPartitioning.cpp
01_mesh/02_partitioning/MeshMetisPartitioning.cpp
01_mesh/02_partitioning/MeshASMPartitioning.cpp
01_mesh/03_refinement/MeshRefinement.cpp
01_mesh/Elem.cpp
01_mesh/Mesh.cpp
01_mesh/MultiLevelMesh.cpp
02_solution/MultiLevelSolution.cpp
02_solution/Solution.cpp
02_solution/01_output/Writer.cpp
02_solution/01_output/VTKWriter.cpp
02_solution/01_output/GMVWriter.cpp
02_solution/01_output/XDMFWriter.cpp
03_equations/assemble/Quantity.cpp
03_equations/assemble/DofMap.cpp
03_equations/assemble/Assemble_jacobian.cpp
03_equations/assemble/Assemble_unknown.cpp
03_equations/assemble/BoundaryConditions.cpp
03_equations/assemble/CurrentElem.cpp
03_equations/assemble/CurrentQuantity.cpp
03_equations/assemble/SystemTwo.cpp
03_equations/ExplicitSystem.cpp
03_equations/ImplicitSystem.cpp
03_equations/LinearImplicitSystem.cpp
03_equations/MonolithicFSINonLinearImplicitSystem.cpp
03_equations/NonLinearImplicitSystem.cpp
03_equations/NonLinearImplicitSystemWithPrimalDualActiveSetMethod.cpp
03_equations/MultiLevelProblem.cpp
03_equations/System.cpp
03_equations/TimeLoop.cpp
03_equations/TransientSystem.cpp
algebra/DenseMatrixBase.cpp
algebra/DenseMatrix.cpp
algebra/DenseSubmatrix.cpp
algebra/DenseVectorBase.cpp
algebra/FieldSplitTree.cpp
algebra/Graph.cpp
algebra/LinearEquation.cpp
algebra/LinearEquationSolver.cpp
algebra/NumericVector.cpp
algebra/LinearEquationSolverPetsc.cpp
algebra/LinearEquationSolverPetscAsm.cpp
algebra/LinearEquationSolverPetscFieldSplit.cpp
algebra/PetscMatrix.cpp
algebra/PetscPreconditioner.cpp
algebra/PetscVector.cpp
algebra/Preconditioner.cpp
algebra/SparseMatrix.cpp
algebra/FunctionBase.cpp
algebra/ParsedFunction.cpp
algebra/SlepcSVD.cpp
fe/ElemType.cpp
fe/ElemType_template.cpp
fe/ElemType_template_base.cpp
fe/Hexaedron.cpp
fe/Line.cpp
fe/Quadrilateral.cpp
fe/Tetrahedral.cpp
fe/Triangle.cpp
fe/Wedge.cpp
fe_quadrature/GaussPoints.cpp
ism/Marker.cpp
ism/PolynomialBases.cpp
ism/Line.cpp
parallel/MyMatrix.cpp
parallel/MyVector.cpp
physics/Fluid.cpp
physics/Material.cpp
physics/Parameter.cpp
physics/Solid.cpp
uq/uq.cpp
uq/sparseGrid.cpp
)

# IF(HAVE_HDF5)
#   SET(femus_src ${femus_src} 
#   01_mesh/input/MED_IO.cpp
#   02_solution/01_output/XDMFWriter.cpp
#   01_mesh/gencase/GenCase.cpp
#   03_equations/TimeLoop.cpp
#   )
# ENDIF(HAVE_HDF5) 

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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(femus)

#############################################################################################
### Source files
#############################################################################################

SET(femus_src 
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
equations/assemble/Quantity.cpp
equations/assemble/DofMap.cpp
equations/assemble/Assemble_jacobian.cpp
equations/assemble/Assemble_unknown.cpp
equations/assemble/BoundaryConditions.cpp
equations/assemble/CurrentElem.cpp
equations/assemble/CurrentQuantity.cpp
equations/assemble/SystemTwo.cpp
equations/ExplicitSystem.cpp
equations/ImplicitSystem.cpp
equations/LinearImplicitSystem.cpp
equations/MonolithicFSINonLinearImplicitSystem.cpp
equations/NonLinearImplicitSystem.cpp
equations/NonLinearImplicitSystemWithPrimalDualActiveSetMethod.cpp
equations/MultiLevelProblem.cpp
equations/System.cpp
equations/TimeLoop.cpp
equations/TransientSystem.cpp
fe/ElemType.cpp
fe/ElemType_template.cpp
fe/ElemType_template_base.cpp
fe/Hexaedron.cpp
fe/Line.cpp
fe/Quadrilateral.cpp
fe/Tetrahedral.cpp
fe/Triangle.cpp
fe/Wedge.cpp
ism/Marker.cpp
ism/PolynomialBases.cpp
ism/Line.cpp
mesh/gencase/Box.cpp
mesh/gencase/Domain.cpp
mesh/gencase/ElemSto.cpp
mesh/gencase/MultiLevelMeshTwo.cpp
mesh/gencase/GenCase.cpp
mesh/geom_elements/GeomElemBase.cpp
mesh/geom_elements/GeomElemHex8.cpp
mesh/geom_elements/GeomElemHex27.cpp
mesh/geom_elements/GeomElemQuad4.cpp
mesh/geom_elements/GeomElemQuad9.cpp
mesh/geom_elements/GeomElemTet4.cpp
mesh/geom_elements/GeomElemTet10.cpp
mesh/geom_elements/GeomElemTri3.cpp
mesh/geom_elements/GeomElemTri6.cpp
mesh/geom_elements/GeomElemTri7.cpp
mesh/geom_elements/GeomElemEdge2.cpp
mesh/geom_elements/GeomElemEdge3.cpp
mesh/Elem.cpp
mesh/Mesh.cpp
mesh/MultiLevelMesh.cpp
mesh/obj_io.cpp
mesh/input/MeshGeneration.cpp
mesh/input/GambitIO.cpp
mesh/input/MED_IO.cpp
mesh/refinement/MeshRefinement.cpp
mesh/partitioning/MeshPartitioning.cpp
mesh/partitioning/MeshMetisPartitioning.cpp
mesh/partitioning/MeshASMPartitioning.cpp
parallel/MyMatrix.cpp
parallel/MyVector.cpp
physics/Fluid.cpp
physics/Material.cpp
physics/Parameter.cpp
physics/Solid.cpp
quadrature/GaussPoints.cpp
solution/MultiLevelSolution.cpp
solution/Solution.cpp
solution/output/Writer.cpp
solution/output/VTKWriter.cpp
solution/output/GMVWriter.cpp
solution/output/XDMFWriter.cpp
utils/FemusInit.cpp
utils/Files.cpp
utils/input_parser/InputParser.cpp
utils/input_parser/JsonInputParser.cpp
uq/uq.cpp
uq/sparseGrid.cpp
)

# IF(HAVE_HDF5)
#   SET(femus_src ${femus_src} 
#   mesh/input/MED_IO.cpp
#   solution/output/XDMFWriter.cpp
#   mesh/gencase/GenCase.cpp
#   equations/TimeLoop.cpp
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

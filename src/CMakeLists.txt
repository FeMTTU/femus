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
equations/DofMap.cpp
equations/BoundaryConditions.cpp
equations/CurrentElem.cpp
equations/CurrentGaussPoint.cpp
equations/CurrentGaussPointBase.cpp
equations/CurrentQuantity.cpp
equations/ExplicitSystem.cpp
equations/ImplicitSystem.cpp
equations/LinearImplicitSystem.cpp
equations/MonolithicFSINonLinearImplicitSystem.cpp
equations/NonLinearImplicitSystem.cpp
equations/NonLinearImplicitSystemWithPrimalDualActiveSetMethod.cpp
equations/MultiLevelProblem.cpp
equations/System.cpp
equations/SystemTwo.cpp
equations/TransientSystem.cpp
equations/FE_convergence.cpp
equations/Assemble_jacobian.cpp
fe/ElemType.cpp
fe/Hexaedron.cpp
fe/Line.cpp
fe/Quadrilateral.cpp
fe/Tetrahedral.cpp
fe/Triangle.cpp
fe/Wedge.cpp
ism/Marker.cpp
ism/PolynomialBases.cpp
ism/Line.cpp
meshGencase/Box.cpp
meshGencase/Domain.cpp
meshGencase/ElemSto.cpp
meshGencase/MultiLevelMeshTwo.cpp
meshGencase/GeomElemBase.cpp
meshGencase/GeomElemHex8.cpp
meshGencase/GeomElemHex27.cpp
meshGencase/GeomElemQuad4.cpp
meshGencase/GeomElemQuad9.cpp
meshGencase/GeomElemTet4.cpp
meshGencase/GeomElemTet10.cpp
meshGencase/GeomElemTri3.cpp
meshGencase/GeomElemTri6.cpp
meshGencase/GeomElemEdge2.cpp
meshGencase/GeomElemEdge3.cpp
mesh/Elem.cpp
mesh/Mesh.cpp
mesh/MultiLevelMesh.cpp
mesh/MeshGeneration.cpp
mesh/GambitIO.cpp
mesh/obj_io.cpp
mesh/MeshRefinement.cpp
mesh/MeshMetisPartitioning.cpp
mesh/MeshPartitioning.cpp
mesh/MeshASMPartitioning.cpp
parallel/MyMatrix.cpp
parallel/MyVector.cpp
physics/Fluid.cpp
physics/Material.cpp
physics/Parameter.cpp
physics/Solid.cpp
quadrature/GaussPoints.cpp
solution/MultiLevelSolution.cpp
solution/Quantity.cpp
solution/Solution.cpp
solution/Writer.cpp
solution/VTKWriter.cpp
solution/GMVWriter.cpp
utils/FemusInit.cpp
utils/Files.cpp
utils/InputParser.cpp
utils/JsonInputParser.cpp
uq/uq.cpp
uq/sparseGrid.cpp
)

IF(HAVE_HDF5)
  SET(femus_src ${femus_src} 
  mesh/MED_IO.cpp
  solution/XDMFWriter.cpp
  meshGencase/GenCase.cpp
  equations/TimeLoop.cpp
  )
ENDIF(HAVE_HDF5) 

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

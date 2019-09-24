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
equations/assemble/CurrentGaussPoint.cpp
equations/assemble/CurrentGaussPointBase.cpp
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
meshGencase/Box.cpp
meshGencase/Domain.cpp
meshGencase/ElemSto.cpp
meshGencase/MultiLevelMeshTwo.cpp
meshGencase/GenCase.cpp
mesh/geom_elements/GeomElemBase.cpp
mesh/geom_elements/GeomElemHex8.cpp
mesh/geom_elements/GeomElemHex27.cpp
mesh/geom_elements/GeomElemQuad4.cpp
mesh/geom_elements/GeomElemQuad9.cpp
mesh/geom_elements/GeomElemTet4.cpp
mesh/geom_elements/GeomElemTet10.cpp
mesh/geom_elements/GeomElemTri3.cpp
mesh/geom_elements/GeomElemTri6.cpp
mesh/geom_elements/GeomElemEdge2.cpp
mesh/geom_elements/GeomElemEdge3.cpp
mesh/Elem.cpp
mesh/Mesh.cpp
mesh/MultiLevelMesh.cpp
mesh/MeshGeneration.cpp
mesh/GambitIO.cpp
mesh/MED_IO.cpp
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
solution/Solution.cpp
solution/Writer.cpp
solution/VTKWriter.cpp
solution/GMVWriter.cpp
solution/XDMFWriter.cpp
utils/FemusInit.cpp
utils/Files.cpp
utils/InputParser.cpp
utils/JsonInputParser.cpp
uq/uq.cpp
uq/sparseGrid.cpp
)

IF (NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib64 CACHE PATH
       "Single output directory for building all libraries.")
  MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
ENDIF (NOT LIBRARY_OUTPUT_PATH)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${femus_src})


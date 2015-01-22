CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(femus)

#############################################################################################
### Source files
#############################################################################################

SET(femus_src 
algebra/AsmPetscLinearEquationSolver.cpp
algebra/DenseMatrixBase.cpp
algebra/DenseMatrix.cpp
algebra/DenseSubmatrix.cpp
algebra/DenseVectorBase.cpp
algebra/Graph.cpp
algebra/LinearEquation.cpp
algebra/LinearEquationSolver.cpp
algebra/NumericVector.cpp
algebra/GmresPetscLinearEquationSolver.cpp
algebra/PetscLinearSolverM.cpp
algebra/PetscMatrix.cpp
algebra/PetscPreconditioner.cpp
algebra/PetscVector.cpp
algebra/Preconditioner.cpp
algebra/SparseMatrix.cpp
algebra/VankaPetscLinearEquationSolver.cpp
algebra/FunctionBase.cpp
algebra/ParsedFunction.cpp
equations/DofMap.cpp
equations/CurrentElem.cpp
equations/CurrentGaussPoint.cpp
equations/CurrentGaussPointBase.cpp
equations/CurrentQuantity.cpp
equations/ExplicitSystem.cpp
equations/ImplicitSystem.cpp
equations/LinearImplicitSystem.cpp
equations/MonolithicFSINonLinearImplicitSystem.cpp
equations/NonLinearImplicitSystem.cpp
equations/MultiLevelProblem.cpp
equations/MultiLevelProblemTwo.cpp
equations/System.cpp
equations/SystemTwo.cpp
equations/TimeLoop.cpp
equations/TransientSystem.cpp
equations/NewmarkTransientSystem.cpp
fe/ElemType.cpp
fe/FEElemBase.cpp
fe/FEHex1.cpp
fe/FEHex27.cpp
fe/FEHex8.cpp
fe/FEQuad1.cpp
fe/FEQuad4.cpp
fe/FEQuad9.cpp
fe/FETet10.cpp
fe/FETet1.cpp
fe/FETet4.cpp
fe/FETri1.cpp
fe/FETri3.cpp
fe/FETri6.cpp
fe/FEEdge3.cpp
fe/FEEdge2.cpp
fe/FEEdge1.cpp
fe/Hexaedron.cpp
fe/Line.cpp
fe/Quadrilateral.cpp
fe/Tetrahedral.cpp
fe/Triangle.cpp
fe/Wedge.cpp
meshGencase/Box.cpp
meshGencase/Domain.cpp
meshGencase/ElemSto.cpp
meshGencase/GeomEl.cpp
meshGencase/MultiLevelMeshTwo.cpp
meshGencase/GenCase.cpp
mesh/Elem.cpp
mesh/Mesh.cpp
mesh/MultiLevelMesh.cpp
mesh/MeshGeneration.cpp
mesh/GambitIO.cpp
mesh/MeshRefinement.cpp
mesh/MeshMetisPartitioning.cpp
mesh/MeshPartitioning.cpp
mesh/MeshASMPartitioning.cpp
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
solution/XDMFWriter.cpp
utils/FemusInit.cpp
utils/Files.cpp
utils/IO.cpp
utils/InputParser.cpp
utils/JsonInputParser.cpp
)

IF (NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib64 CACHE PATH
       "Single output directory for building all libraries.")
  MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
ENDIF (NOT LIBRARY_OUTPUT_PATH)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${femus_src})


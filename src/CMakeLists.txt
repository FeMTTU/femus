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
algebra/LinearSolverM.cpp
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
equations/EqnBase.cpp
equations/EquationsMap.cpp
equations/ExplicitSystem.cpp
equations/ImplicitSystem.cpp
equations/LinearImplicitSystem.cpp
equations/MonolithicFSINonLinearImplicitSystem.cpp
equations/NonLinearImplicitSystem.cpp
equations/MultiLevelProblem.cpp
equations/Quantity.cpp
equations/QuantityLocal.cpp
equations/System.cpp
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
fe/Hexaedron.cpp
fe/Line.cpp
fe/Quadrilateral.cpp
fe/Tetrahedral.cpp
fe/Triangle.cpp
fe/Wedge.cpp
meshGencase/Box.cpp
meshGencase/CurrElem.cpp
meshGencase/Domain.cpp
meshGencase/ElemSto.cpp
meshGencase/GeomEl.cpp
meshGencase/MeshTwo.cpp
mesh/Elem.cpp
mesh/Mesh.cpp
mesh/MultiLevelMesh.cpp
mesh/MeshGeneration.cpp
mesh/GambitIO.cpp
physics/Fluid.cpp
physics/Material.cpp
physics/Parameter.cpp
physics/Physics.cpp
physics/Solid.cpp
quadrature/CurrGaussPoint.cpp
quadrature/CurrGaussPointBase.cpp
quadrature/GaussPoints.cpp
quadrature/QRule.cpp
solution/MultiLevelSolution.cpp
solution/Solution.cpp
solution/Writer.cpp
solution/VTKWriter.cpp
solution/GMVWriter.cpp
solution/XDMFWriter.cpp
utils/FemTTUInit.cpp
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


CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(femlib)

SET(femlib_src 
algebra/DenseMatrixBase.cpp
algebra/DenseMatrix.cpp
algebra/DenseSubmatrix.cpp
algebra/DenseVectorBase.cpp
algebra/Graph.cpp
algebra/LinearEquation.cpp
algebra/LinearEquationSolver.cpp
algebra/LinearSolverM.cpp
algebra/NumericVector.cpp
algebra/PetscLinearEquationSolver.cpp
algebra/PetscLinearSolverM.cpp
algebra/PetscMatrix.cpp
algebra/PetscPreconditioner.cpp
algebra/PetscVector.cpp
algebra/Preconditioner.cpp
algebra/SparseMatrix.cpp
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
equations/Solution.cpp
equations/System.cpp
equations/TimeLoop.cpp
equations/TransientSystem.cpp
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
fe/Hexaedron.cpp
fe/Line.cpp
fe/Quadrilateral.cpp
fe/Tetrahedral.cpp
fe/Triangle.cpp
fe/Wedge.cpp
mesh/Box.cpp
mesh/CurrElem.cpp
mesh/Domain.cpp
mesh/Elem.cpp
mesh/ElemSto.cpp
mesh/GeomEl.cpp
mesh/Mesh.cpp
mesh/MeshTwo.cpp
physics/Fluid.cpp
physics/Material.cpp
physics/Parameter.cpp
physics/Physics.cpp
physics/Solid.cpp
quadrature/CurrGaussPoint.cpp
quadrature/CurrGaussPointBase.cpp
quadrature/QRule.cpp
utils/FemTTUInit.cpp
utils/FemusInit.cpp
utils/Files.cpp
utils/Utils.cpp
)

IF (NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib64 CACHE PATH
       "Single output directory for building all libraries.")
  MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
ENDIF (NOT LIBRARY_OUTPUT_PATH)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${femlib_src})


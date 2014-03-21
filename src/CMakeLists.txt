CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(femlib)

SET(femlib_src 
mesh/Box.cpp
mesh/CurrElem.cpp
quadrature/CurrGaussPointBase.cpp
quadrature/CurrGaussPoint.cpp
algebra/DenseMatrixBase.cpp
algebra/DenseMatrix.cpp
algebra/DenseSubmatrix.cpp
algebra/DenseVectorBase.cpp
mesh/Domain.cpp
mesh/Elem.cpp
mesh/ElemSto.cpp
fe/ElemType.cpp
equations/EqnBase.cpp
equations/EquationsMap.cpp
equations/ExplicitSystem.cpp
fe/FEElemBase.cpp
fe/FEHex1.cpp
fe/FEHex27.cpp
fe/FEHex8.cpp
utils/FemTTUInit.cpp
utils/FemusInit.cpp
fe/FEQuad1.cpp
fe/FEQuad4.cpp
fe/FEQuad9.cpp
fe/FETet10.cpp
fe/FETet1.cpp
fe/FETet4.cpp
fe/FETri1.cpp
fe/FETri3.cpp
fe/FETri6.cpp
utils/Files.cpp
physics/Fluid.cpp
mesh/GeomEl.cpp
algebra/Graph.cpp
fe/Hexaedron.cpp
equations/ImplicitSystem.cpp
equations/LinearImplicitSystem.cpp
algebra/LinearSolver.cpp
algebra/LinearSolverM.cpp
fe/Line.cpp
equations/LinSysPde.cpp
physics/Material.cpp
mesh/mesh.cpp
mesh/Mesh.cpp
equations/NonLinearImplicitSystem.cpp
equations/MultiLevelProblem.cpp
equations/NonLinearTimeDependentMultiLevelProblem.cpp
algebra/NumericVector.cpp
physics/Parameter.cpp
algebra/PetscLinearSolver.cpp
algebra/PetscLinearSolverM.cpp
algebra/PetscMatrix.cpp
algebra/PetscPreconditioner.cpp
algebra/PetscVector.cpp
physics/Physics.cpp
algebra/Preconditioner.cpp
quadrature/QRule.cpp
fe/Quadrilateral.cpp
equations/Quantity.cpp
equations/QuantityLocal.cpp
physics/Solid.cpp
equations/Solution.cpp
algebra/SparseMatrix.cpp
equations/System.cpp
fe/Tetrahedral.cpp
equations/TimeLoop.cpp
equations/TransientSystem.cpp
fe/Triangle.cpp
utils/Utils.cpp
fe/Wedge.cpp
)

IF (NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib64 CACHE PATH
       "Single output directory for building all libraries.")
  MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
ENDIF (NOT LIBRARY_OUTPUT_PATH)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${femlib_src})


CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(femlib)

SET(femlibsrc 
Box.cpp
CurrElem.cpp
CurrGaussPointBase.cpp
CurrGaussPoint.cpp
DenseMatrixBase.cpp
DenseMatrix.cpp
DenseSubmatrix.cpp
DenseVectorBase.cpp
Domain.cpp
Elem.cpp
ElemSto.cpp
ElemType.cpp
EqnBase.cpp
EquationsMap.cpp
ExplicitSystem.cpp
FEElemBase.cpp
FEHex1.cpp
FEHex27.cpp
FEHex8.cpp
FemTTUInit.cpp
FemusInit.cpp
FEQuad1.cpp
FEQuad4.cpp
FEQuad9.cpp
FETet10.cpp
FETet1.cpp
FETet4.cpp
FETri1.cpp
FETri3.cpp
FETri6.cpp
Files.cpp
Fluid.cpp
GeomEl.cpp
Graph.cpp
Hexaedron.cpp
ImplicitSystem.cpp
LinearImplicitSystem.cpp
LinearSolver.cpp
LinearSolverM.cpp
Line.cpp
LinSysPde.cpp
Material.cpp
mesh.cpp
Mesh.cpp
NonLinearImplicitSystem.cpp
NonLinearMultiLevelProblem.cpp
NonLinearTimeDependentMultiLevelProblem.cpp
NumericVector.cpp
Parameter.cpp
PetscLinearSolver.cpp
PetscLinearSolverM.cpp
PetscMatrix.cpp
PetscPreconditioner.cpp
PetscVector.cpp
Physics.cpp
Preconditioner.cpp
QRule.cpp
Quadrilateral.cpp
Quantity.cpp
QuantityLocal.cpp
Solid.cpp
Solution.cpp
SparseMatrix.cpp
System.cpp
Tetrahedral.cpp
TimeLoop.cpp
TransientSystem.cpp
Triangle.cpp
Utils.cpp
Wedge.cpp
)

IF (NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib64 CACHE PATH
       "Single output directory for building all libraries.")
  MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
ENDIF (NOT LIBRARY_OUTPUT_PATH)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${femlibsrc})


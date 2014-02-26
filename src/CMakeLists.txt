CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(femlib)

SET(femlibsrc 
DenseMatrix.cpp
DenseMatrixBase.cpp
DenseSubmatrix.cpp
DenseVectorBase.cpp
# distributed_vector.cpp
Elem.cpp
ElemType.cpp
FemTTUInit.cpp
Fluid.cpp
Hexaedron.cpp
Line.cpp
LinearSolver.cpp
LinSysPde.cpp
Material.cpp
Mesh.cpp
NonLinearMultiLevelProblem.cpp
NonLinearTimeDependentMultiLevelProblem.cpp
NumericVector.cpp
Parameter.cpp
PetscLinearSolver.cpp
PetscMatrix.cpp
PetscPreconditioner.cpp
PetscVector.cpp
Preconditioner.cpp
Quadrilateral.cpp
Solid.cpp
Solution.cpp
SparseMatrix.cpp
Tetrahedral.cpp
Triangle.cpp
Wedge.cpp
)

IF (NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib64 CACHE PATH
       "Single output directory for building all libraries.")
  MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
ENDIF (NOT LIBRARY_OUTPUT_PATH)

ADD_LIBRARY(${PROJECT_NAME} SHARED ${femlibsrc})


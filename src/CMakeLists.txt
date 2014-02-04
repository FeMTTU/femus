CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(femlib)

SET(femlibsrc 
DenseMatrix.cpp
DenseMatrixBase.cpp
DenseSubmatrix.cpp
DenseVectorBase.cpp
distributed_vector.cpp
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
PetscRectangularMatrix.cpp
PetscVector.cpp
Preconditioner.cpp
Quadrilateral.cpp
Solid.cpp
Solution.cpp
SparseMatrix.cpp
SparseRectangularMatrix.cpp
Tetrahedral.cpp
Triangle.cpp
Wedge.cpp
)

#SET(CMAKE_BUILD_TYPE Release)
SET(CMAKE_BUILD_TYPE Debug)


ADD_LIBRARY(${PROJECT_NAME} SHARED ${femlibsrc})


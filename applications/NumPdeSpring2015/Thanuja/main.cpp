// #include "MultiLevelProblem.hpp"
// #include "TransientSystem.hpp"
// #include "NumericVector.hpp"
// #include "Fluid.hpp"
// #include "Parameter.hpp"

// #include "SparseMatrix.hpp"
// #include "VTKWriter.hpp"
// #include "NonLinearImplicitSystem.hpp"
// #include "FElemTypeEnum.hpp"
// #include "Files.hpp"

// using std::cout;
// using std::endl;


#include "FemusInit.hpp"
using namespace femus;

int main(int argc,char **args) {
  
  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);
 
  return 0;
}


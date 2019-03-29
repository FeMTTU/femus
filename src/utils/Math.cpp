#include "Math.hpp"
#include "MultiLevelSolution.hpp"
#include "System.hpp"


namespace femus {

   
 void UnknownLocal::initialize(const Unknown & unknown, const MultiLevelSolution * ml_sol, System * mlPdeSys) {
   
   Solname     = unknown._name;
   SolPdeIndex = mlPdeSys->GetSolPdeIndex(Solname.c_str());
   SolIndex    = ml_sol->GetIndex        (Solname.c_str());
   SolFEType   = ml_sol->GetSolutionType(SolIndex);
  
  }


}

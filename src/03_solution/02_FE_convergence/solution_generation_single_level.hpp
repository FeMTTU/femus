#ifndef __femus_FE_convergence_sol_generation_hpp__
#define __femus_FE_convergence_sol_generation_hpp__


#include "MultiLevelSolution.hpp"


namespace femus {

class MultiLevelProblem;
class MultiLevelMesh;



class Solution_generation_single_level {
    
public: 

// What this function does is to yield a MultiLevelSolution object at each level, to construct a level hierarchy to be used for error analysis 
virtual const MultiLevelSolution  run_on_single_level(MultiLevelProblem & ml_prob,
                                                      MultiLevelMesh & ml_mesh,
                                                      const unsigned lev,
                                                      const std::vector< Unknown > & unknowns,
                                                      const std::vector< Math::Function< double > * > &  exact_sol,
                                                      const MultiLevelSolution::InitFuncMLProb SetInitialCondition,
                                                      const MultiLevelSolution::BoundaryFuncMLProb  SetBoundaryCondition,
                                                  const bool equation_solve
                                                     ) const = 0;
                                                     
};


}


#endif

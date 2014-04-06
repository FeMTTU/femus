#include "MultiLevelProblem.hpp"
#include "MultiLevelMesh.hpp"
#include "RunTimeMap.hpp"
#include "main.hpp"

class MyMultiGrid: public MultiLevelProblem {
  
public:
  
   ///configuration / physics
  RunTimeMap<double>* _runtime_double;
  void AddParameters(RunTimeMap<double>* runtimein);
  
  unsigned _qty_idx[N_QTIES];
  unsigned   _qty_ncomps[N_QTIES];

  
  MyMultiGrid(MultiLevelMesh * mlmesh_in);
  double Error();
  
  void   BuildSparsityPattern();
  double IncreaseStep(unsigned i);
  
  double ComputeProductivityIndexNum(int bd);
  std::vector<double> ComputeProductivityIndexDen(unsigned nints);

  
  
};
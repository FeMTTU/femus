#include "MultiLevelProblem.hpp"
#include "MultiLevelMesh.hpp"
#include "FemusInputParser.hpp"
#include "main.hpp"

using namespace femus;

class MyMultiGrid: public MultiLevelProblem {
  
public:
  
   ///configuration / physics
  FemusInputParser<double>* _runtime_double;
  void AddParameters(FemusInputParser<double>* runtimein);
  
  unsigned _qty_idx[N_QTIES];
  unsigned   _qty_ncomps[N_QTIES];

  
  MyMultiGrid(MultiLevelMesh * mlmesh_in);
  double Error();
  
  void   BuildSparsityPattern();
  double IncreaseStep(unsigned i);
  
  double ComputeProductivityIndexNum(int bd);
  std::vector<double> ComputeProductivityIndexDen(unsigned nints);

  
  
};
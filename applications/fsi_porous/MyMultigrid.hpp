#include "MultiLevelProblem.hpp"

#include "RunTimeMapSingleton.hpp"
#include "main.hpp"

class MyMultiGrid: public MultiLevelProblem {
  
public:
  
   ///configuration / physics
  RunTimeMap<double>* _runtime_double;
  void AddParameters(RunTimeMap<double>* runtimein);
  
  const unsigned _qty_idx[N_QTIES];
  const int   _qty_ncomps[N_QTIES];

  
  MyMultiGrid(const unsigned short &igridn,const unsigned short &igridr,
              const char mesh_file[], const char GaussOrder[],
              const double Lref=1.);
  double Error();
  
  void   BuildSparsityPattern();
  double IncreaseStep(unsigned i);
  
  double ComputeProductivityIndexNum(int bd);
  std::vector<double> ComputeProductivityIndexDen(unsigned nints);

  
  
};
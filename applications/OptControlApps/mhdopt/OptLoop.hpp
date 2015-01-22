#ifndef __optloop_h__
#define __optloop_h__

#include "Typedefs.hpp"
#include "FemusInputParser.hpp"
#include "SystemTwo.hpp"
#include "TimeLoop.hpp"



namespace femus {


// Forward class
class Files;
class SystemTwo;


class OptLoop  : public TimeLoop {

public:


  OptLoop(Files& files_in);
 ~OptLoop();

  void optimization_loop(MultiLevelProblemTwo& e_map_in);

  void init_equation_data(const SystemTwo* eqn);
 
  //====data  
    std::vector<NumericVector *> _x_oldopt;  //old optimization step

};

//prototypes that can stay outside of a class
  double ComputeIntegral (const uint Level, const MultiLevelMeshTwo* mesh, const SystemTwo* eqn);

  int ElFlagControl(const std::vector<double> el_xm, const MultiLevelMeshTwo* mesh);



} //end namespace femus



#endif
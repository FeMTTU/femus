#ifndef __toptloop_h__
#define __toptloop_h__

#include "Typedefs.hpp"
#include "FemusInputParser.hpp"
#include "SystemTwo.hpp"
#include "TimeLoop.hpp"



namespace femus {


// Forward class
class Files;



class OptLoop  : public TimeLoop {

public:


  OptLoop(Files& files_in);


void optimization_loop(MultiLevelProblemTwo& e_map_in);



};

 double ComputeIntegral (const uint Level,const MultiLevelMeshTwo* mesh, const SystemTwo* eqn);

 double ComputeNormControl (const uint Level, const MultiLevelMeshTwo* mesh, const SystemTwo* eqn, const uint reg_ord );

 int ElFlagControl(const std::vector<double> el_xm, const MultiLevelMeshTwo* mesh);


} //end namespace femus



#endif
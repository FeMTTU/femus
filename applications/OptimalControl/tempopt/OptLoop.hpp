#ifndef __femus_tempopt_OptLoop_hpp__
#define __femus_tempopt_OptLoop_hpp__

#include "Typedefs.hpp"
#include "FemusInputParser.hpp"
#include "SystemTwo.hpp"
#include "TimeLoop.hpp"
#include "CurrentElem.hpp"
#include "CurrentQuantity.hpp"


namespace femus {


// Forward class
class Files;



class OptLoop  : public TimeLoop {

public:


  OptLoop(Files& files_in, const FemusInputParser<double> & map_in);


void optimization_loop( MultiLevelProblem & e_map_in );



};

 double ComputeIntegral (const uint Level,  const MultiLevelMeshTwo* mesh, const SystemTwo* eqn, const std::string output_time);

 double ComputeNormControl (const uint Level, const MultiLevelMeshTwo* mesh, const SystemTwo* eqn, const uint reg_ord );

 int ElFlagControl(const std::vector<double> el_xm, const MultiLevelMesh* mesh);

 void TempDesired(CurrentQuantity& myvect, const CurrentElem<double> & currelem);


 bool  SetBoundaryCondition(const MultiLevelProblem * ml_prob, const std::vector <double> &xp, const char name[], double &value, const int facename, const double time);

double  SetInitialCondition(const MultiLevelProblem * ml_prob, const std::vector <double> &xp, const char name[]);

} //end namespace femus



//********* SPACE DIMENSION ************

#define DIMENSION    2
//   #define DIMENSION    3
// **************************************


#endif

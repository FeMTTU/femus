#ifndef __optloop_h__
#define __optloop_h__

#include "Typedefs.hpp"
#include "RunTimeMap.hpp"
#include "EqnBase.hpp"
#include "TimeLoop.hpp"



namespace femus {


// Forward class
class Files;



class OptLoop  : public TimeLoop {

public:


  OptLoop(Files& files_in);


void optimization_loop(EquationsMap& e_map_in);



};



} //end namespace femus



#endif
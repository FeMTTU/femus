#ifndef __mgtimeloop_h__
#define __mgtimeloop_h__

#include "Typedefs.hpp"
#include "FemusInputParser.hpp"
#include "SystemTwo.hpp"

namespace femus {



// Forward class
class Files;

// ===============================================
//                  TimeLoop class
// ===============================================
class TimeLoop {

/*  protected:*///TODO make protected use get/set
public:
    // Data ---------------------------
    Files&              _files; ///< Utils pointer
    FemusInputParser<double>  _timemap; 

  uint      _t_idx_in;  //initial time step index
  double    _time_in;  //initial time absolute value
  uint      _t_idx_final;  
  double    _time_final; 

    uint     _curr_t_idx; 
    double   _curr_time;

  TimeLoop(Files& files_in);

  ~TimeLoop(){};
   
// print -------------------------------------------
   //this function is ok here because it doesn't involve the map 
   //of the EQUATIONS, it is just a print of a time sequence to a .xmf file
   void transient_print_xmf ( const uint t_idx_in,const uint t_idx_final, const uint nolevels) const;

 // i did it "static" so that it can be used regardless of the specific instantiation;
 // since it is static it cannot act on the class runtime map which is not static datum;
 // so i have to pass the "unconstrained" runtime map explicitly  
  static void check_time_par(FemusInputParser<double>&  time_in);

  /////< MG time step solver (backward Euler)
  double MGTimeStep(const uint iter, SystemTwo * eqn) const;   
  
  void OneTimestepEqnLoop(const uint delta_t_step_in, const MultiLevelProblemTwo & eqnmap) const;

  void TransientLoop(const MultiLevelProblemTwo & eqnmap);   //a standard transient loop in alphabetical order
  void TransientSetup(const MultiLevelProblemTwo & eqnmap);  //initialization of all the equations in the map

};


} //end namespace femus



#endif //  ---------  end header -------------------------

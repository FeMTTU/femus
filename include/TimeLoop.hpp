#ifndef __mgtimeloop_h__
#define __mgtimeloop_h__

// configure files  --------------------
#include "Typedefs_conf.hpp"

// Forward class
class Utils;
#include "RunTimeMap.hpp"

// ===============================================
//                  TimeLoop class
// ===============================================
class TimeLoop  {

/*  protected:*///TODO make protected use get/set
public:
    // Data ---------------------------
    Utils&              _utils; ///< Utils pointer
    RunTimeMap<double>  _timemap; 

  uint      _t_idx_in;  //initial time step index
  double    _time_in;  //initial time absolute value
  uint      _t_idx_final;  
  double    _time_final; 

    uint     _curr_t_idx; 
    double   _curr_time;

  
// public:
  // Constructor -------------------------------
//  TimeLoop(Utils& mgutils_in, EquationsMap& mgeqmap_in);///< Constructor
  TimeLoop(Utils& mgutils_in);///< Constructor

  ~TimeLoop(){}; ///< Destructor
   
// print -------------------------------------------
   //this function is ok here because it doesn't involve the map 
   //of the EQUATIONS, it is just a print of a time sequence to a .xmf file
   void transient_print_xmf ( const uint t_idx_in,const uint t_idx_final);

 // i did it "static" so that it can be used regardless of the specific instantiation;
 // since it is static it cannot act on the class runtime map which is not static datum;
 // so i have to pass the "unconstrained" runtime map explicitly  
static void check_time_par(RunTimeMap<double>&  time_in);

};

#endif //  ---------  end header -------------------------

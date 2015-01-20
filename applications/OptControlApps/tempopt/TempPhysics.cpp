#include <cmath>

//library
#include "EquationsMap.hpp"
#include "TimeLoop.hpp"
#include "MeshTwo.hpp"
#include "Box.hpp"

//application
#include "Temp_conf.hpp"
#include "TempPhysics.hpp"
#include "EqnT.hpp"

namespace femus {

TempPhysics::TempPhysics( RunTimeMap<double> & map_in):
  Physics(map_in) { 
  
  
} 

 // ========================================================
 //the Box must already be initialized here
 //it must receive an el_xm already in the REFBox frame

///AAA questa funzione NON lavora se tu fai solo DUE SUDDIVISIONI PER LATO e nolevels=1 !!!

 int TempPhysics::ElFlagControl(const std::vector<double> el_xm)  const {

  Box* box= static_cast<Box*>(_mesh->GetDomain());
   
     int el_flagdom=0;

///optimal control
  #if DIMENSION==2
//============== OUTFLOW     
//        if (   el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])
// 	   && el_xm[0] < 0.75*(box->_le[0] - box->_lb[0])
// 	   && el_xm[1] > 0.75*(box->_le[1] - box->_lb[1]) ) {
//                  el_flagdom=1;
//              }
//============== RIGHT
//      if (   el_xm[0] > 0.85*(box->_le[0] - box->_lb[0])
// 	   && el_xm[1] > 0.25*(box->_le[1] - box->_lb[1])
// 	   && el_xm[1] < 0.75*(box->_le[1] - box->_lb[1]) ) {
//                  el_flagdom=1;
//              }
//============== INFLOW
// facciamo 0.2-0.4; 0.4-0.6; 0.6-0.8     
     if (     el_xm[1] > 0.4*(box->_le[1] - box->_lb[1])
           && el_xm[1] < 0.6*(box->_le[1] - box->_lb[1])
	   && el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])
	   && el_xm[0] < 0.75*(box->_le[0] - box->_lb[0]) ) {
                 el_flagdom=1;
             }
// //============== CENTER
//      if (  /*   el_xm[1] > 0.04*(box->_le[1] - box->_lb[1])*/
//           /* &&*/ el_xm[1] < 0.06*(box->_le[1] - box->_lb[1])
// 	   && el_xm[0] > 0.4*(box->_le[0] - box->_lb[0])
// 	   && el_xm[0] < 0.6*(box->_le[0] - box->_lb[0]) ) {
//                  el_flagdom=1;
//              }
     
 #else
      if ( el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])  
	&& el_xm[0] < 0.75*(box->_le[0] - box->_lb[0]) 
	&& el_xm[1] > 0.75*(box->_le[1] - box->_lb[1])
	&& el_xm[2] > 0.25*(box->_le[2] - box->_lb[2]) 
	&& el_xm[2] < 0.75*(box->_le[2] - box->_lb[2]) ) {
	el_flagdom=1;
        }
 #endif

return el_flagdom; 
}



// ====================================================
/// This function sets the nondimensional groups
void TempPhysics::set_nondimgroups() {

  const double rhof   = _physrtmap.get("rho0");
  const double Uref   = _physrtmap.get("Uref");
  const double Lref   = _physrtmap.get("Lref");
  const double  muf   = _physrtmap.get("mu0");

  _pref = rhof*Uref*Uref;
  _Re  = (rhof*Uref*Lref)/muf;
  _Fr  = (Uref*Uref)/(9.81*Lref);
  _Pr=muf/rhof;    ///< Prandl Number


  return;
}


 //=================
    void TempPhysics::transient_loopPlusJ(EquationsMap & eqmap_in)  {

    //  parameters
    double         dt = eqmap_in._timeloop._timemap.get("dt");
    int    print_step = eqmap_in._timeloop._timemap.get("printstep");

    double curr_time = eqmap_in._timeloop._time_in;  //initialize current time
    
 for (uint curr_step = eqmap_in._timeloop._t_idx_in + 1; curr_step <= eqmap_in._timeloop._t_idx_final; curr_step++) {

   curr_time += dt;

#if DEFAULT_PRINT_TIME==1
      std::clock_t  start_time=std::clock();
#endif

       eqmap_in._timeloop._curr_t_idx = curr_step;
       eqmap_in._timeloop._curr_time  = curr_time;

      std::cout << "\n  ** Solving time step " << eqmap_in._timeloop._curr_t_idx
                << ", time = "                 << eqmap_in._timeloop._curr_time   << " ***" << std::endl;

       
	const uint delta_t_step = curr_step - eqmap_in._timeloop._t_idx_in;

      //  time step for each system, without printing (good)
      eqmap_in._timeloop.OneTimestepEqnLoop(delta_t_step,eqmap_in);

#if DEFAULT_PRINT_TIME==1
      std::clock_t    end_time=std::clock();
#endif 

      // print solution
      if (delta_t_step%print_step == 0) eqmap_in.PrintSol(curr_step,curr_time);   //print sol.N.h5 and sol.N.xmf
    

#if DEFAULT_PRINT_TIME==1
      std::clock_t    end_time2=std::clock();
      std::cout << " Time solver ----->= "   << double(end_time- start_time)/ CLOCKS_PER_SEC
                << " Time printing ----->= " << double(end_time2- end_time) / CLOCKS_PER_SEC <<
                std::endl;
#endif 

//=====functional evaluations=======

#if T_EQUATIONS==1
		EqnT* eqnT = static_cast<EqnT*>(eqmap_in.get_eqs("Eqn_T"));

		
     double J = 0.;
J = eqnT->ComputeIntegral    ( eqmap_in._mesh._NoLevels - 1);
J = eqnT->ComputeNormControl ( eqmap_in._mesh._NoLevels - 1,0 );
J = eqnT->ComputeNormControl ( eqmap_in._mesh._NoLevels - 1,1 );
//=====functional evaluations =======

#endif


    }   // end time loop

    return;
  }
  


} //end namespace femus


  


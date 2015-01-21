#include "OptLoop.hpp"

#include "EquationsMap.hpp"

namespace femus {
  
  
 OptLoop::OptLoop(Files& files_in): TimeLoop(files_in) { }

 //=================
    void OptLoop::optimization_loop(EquationsMap & eqmap_in)  {

    //  parameters
    int    print_step = _timemap.get("printstep");

    double curr_time = _time_in;  //initialize current time
    
 for (uint curr_step = _t_idx_in + 1; curr_step <= _t_idx_final; curr_step++) {

   curr_time += 1.;

#if DEFAULT_PRINT_TIME==1
      std::clock_t  start_time=std::clock();
#endif

       _curr_t_idx = curr_step;
       _curr_time  = curr_time;

      std::cout << "\n  ** Solving time step " << _curr_t_idx
                << ", time = "                 << _curr_time   << " ***" << std::endl;

       
	const uint delta_t_step = curr_step - _t_idx_in;

      //  time step for each system, without printing (good)
      OneTimestepEqnLoop(delta_t_step,eqmap_in);

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


  
  
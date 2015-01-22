// ==================================================================
//                  Class TimeLoop
// ==================================================================
// lib include
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>


#include "TimeLoop.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "Files.hpp"
#include "MeshTwo.hpp"
#include "FemusInputParser.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"

#include "paral.hpp"


namespace femus {




 void TimeLoop::check_time_par(FemusInputParser<double>& time_in) {
  
  if (time_in.get("initial_step") < 0.) {std::cout << " negative restart ;;;;;;;;;;;;;;;;;;;" << std::endl  ; abort();}
 
  return;
}



// ==================================================================
/// Constructor
 
 TimeLoop::TimeLoop(Files& files_in):
 _files(files_in),
 _timemap("TimeLoop",_files._output_path)   {

 //inizialize to zero   
   _t_idx_in  = 0;  
    _time_in  = 0.;  
 _t_idx_final = 0;
  _time_final = 0.; 
 _curr_t_idx  = 0;
  _curr_time  = 0.;
   
}



// =================================================================
/// Xdmf transient  print
//print from time t_idx_in to t_idx_final
//I will print a separate time sequence for each LEVEL
//TODO see if there is a way to read multiple sol collections, defined on different grids, in the SAME time file
      void TimeLoop::transient_print_xmf(const uint t_idx_in,const uint t_idx_final, const uint nolevels_in) const {

	//multigrid
	uint NoLevels = nolevels_in;

        // time parameters
        const uint ndigits     = DEFAULT_NDIGITS;
        const int print_step   = _timemap.get("printstep");
        // dir names
        std::string    basetime     = DEFAULT_BASETIME;
        std::string    ext_xdmf     = DEFAULT_EXT_XDMF;
        std::string    basesol      = DEFAULT_BASESOL;
        std::string    aux_xdmf     = DEFAULT_AUX_XDMF;

// =================================
// ============= LEVELS ============
// =================================
	
	for (uint l=0; l < NoLevels; l++) {

	// file 
        std::ostringstream Name;
        Name << _files._output_path << "/" << basetime << "."
             << std::setw(ndigits) << std::setfill('0') << t_idx_in << "-"
             << std::setw(ndigits) << std::setfill('0') << t_idx_final  << "_l" << l
             << ext_xdmf;
        std::ofstream out(Name.str().c_str());
	if (out.fail()) {std::cout << "transient_print_xmf: cannot print timeN.xmf" << std::endl; abort();}
	
        uint nprt=1;
        std::string gname[3];  gname[0]=basesol;

        //   Mesh -----------------------------------
        out << "<?xml version=\"1.0\" ?> \n";
        out << "<!DOCTYPE Xdmf SYSTEM "
        <<  "\"" << _files._output_path << "/" << aux_xdmf << "\"" << "[]>\n";
        out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> \n";
        out << "<Domain> \n";
        for (uint kp=0;kp< nprt; kp++)    {
          out << "<Grid Name=\""<< gname[kp].c_str() <<"\"  GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
          // time loop for grid sequence
          for (uint it = t_idx_in; it <= t_idx_final; it++) if (it%print_step ==0)   {
              out << "<xi:include href=\""
                  << basesol << "." << std::setw(ndigits) << std::setfill('0') << it << "_l" << l <<  ext_xdmf 
                  << "\"" << " xpointer=\"xpointer(//Xdmf/Domain/Grid["<< kp+1 <<"])\" >\n";
              out << "<xi:fallback />\n";
              out << " </xi:include>\n";
            }
          out << "</Grid> \n";
        }
        // Grid Collection end
        out << "</Domain> \n";
        out << "</Xdmf> \n";
        out.close();
	
	} //end levels
	
      return;
    }

    
// ======================================================
/// This function controls the time step operations:
//first, prepare the A matrices for all levels
// and the b rhs for the fine level
//then call the multigrid solver
//eventually, update the old solution

//this is the NONLINEAR loop
//it is settled as a multigrid, but for 0 levels it is a normal single grid
//it returns the l2 norm of the difference between two nonlinear steps, x - x_old

//time is passed only because the matrices and rhs can depend on time
//   its value is given by TimeLoop

/// AAAAAAA TODO I have to use this for the PENALTY METHOD!!!
/// In general, it is safer to call both of them because
/// you may fill both Ke and Fe in both and they are needed at all levels
//     GenMatRhs(time,Level,0); // matrix
//     GenRhsB(time,Level);
//yes with NESTED! //with nonnested the rhs is overwritten,you're solving for another problem in the other levels
//with non-nested you do not fill the rhs (mode=0) nor you fill with the boundary integrals (GenRhsB commented)
//pay attention to this!
//i dont wanna use mode any longer!
//mode is used to avoid filling the RHS
//because the multigrid has another RHS at the non-fine levels.
//But actually the RHS is AUTOMATICALLY OVERWRITTEN for the non-fine levels,
//because it is always the RESTRICTION of the RESIDUAL from the UPPER LEVEL

//TODO Pay attention here, we are filling BOTH A and B at ALL LEVELS,
// But then in the multigrid algorithm b for some levels will be OVERWRITTEN!!!
//The point is that you still have to fill A, also with the boundary conditions if necessary (SEE PENALTY DIRICHLET BC)
// Of course, the boundary conditions must be satisfied at all levels

/// Basically, the single time step consists in  ASSEMBLING + SOLVING

double TimeLoop::MGTimeStep(const uint iter, SystemTwo * eqn_in) const {

    std::cout  << std::endl << " Solving " << eqn_in->_eqname << " , step " << iter << std::endl;

    ///A0) Put x_old into x_oold
    *(eqn_in->_x_oold[eqn_in->_NoLevels-1]) = *(eqn_in->_x_old[eqn_in->_NoLevels-1]);

    /// A) Assemblying 
#if  DEFAULT_PRINT_TIME==1
    std::clock_t start_time=std::clock();
#endif

    for (uint Level = 0 ; Level < eqn_in->_NoLevels; Level++) {

        eqn_in->_A[Level]->zero();
        eqn_in->_b[Level]->zero();

        eqn_in->GenMatRhs(Level);

#ifdef DEFAULT_PRINT_INFO
        eqn_in->_A[Level]->close();
        double ANorm = eqn_in->_A[Level]->l1_norm();
	
//	_A[Level]->print_graphic(true); TODO should pass this true or false as a parameter
	
        std::cout << " ANorm l1 " << Level << " "  << ANorm  << std::endl;
#endif
    }

#if    DEFAULT_PRINT_TIME==1
    std::clock_t end_time=std::clock();
    std::cout << " ================ Assembly time = " << double(end_time- start_time) / CLOCKS_PER_SEC
              << " s " << std::endl;
#endif

///std::cout << " $$$$$$$$$ Prepared A for all levels and b for the fine level $$$$$$" << std::endl;
// The matrices R and P for all levels were already prepared at read-time
// and they do not depend on  time (because no adaptive refinement is performed)
// They only depend on the nodes (dofs!)

/// C) Reinitialization of  the ILU (if ILUINIT defined),
// #ifdef ILUINIT
//   if (iter%ILUINIT==0) {
//     std::cout << " ILU Reinitialized: " << std::endl;
//     for (int Level = NoLevels - 1; Level > 0; Level--) *L[Level].ILUExists = _LPFalse;
//   }
// #endif

    /// D) Solution of the linear MGsystem
#ifdef DEFAULT_PRINT_TIME
        std::clock_t start_time_sol = std::clock();
#endif
      
        eqn_in->MGSolve(DEFAULT_EPS_LSOLV, DEFAULT_MAXITS_LSOLV);
	
// //     for (uint Level = 0 ; Level < _NoLevels; Level++)  { 
// //       _solver[Level]->solve(*_A[Level],*_x[Level],*_b[Level],1.e-6,40);
// //            _x[Level]->localize(*_x_old[Level]);   // x_old = x
// //     }

#if    DEFAULT_PRINT_TIME==1
    std::clock_t end_time_sol = std::clock();
    std::cout << " ================ Solver time = " << double(end_time_sol- start_time_sol) / CLOCKS_PER_SEC
              << " s "<< std::endl;
#endif

      
/// std::cout << "$$$$$$$$$ Computed the x with the MG method $$$$$$$" << std::endl;

    /// E) Update of the old solution at the top Level
    eqn_in->_x[eqn_in->_NoLevels-1]->localize(*(eqn_in->_x_old[eqn_in->_NoLevels-1]));   // x_old = x
#ifdef DEFAULT_PRINT_INFO
    std::cout << "$$$$$$$$$ Updated the x_old solution $$$$$$$$$" << std::endl;
#endif
/// std::cout << "$$$$$$$$$ Check the convergence $$$$$$$" << std::endl;

    eqn_in->_x_tmp[eqn_in->_NoLevels-1]->zero();
    eqn_in->_x_tmp[eqn_in->_NoLevels-1]->add(+1.,*(eqn_in->_x_oold[eqn_in->_NoLevels-1]));
    eqn_in->_x_tmp[eqn_in->_NoLevels-1]->add(-1.,*(eqn_in->_x_old[eqn_in->_NoLevels-1]));
    // x_oold -x_old =actually= (x_old - x)
    //(x must not be touched, as you print from it)
    //x_oold must not be touched ! Because it's used later for UPDATING Becont!
    //so you must create a temporary vector necessarily.

    eqn_in->_x_tmp[eqn_in->_NoLevels-1]->close();
    double deltax_norm = eqn_in->_x_tmp[eqn_in->_NoLevels-1]->l2_norm();
    std::cout << " $$$$$$ " << eqn_in->_eqname << " error l2 " << deltax_norm << std::endl;
    std::cout << " $$$$$$ " << eqn_in->_eqname << " error linfty " << eqn_in->_x_tmp[eqn_in->_NoLevels-1]->linfty_norm() << std::endl;
//AAA when the vectors have nan's, the norm becomes zero!
//when the residual norm in pre and post smoothing is too big,
//then it doesnt do any iterations, the solver doesnt solve anymore, so the solution remains frozen

    return deltax_norm;  //TODO do we have to be based on l2norm or linfty norm???
}
    
    
// ==========================================================================================
/// This function performes all the Physics time step routines
void TimeLoop::OneTimestepEqnLoop(
    const uint delta_t_step_in,     // integer time
    const MultiLevelProblemTwo & eqnmap) const {
    // loop for time steps
    for (MultiLevelProblemTwo::const_iterator eqn = eqnmap.begin(); eqn != eqnmap.end(); eqn++)  {
        SystemTwo* equation = eqn->second;
        MGTimeStep(delta_t_step_in,equation);
    }
    return;
}


////////////////////////////////
////////////////////////////////

void TimeLoop::TransientSetup(const MultiLevelProblemTwo & eqnmap)  {

    const uint initial_step = _timemap.get("initial_step");
    const uint ndigits      = DEFAULT_NDIGITS;

    std::string   lastrun_f = DEFAULT_LAST_RUN;
    std::string     basesol = DEFAULT_BASESOL;
    std::string    ext_xdmf = DEFAULT_EXT_XDMF;
    std::string      ext_h5 = DEFAULT_EXT_H5;

    std::string    basecase = DEFAULT_BASECASE;
    std::string    basemesh = DEFAULT_BASEMESH;

    std::string  aux_xdmf   = DEFAULT_AUX_XDMF;
    std::string  connlin    = DEFAULT_CONNLIN;


//now, every run, restart or not, has a new output dir.
//So, if you restart, you have to copy sol.N.h5 and sol.N.xmf
//to the new output directory
//just a raw copy, nothing more, because the file paths in sol.N.xmf
// DO NOT DEPEND ON THE OUTPUTDIR, only on the INPUT_DIR for now.
//the problem is that you have to know the PREVIOUS output_dir
// to automatically copy the files...
//let us just try by hand now...no we cant, because we will not know
// the NEW output_dir...
//We have to find a way to keep track of the PREVIOUS output dir.
//--- use a shell variable. $femus_last_run
//SORRY, BUT YOU CANT DO THAT!
// A process cannot export to parent processes.
// So the solution will be to print a very small file, called
// femus_last_run,that contains the name of the last output dir.

//We know the NEW one at this point because we constructed it.
// Well, we put them in the main output and read from there.
//So you see once more that it is important to make the difference
// between the BASE_OUTPUT and the OVERALL_OUTPUT paths
//We will distinguish them later.

    //------initial data
    if (_files._restart_flag) {
        _t_idx_in  = initial_step;
        std::cout << "We wish to restart from time step " << _t_idx_in << std::endl;
        std::cout << "\n *+*+* TimeLoop::transient_setup: RESTART  " << std::endl;


        if (paral::get_rank() == 0) {

            stringstream tidxin;
            tidxin << setw(ndigits) << setfill('0') << _t_idx_in;
            std::cout << " Restarting from run: " << _files._input_path << std::endl;
	    std::ostringstream cp_src_xmf_stream;
	    cp_src_xmf_stream << _files._input_path << "/" << basesol << "." << tidxin.str() << "_l" << (eqnmap._mesh._NoLevels - 1) << ext_xdmf;
            std::string cp_src_xmf = cp_src_xmf_stream.str();
            fstream file_cp_xmf(cp_src_xmf.c_str());
            if (!file_cp_xmf.is_open()) {
                std::cout << "No xmf file" << std::endl;
                abort();
            }

	    std::ostringstream cp_src_h5_stream;
	    cp_src_h5_stream << _files._input_path << "/" << basesol << "." << tidxin.str() /*<< "_l" << (_mesh._NoLevels - 1)*/ << ext_h5;
            std::string cp_src_h5 = cp_src_h5_stream.str();
            fstream file_cp_h5(cp_src_h5.c_str());
            if (!file_cp_h5.is_open()) {
                std::cout << "No h5 file" << std::endl;
                abort();
            }

            std::string cp_dest_dir = _files._output_path;

            std::string cp_cmd_xmf = "cp " + cp_src_xmf  + " " +  cp_dest_dir;
            std::string cp_cmd_h5  = "cp " + cp_src_h5   + " " +  cp_dest_dir;

            std::cout << "Copying the two restart files to the NEW output directory created before" << std::endl;
            //you should first check that sol.N.xmf and sol.N.h5 are there

            std::cout <<cp_cmd_xmf	 << std::endl;
            std::cout <<cp_cmd_h5  << std::endl;
            system(cp_cmd_xmf.c_str() );
            system(cp_cmd_h5.c_str() );

	    //TODO why did I not the CopyFile function from Files?
        }

//here you should wait so that you are sure that sol.N.xmf and sol.N.h5 have been copied to the new directory
        std::cout << "***** Barrier so that all the processors have the sol.N.xmf and sol.N.h5 to read from" << std::endl;
#ifdef HAVE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif

//now the question arises: if to restart,you need to re-read the sol. files,
//you will also need to
// - COPY the mesh.h5 and mesh_conn_lin.h5 (you can avoid regenerating them) (and also mesh.xmf to be able to read them)
// - COPY the MG OPERATORS (if you have them in your run) if you keep the same number of levels
//The rule is: what is connected to the PROBLEM to be solved (DOMAIN, BC's) should be:
        // LEVEL INDEPENDENT
        //PROCESSOR INDEPENDENT
//An appropriate RESTART of the solution of a problem should be
//INDEPENDENT of the METHOD to SOLVE the PROBLEM (multigrid or not, parallel computing or not...)
//So the files for the DOMAIN and the BOUNDARY CONDITIONS should reflect that.
// e.g. consider mesh.h5 and case.h5
//given some FINE DISCRETIZATION (I dont consider the case of RESTARTING with a DIFFERENT FINE discretization...
//then during the run one might do adaptive mesh, but this is another thing... ),
//some parts of them are good for a SPECIFIC NOLEVELS or a SPECIFIC NO_PROCESSORS
//while others, IN PARTICULAR THOSE THAT ARE NECESSARY FOR RESTART, must be INDEPENDENT OF THAT.

        eqnmap.ReadSol(_t_idx_in,_time_in); //read  sol.N.h5 and sol.N.xmf
        //AAA: here _t_idx_in is intent in, _time_in is intent out
        //reading files can be done in parallel with no problems
        //well, not really... reading can be done if you are sure that the file is there at the moment you call to read it
        //so let's put this inside procid=0 .... NO, THIS READ IS DONE IN PARALLEL!

    }


    else {
        //Set up initial time step index and time value
        _t_idx_in = 0;                            //time step index
        _time_in = 0.;                           //time absolute value

        eqnmap.PrintSol(_t_idx_in,_time_in);  //print sol.0.h5 and sol.0.xmf
        //AAA: here _t_idx_in is intent-in, and also _time_in is intent-in
    }

    std::cout << "\nInitial time index: " << _t_idx_in
              << ", Initial time value: " << _time_in << std::endl;

    //-------final data
    const int nsteps       = _timemap.get("nsteps");
    const double dt        = _timemap.get("dt");

    _t_idx_final = _t_idx_in + nsteps;
    _time_final  = _time_in + nsteps*dt;

    std::cout << "\nFinal time index: " << _t_idx_final
              << ", Final time value: " << _time_final
              << std::endl;

    //now you can update last_run with new_run for a following run
    //well,actually before putting the last_run you should be sure that this run was completely finished.
    //That is why I'd better put this call at the end of the main program
//   _utils._files.PrintRunForRestart(DEFAULT_LAST_RUN);

//------- print
    //this happens when the output dir is already set
    //at this point this is already true
    eqnmap.PrintCase(_t_idx_in);       //print caseN.xmf&h5 = IC + BC flags
    transient_print_xmf(_t_idx_in,_t_idx_final,eqnmap._mesh._NoLevels); //print timeN.xmf

    return;
}



///////////////////////////////////////////////////////
/*standard time loop over the map equations.
 The equations are solved in alphabetical order given by the map*/
void TimeLoop::TransientLoop(const MultiLevelProblemTwo & eqnmap)  {

    //  parameters
    double         dt = _timemap.get("dt");
    int print_step =    _timemap.get("printstep");

    double curr_time = _time_in;  //initialize current time

    for (uint curr_step = _t_idx_in + 1; curr_step <= _t_idx_final; curr_step++) {

        curr_time += dt;
        _curr_time  = curr_time;

#if DEFAULT_PRINT_TIME==1 // only for cpu time check --------
        std::clock_t  start_time=std::clock();
#endif // ------------------------------------------- 

        // set up the time step
        std::cout << "\n  ** Solving time step " << curr_step
                  << ", time = "                 << curr_time  << " ***" << std::endl;

        const uint delta_t_step = curr_step - _t_idx_in;


        //  time step for each system, without printing (good)
        OneTimestepEqnLoop(delta_t_step,eqnmap);

#if DEFAULT_PRINT_TIME==1 // only for cpu time check --------
        std::clock_t    end_time=std::clock();
#endif  // ------------------------------------------

        // print solution
        if (delta_t_step%print_step == 0) eqnmap.PrintSol(curr_step,curr_time);   //print sol.N.h5 and sol.N.xmf


#if DEFAULT_PRINT_TIME==1 // only for cpu time check --------
        std::clock_t    end_time2=std::clock();
        std::cout <<" Time solver ----->= "   << double(end_time- start_time)/ CLOCKS_PER_SEC
                  <<" Time printing ----->= " << double(end_time2- end_time) / CLOCKS_PER_SEC <<
                  std::endl;
#endif  // ------------------------------------------


    }   // end time loop

    return;
}

    
    
} //end namespace femus


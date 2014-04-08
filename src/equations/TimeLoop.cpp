// ==================================================================
//                  Class TimeLoop
// ==================================================================
// lib include
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath> //for fabs


// local include
#include "TimeLoop.hpp"

// class include
#include "Files.hpp"
#include "RunTimeMap.hpp"

#include "paral.hpp"


 void TimeLoop::check_time_par(RunTimeMap<double>& time_in) {
  
  if (time_in.get("restart") < 0.) {std::cout << " negative restart ;;;;;;;;;;;;;;;;;;;" << std::endl  ; abort();}
 
  return;
}



// ==================================================================
/// Constructor
 
 TimeLoop::TimeLoop(Files& files_in):
 _files(files_in),
 _timemap("TimeLoop",_files.get_basepath())   {

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
        const uint ndigits     = _timemap.get("ndigits");
        const int print_step   = _timemap.get("printstep");
        // dir names
        std::string    basepath     = _files.get_basepath();
        std::string    output_dir   = DEFAULT_OUTPUTDIR;
        std::string    outtime_dir  = _files.get_frtmap().get("OUTTIME_DIR");
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
        Name << basepath << "/" << output_dir << outtime_dir << basetime << "."
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
        <<  "\"" << basepath << "/" << output_dir << outtime_dir << "/" << aux_xdmf << "\"" << "[]>\n";
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

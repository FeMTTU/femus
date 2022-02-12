// run_advection_benchmark.cpp - Main benchmark program for advection benchmark

// This code is a comment-free zone - sorry about this...

#include <cmath>
#include <vector>

#include "Timer.h"
#include "adept.h"

/* This file can access four different automatic differentiation
   libraries by including a header file three times, each time with
   "areal" defined differently */

#define adouble double
#include "advection_schemes.h"
#undef adouble

#ifdef CPPAD
#include "cppad/cppad.hpp"
using CppAD::AD;
#define adouble AD<double>
#include "advection_schemes.h"
#undef adouble
typedef AD<double> cppad_d;
#endif

#ifdef SACADO
#include "Sacado.hpp"
template<> int Sacado::Rad::ADmemblock<double>::n_blocks = 0;
#define adouble Sacado::Rad::ADvar<double>
#include "advection_schemes.h"
#undef adouble
typedef Sacado::Rad::ADvar<double> sacado_d;
#endif

#ifdef SACADO_FAD
#define adouble Sacado::ELRFad::DFad<double>
#include "advection_schemes.h"
#undef adouble
typedef Sacado::ELRFad::DFad<double> sacado_fad_d;
#endif

#define adouble adept::adouble
#include "advection_schemes.h"
#undef adouble

#ifdef ADOLC
#include "adolc/adolc.h"
#include "advection_schemes.h"
typedef adouble adolc_d;
#endif

#include "advection_schemes_AD.h"

typedef adept::adouble adept_d;

#define DEBUG std::cerr << "DEBUG line " << __LINE__ << "\n"

Timer timer;

int
main(int argc, char** argv)
{


  double pi = 4.0*atan(1.0);
  //  int NX = 102;
  //  int NX = 18;
  double q_init[NX];
  double q[NX];
  double q_AD[NX];
  double q_init_AD[NX];
  //  for (int i = NX/4; i < NX*3/4; i++) q_init[i] = 1.0;

  int nt = 100;
  int nr = 100;
  double dt = 0.125;

  bool is_lax_wendroff = true;
  bool verbose = false;

  std::string scheme;

  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " SCHEME [NUM_TIMESTEPS [NUM_REPEATS [COURANT_NUMBER [jacobian]]]]\n";
    std::cout << "  where SCHEME is either \"lax_wendroff\" or \"toon\"\n"
	      << "        NUM_TIMESTEPS has default value " << nt << "\n"
	      << "        NUM_REPEATS has default value " << nr << "\n"
	      << "        COURANT_NUMBER has default value " << dt << "\n";
    exit(1);
  }
  else {
    scheme = argv[1];
    if (scheme == "lax_wendroff") {
      is_lax_wendroff = true;
    }
    else if (scheme == "toon") {
      is_lax_wendroff = false;
    }
    else {
      std::cout << "SCHEME must be either \"lax_wendroff\" or \"toon\"\n";
      exit(1);      
    }
    if (argc > 2) {
      nt = atoi(argv[2]);
      if (argc > 3) {
	nr = atoi(argv[3]);
	if (argc > 4) {
	  dt = atof(argv[4]);
	}
      }
    }
  }

  std::cout << "------------------------------\n";
  std::cout << "Scheme: " << scheme << "\n";
  std::cout << "Number of timesteps: " << nt << "\n";
  std::cout << "Number of repeats: " << nr << "\n";
  std::cout << "Courant number: " << dt << "\n";

  
  int do_all = 1;
  bool do_adept = 0;
  bool do_cppad = 0;
  bool do_sacado = 0;
  bool do_sacado_fad = 0;
  bool do_adolc = 0;
  bool do_jacobian = 0;
  int force_jacobian = 0; // -1 to force reverse, +1 to force forward

  if (argc > 5) {
    if (std::string(argv[5]) == "adept") {
      do_adept = 1; do_all = 0;
    }
    if (std::string(argv[5]) == "adept-jacobian") {
      do_adept = 1; do_all = 0; do_jacobian = true;
    }
    else if (std::string(argv[5]) == "adolc") {
      do_adolc = 1; do_all = 0;
    }
    else if  (std::string(argv[5]) == "sacado") {
      do_sacado = 1; do_all = 0;
    }
    else if  (std::string(argv[5]) == "sacado-jacobian") {
      do_sacado_fad = 1; do_all = 0;
    }
    else if  (std::string(argv[5]) == "adolc-jacobian") {
      do_adolc = 1; do_all = 0; do_jacobian = true;
    }
    else if  (std::string(argv[5]) == "cppad-jacobian") {
      do_cppad = 1; do_all = 0; do_jacobian = true;
    }
    else if  (std::string(argv[5]) == "cppad") {
      do_cppad = 1; do_all = 0;
    }
    else if  (std::string(argv[5]) == "jacobian") {
      do_cppad = 1; 
      do_adept = 1;
      do_adolc = 1;
      do_sacado_fad = 1;
      do_jacobian = true; do_all = 0;
    }
    else if (std::string(argv[5]) == "forward-jacobian") {
      do_cppad = 1; 
      do_adept = 1;
      do_adolc = 1;
      do_sacado_fad = 1;
      do_jacobian = true; do_all = 0;
      force_jacobian = +1;
    }
    else if (std::string(argv[5]) == "reverse-jacobian") {
      do_cppad = 1; 
      do_adept = 1;
      do_adolc = 1;
      do_sacado_fad = 1;
      do_jacobian = true; do_all = 0;
      force_jacobian = -1;
    }
    else {
      std::cout << argv[5] << " not understood\n";
      exit(1);
    }
  }

  if (do_jacobian) {
    if (force_jacobian > 0) {
      std::cout << "Computing Jacobian: forcing forward mode\n";
    }
    else if (force_jacobian < 0) {
      std::cout << "Computing Jacobian: forcing reverse mode\n";
    }
    else {
      std::cout << "Computing Jacobian\n";
    }
  }
  std::cout << "------------------------------\n";

  for (int i = 0; i < NX; i++) q_init[i] = (0.5+0.5*sin((i*2.0*pi)/(NX-1.5)))+0.0001;
  for (int i = 0; i < NX; i++) q_AD[i] = 0.0;
  q_AD[NX/2]=1.0;

  //  lax_wendroff(100, 0.5, q_init, q);

  if (verbose) {
    std::cout << "Initial values: ";
    for (int i = 0; i < NX; i++) { std::cout << " " << q_init[i]; }
    std::cout << "\n";
  }


  int original_id = timer.new_activity("Original function");

  for (int j = 0; j < nr; j++)
  {
    timer.start(original_id);
    if (is_lax_wendroff) {
      lax_wendroff(nt, dt, q_init, q);
    }
    else {
      toon(nt, dt, q_init, q);
    }
    timer.stop();
    if (j == nr-1) {
      std::cout << "Original function\n";
      if (verbose) {
	std::cout << "Values:";
	for (int i = 0; i < NX; i++) { std::cout << " " << q[i]; }
	std::cout << "\n";
      }
    }
  }



  int hand_coded_id = timer.new_activity("Hand-coded adjoint");

  for (int j = 0; j < nr; j++)
  {
    timer.start(hand_coded_id);
    if (is_lax_wendroff) {
      lax_wendroff_AD(nt, dt, q_init, q, q_AD, q_init_AD);
    }
    else {
      toon_AD(nt, dt, q_init, q, q_AD, q_init_AD);
    }
    timer.stop();
  
    if (j == nr-1) {
      std::cout << "Hand-coded adjoint\n";
      if (verbose) {
	std::cout << "Values:";
	for (int i = 0; i < NX; i++) { std::cout << " " << q[i]; }
	std::cout << "\nAdjoints";
	for (int i = 0; i < NX; i++) { std::cout << " " << q_init_AD[i]; }
	std::cout << "\n";
      }
    }

  }
 

  int adept_record_id = timer.new_activity("Adept record");
  int adept_adjoint_id = timer.new_activity("Adept adjoint");
  int adept_jacobian_id = timer.new_activity("Adept Jacobian");

  if (do_all || do_adept) 
  {
    bool jac_is_printed = false;
    adept::Stack adept_stack;
      std::cout << "Adept\n";

  for (int j = 0; j < nr; j++)
  {
    //    adept_stack.start(1000);
    adept_d adept_q_init[NX];
    adept_d adept_q[NX];
    adept::set_values(adept_q_init, NX, q_init);

    adept_stack.new_recording();

    timer.start(adept_record_id);
    if (is_lax_wendroff) {
      lax_wendroff(nt, dt, adept_q_init, adept_q);
    }
    else {
      toon(nt, dt, adept_q_init, adept_q);
    }

    if (!do_jacobian) {
    timer.start(adept_adjoint_id);

    set_gradients(adept_q, NX, q_AD);
    adept_stack.compute_adjoint();
    adept::get_gradients(adept_q_init, NX, q_init_AD);
    timer.stop();

    if (verbose) {
      if (j == nr-1) {
	std::cout << "Values:";
	for (int i = 0; i < NX; i++) { std::cout << " " << adept_q[i]; }
	std::cout << "\nAdjoints";
	for (int i = 0; i < NX; i++) { std::cout << " " << q_init_AD[i]; }
	std::cout << "\n";
	std::cout << adept_stack;
      }
    }


    }
    else {
      adept_stack.independent(adept_q_init, NX);
      adept_stack.dependent(adept_q, NX);
      double jacobian[NX*NX];
      timer.start(adept_jacobian_id);
      if (force_jacobian > 0) {
	adept_stack.jacobian_forward(jacobian);
      }
      else if (force_jacobian < 0) {
	adept_stack.jacobian_reverse(jacobian);
      }
      else {
	adept_stack.jacobian(jacobian);
      }
      timer.stop();
      if (verbose) {
	if (!jac_is_printed) {
	  std::cerr << "Adept Jacobian\n";
	  for (int i = 0; i < NX; i++) {  
	    for (int k = 0; k < NX; k++) {
	      std::cerr << " " << jacobian[i*NX + k];
	    }
	    std::cerr << "\n";
	  }
	}                                     
      }
      jac_is_printed = true; 
    }
  }
  }




 
#ifdef ADOLC
  int adolc_record_id = timer.new_activity("ADOL-C record");
  int adolc_adjoint_id = timer.new_activity("ADOL-C adjoint");
  int adolc_jacobian_id = timer.new_activity("ADOL-C Jacobian");
  if (do_all || do_adolc) {
    bool jac_is_printed = false;
    std::cout << "ADOL-C\n";

    for (int j = 0; j < nr; j++)
      {
	adolc_d adolc_q_init[NX];
	adolc_d adolc_q[NX];
	timer.start(adolc_record_id);
	
	trace_on(1,1);
	for (int i = 0; i < NX; i++) {
	  adolc_q_init[i] <<= q_init[i];
	}
	
	if (is_lax_wendroff) {
	  lax_wendroff(nt, dt, adolc_q_init, adolc_q);
	}
	else {
	  toon(nt, dt, adolc_q_init, adolc_q);
	}
	
	for (int i = 0; i < NX; i++) {
	  adolc_q[i] >>= q[i];
	}
	trace_off();
	
	timer.stop();
	if (j == nr-1) 
	  {
	    size_t counts[11];
	    tapestats(1, counts);
	    std::cout << "Counts: "
		      << counts[0] << " "
		      << counts[1] << " "
		      << counts[2] << " "
		      << counts[3] << " "
		      << counts[4] << " "
		      << counts[5] << "\n";    
	  }
	
	
	if (!do_jacobian) {

	  timer.start(adolc_adjoint_id);
	  reverse(1,NX, NX,0, q_AD, q_init_AD);
	  timer.stop();
      
	  if (verbose) {
	    if (j == nr-1) {
	      std::cout << "Values:";
	      for (int i = 0; i < NX; i++) { std::cout << " " << q[i]; }
	      std::cout << "\nAdjoints";
	      for (int i = 0; i < NX; i++) { std::cout << " " << q_init_AD[i]; }
	      std::cout << "\n";
	    }
	  }
      
	}
	else {
	  static double** jac = myalloc2(NX,NX);
	  timer.start(adolc_jacobian_id);
	  if (force_jacobian < 0) {
	    double** I = myallocI2(NX);
	    double* result = myalloc1(NX);
	    int rc = zos_forward(1, NX, NX, 1, q_init, result);
	    if (rc < 0) {
	      std::cerr << "ERROR OCCURRED IN ADOL-C's zof_forward()\n";
	      exit(rc);
	    }
	    MINDEC(rc,fov_reverse(1, NX, NX, NX, I, jac));
	    myfreeI2(NX, I);
	    myfree1(result);
	  }
	  else if (force_jacobian > 0) {
	    double* result = myalloc1(NX);
	    double** I = myallocI2(NX);
	    int rc = fov_forward(1, NX, NX, NX, q_init, I, result, jac);
	    myfreeI2(NX, I);
	    myfree1(result);
	  }
	  else {
	    jacobian(1, NX, NX, q_init, jac);
	  }
	  timer.stop();
	  if (verbose) {
	    if (!jac_is_printed) {
	      std::cerr << "ADOL-C: Jacobian\n";
	      for (int i = 0; i < NX; i++) {
		for (int k = 0; k < NX; k++) {
		  std::cerr << " " << jac[k][i];
		}    
		std::cerr << "\n";
	      } 
	      jac_is_printed = true;
	    }
	  }                                                                              
	}                 
      }             

  }

#endif

#ifdef CPPAD
  int cppad_record_id = timer.new_activity("CppAD record");
  int cppad_process_id = timer.new_activity("CppAD process");
  int cppad_adjoint_id = timer.new_activity("CppAD adjoint");
  int cppad_jacobian_id = timer.new_activity("CppAD Jacobian");
  if (do_all || do_cppad) {
    bool jac_is_printed = false;
    CppAD::thread_alloc::hold_memory(true);

    std::cout << "CppAD\n";

    for (int j = 0; j < nr; j++)
      {
	std::vector<cppad_d> cppad_q_init(NX);

	for (int i = 0; i < NX; i++) {
	  cppad_q_init[i] = q_init[i];
	}
	CppAD::Independent(cppad_q_init);
	std::vector<cppad_d> cppad_q(NX);
	timer.start(cppad_record_id);
	if (is_lax_wendroff) {
	  lax_wendroff(nt, dt, &cppad_q_init[0], &cppad_q[0]);
	}
	else {
	  toon(nt, dt, &cppad_q_init[0], &cppad_q[0]);
	}
	
	timer.stop();
	if (j == nr-1) {
       // I needed to hack the CppAD library for the following line to work:
	  size_t mem = cppad_q[0].tape_this()->memory();
	  std::cout << "Memory used by tape: " << mem << " bytes\n";
	}
	timer.start(cppad_process_id);
	CppAD::ADFun<double> f(cppad_q_init, cppad_q);
	if (j == nr-1) {
	  std::cout << "Memory used by function: " << f.size_op_seq() << " bytes\n";                   
	}
	
	if (!do_jacobian) {
	  std::vector<double> cppad_q_AD(NX), cppad_q_init_AD(NX);
	  for (int i = 0; i < NX; i++) {
	    cppad_q_AD[i] = q_AD[i];
	  }
	  timer.start(cppad_adjoint_id);
 
	  cppad_q_init_AD = f.Reverse(1, cppad_q_AD);
	  timer.stop();
	  if (verbose) {
	    if (j == nr-1) {
	      std::cout << "Values:";
	      for (int i = 0; i < NX; i++) { std::cout << " " << CppAD::Value(cppad_q[i]); }
	      std::cout << "\nAdjoints";
	      for (int i = 0; i < NX; i++) { std::cout << " " << cppad_q_init_AD[i]; }
	      std::cout << "\n";
	    }
	  }
	}
	else {
	  static std::vector<double> jacobian(NX*NX);
	  static std::vector<double> x(NX);                            
	  for (int i = 0; i < NX; i++) x[i] = q_init[i];

	  timer.start(cppad_jacobian_id);
  
	  if (force_jacobian < 0) {
	    CppAD::JacobianRev(f, x, jacobian); 
	  }
	  else if (force_jacobian > 0) {
	    CppAD::JacobianFor(f, x, jacobian); 
	  }
	  else {
	    jacobian = f.Jacobian(x);
	  }
	  if (!jac_is_printed) {
	    if (verbose) {
	      std::cerr << "CppAD: Jacobian\n";
	      for (int i = 0; i < NX; i++) {                             
		for (int k = 0; k < NX; k++) {                          
		  std::cerr << " " << jacobian[i + k*NX];                
		}                                                         
		std::cerr << "\n";                                        
	      }                                                           
	      jac_is_printed = true;
	    }
	  }                                                             
	}

      }
  }
#endif 

#ifdef SACADO
  int sacado_record_id = timer.new_activity("Sacado record");
  int sacado_adjoint_id = timer.new_activity("Sacado adjoint");
  if (do_all || do_sacado) {
      std::cout << "Sacado\n";
  for (int j = 0; j < nr; j++)
  {
    sacado_d sacado_q_init[NX];
    sacado_d sacado_q[NX];
    timer.start(sacado_record_id);

    for (int i = 0; i < NX; i++) {
      sacado_q_init[i] = q_init[i];
    }
    
    if (is_lax_wendroff) {
      lax_wendroff(nt, dt, sacado_q_init, sacado_q);
    }
    else {
      toon(nt, dt, sacado_q_init, sacado_q);
    }
    timer.start(sacado_adjoint_id);

    sacado_d objective_func = 0.0;
    for (int i = 0; i < NX; i++) {
      objective_func += sacado_q[i] * q_AD[i];
    }

    Sacado::Rad::ADvar<double>::Gradcomp();
    for (int i = 0; i < NX; i++) { q_init_AD[i] = sacado_q_init[i].adj(); }
    timer.stop();

    if (j == nr-1) {
      if (verbose) {
	std::cout << "Values:";
	for (int i = 0; i < NX; i++) { std::cout << " " << sacado_q[i].val(); }
	std::cout << "\nAdjoints";
	for (int i = 0; i < NX; i++) { std::cout << " " << q_init_AD[i]; }
	std::cout << "\n";
      }
      std::cout << "Memory usage: " << Sacado::Rad::ADmemblock<double>::num_blocks()*8 << " kilobytes\n";
    }

  }
  }
#endif

#ifdef SACADO_FAD
  int sacado_fad_jacobian_id = timer.new_activity("Sacado FAD Jacobian");

  if (do_sacado_fad) {
    bool jac_is_printed = false;
  for (int j = 0; j < nr; j++)
  {
    sacado_fad_d sacado_q_init[NX];
    sacado_fad_d sacado_q[NX];

    for (int i = 0; i < NX; i++) {
      sacado_q_init[i] = q_init[i];
      sacado_q_init[i].resize(NX);
      sacado_q[i].resize(NX);
      sacado_q_init[i].fastAccessDx(i) = 1.0;
    }


    timer.start(sacado_fad_jacobian_id);

    if (is_lax_wendroff) {
      lax_wendroff(nt, dt, sacado_q_init, sacado_q);
    }
    else {
      toon(nt, dt, sacado_q_init, sacado_q);
    }
    timer.stop();

    if (!jac_is_printed) {
      if (verbose) {
	std::cerr << "Sacado: Jacobian forward\n";
	for (int i = 0; i < NX; i++) {                               
	  for (int k = 0; k < NX; k++) {                            
	    std::cerr << " " << sacado_q[k].dx(i);
	  }
	  std::cerr << "\n";
	}
	jac_is_printed = true;
      }
    }                                     
  }
  }
#endif


  int original2_id = timer.new_activity("Original function (2)");

  for (int j = 0; j < nr; j++)
  {
    timer.start(original2_id);
    if (is_lax_wendroff) {
      lax_wendroff(nt, dt, q_init, q);
    }
    else {
      toon(nt, dt, q_init, q);
    }
    timer.stop();
    if (j == nr-1) {
      std::cout << "Original function (2)\n";
      if (verbose) {
	std::cout << "Values:";
	for (int i = 0; i < NX; i++) { std::cout << " " << q[i]; }
	std::cout << "\n";
      }
    }
  }


  const std::vector<double>& t = timer.timings();
  double algo_time = t[original2_id];
  std::cout << "------------------------------\n";
  std::cout << "Algorithm base time (warm-up): " << t[original_id] << " s\n";
  std::cout << "Algorithm base time: " << algo_time << " s\n";
  std::cout << "Handcoded relative time: " << t[hand_coded_id]/algo_time << " s\n";
  if (!do_jacobian) {
    std::cout << "Adept: " << (t[adept_record_id]+t[adept_adjoint_id])/algo_time
	      << " (" << t[adept_record_id]/algo_time
	      << ", " << t[adept_adjoint_id]/algo_time << ")\n";
#ifdef SACADO
    std::cout << "Sacado: " << (t[sacado_record_id]+t[sacado_adjoint_id])/algo_time
	      << " (" << t[sacado_record_id]/algo_time
	    << ", " << t[sacado_adjoint_id]/algo_time << ")\n";
#endif
#ifdef ADOLC
    std::cout << "ADOL-C: " << (t[adolc_record_id]+t[adolc_adjoint_id])/algo_time
	      << " (" << t[adolc_record_id]/algo_time
	      << ", " << t[adolc_adjoint_id]/algo_time << ")\n";
#endif
#ifdef CPPAD
    std::cout << "CppAD: " << (t[cppad_record_id]+t[cppad_process_id]+t[cppad_adjoint_id])/algo_time
	      << " (" << t[cppad_record_id]/algo_time
	      << ", " << t[cppad_process_id]/algo_time
	      << ", " << t[cppad_adjoint_id]/algo_time << ")\n";;
#endif
  } 
  else {
    std::cout << "JACOBIANS:\n";
    std::cout << "Adept: " << (t[adept_record_id]+t[adept_jacobian_id])/algo_time
	      << " (" << t[adept_record_id]/algo_time
	      << ", " << t[adept_jacobian_id]/algo_time << ")\n";
#ifdef SACADO_FAD
    std::cout << "Sacado FAD: " << (t[sacado_fad_jacobian_id])/algo_time << "\n";
#endif
#ifdef ADOLC
    std::cout << "ADOL-C: " << (t[adolc_record_id]+t[adolc_jacobian_id])/algo_time
	      << " (" << t[adolc_record_id]/algo_time
	      << ", " << t[adolc_jacobian_id]/algo_time << ")\n";
#endif
#ifdef CPPAD
    std::cout << "CppAD: " << (t[cppad_record_id]+t[cppad_process_id]+t[cppad_jacobian_id])/algo_time
	      << " (" << t[cppad_record_id]/algo_time
	      << ", " << t[cppad_process_id]/algo_time
	      << ", " << t[cppad_jacobian_id]/algo_time << ")\n";;
#endif
  }
  std::cout << "------------------------------\n";


  //  upwind_burgers(200, q_init, q);
}

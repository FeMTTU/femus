#ifndef _femusdefault_h
#define _femusdefault_h


// This file is meant to hold all the DEFAULT configuration of the library.
// The idea is that EACH PARAMETER HER CAN EVENTUALLY BE OVERRIDDEN at RUNTIME,
// either from a file to be read at runtime, or from the main, or from command line, or something

// ALWAYS REMEMBER THAT THE IFDEFS ARE DANGEROUS BECAUSE YOU HAVE TO MAKE SURE 
// THAT THIS HEADER IS INCLUDED IN THE CORRECT PLACES
// FOR EVERY MACRO VARIABLE YOU SHOULD DO A GREP TO FIND ALL THE FILES IN WHICH IT IS NEEDED
// AND THEN INCLUDE THIS HEADER EXPLICITLY, NOT THROUGH INDIRECT INCLUDES!


//*********************************************
//************** FILES ************************
#define DEFAULT_EXT_H5         ".h5"
#define DEFAULT_EXT_XDMF       ".xmf" 
#define DEFAULT_EXT_IN         ".in" 
#define DEFAULT_EXT_LOG        ".log"
#define DEFAULT_AUX_XDMF       "Xdmf.dtd" 
// TODO ".gam" and ".med" for gambit and salome should be known to libmesh but we don't want to put OUR OWN includes in LIBMESH!!! 
// // // // # ----  Mesh class ----
#define DEFAULT_BASEMESH      "mesh" 
#define DEFAULT_MULTIMESH     "multimesh" 
#define DEFAULT_CONNLIN       "_conn_lin" 
// // // // # ----  Multigrid  --------------
#define DEFAULT_F_MATRIX     "Matrix" 
#define DEFAULT_F_PROL       "Prol" 
#define DEFAULT_F_REST       "Rest" 
// // // // # ----  Time, sort of  --------------
#define DEFAULT_BASESOL       "sol" 
#define DEFAULT_BASECASE      "case" 
#define DEFAULT_BASETIME      "time" 
// // // // # ---- RESTART -------- // # if a run reaches the end, then we write it as a "default restart" run
#define DEFAULT_LAST_RUN      "run_to_restart_from"
#define DEFAULT_CASE_DATA     "case.txt"
// // // // # ----  BCHandling  -------------- // # external ibc file
#define DEFAULT_IBC           "case" 
#define DEFAULT_BDRY_SUFFIX   "_bd"
// // // // # ----  log  --------------
#define DEFAULT_RUN_LOG       "run"
 

#define DEFAULT_BASEPATH     "./"
#define DEFAULT_FEMDIR       "fem/"
#define DEFAULT_CONTRIBDIR   "contrib/"
#define DEFAULT_CONFIGDIR    "input/"
#define DEFAULT_CASEDIR      "case/"
#define DEFAULT_OUTPUTDIR    "output/"  //we have to hardcode it here otherwise we cannot perform restart without logical inconsistencies
#define DEFAULT_RUNTIMECONF  "femus_conf.in"
//*********************************************


//*********************************************
//********** PRINT INFO ***********************
// This is for configuring the "verbosity" of FEMuS:
// important for a log of the program.
// Deactivate it if you want to make faster runs.
// // // #define DEFAULT_FEMUS_PRINT_GRAPHIC 0 //this will be turned into runtime
#define DEFAULT_PRINT_INFO 1
#define DEFAULT_PRINT_TIME 1
#define DEFAULT_PRINT_CONV 1

#if DEFAULT_PRINT_TIME==1
#include <ctime>


namespace femus {


#endif
//*********************************************

//*********************************************
//********** BOUNDARY CONDITIONS **************
#define DEFAULT_BC_FLAG 1 //=0 if you put a function on the RHS //you also have to comment the bc_read
//*********************************************

//****************************************
//********** NUMERIC VECTOR **************
// Default tolerance used when comparing two NumericVectors
// Actually in libmesh it is much less! 
#define DEFAULT_NUMVEC_TOLERANCE 1.e-20
//****************************************

//**************************************************************************************
//************************EQNBASE - MULTIGRID ******************************************
//**************************************************************************************
#define DEFAULT_EPS_LSOLV  1.e-6 //1.e-20
#define DEFAULT_EPS_LSOLV_C  1.e-20//1.e-10
#define DEFAULT_EPS_PREPOST 1.e-20
#define DEFAULT_MAXITS_LSOLV  40
// #define DEFAULT_MAXITS_LSOLV_C  40  //it's not there now!!! is it required a special one for the coarsest?
#define DEFAULT_MG_GAMMA 1
#define DEFAULT_NC_PRE    8 //16
#define DEFAULT_NC_COARSE 40
#define DEFAULT_NC_POST   8 //16 
#define DEFAULT_REST_SIMPLE  0   //no simple restrictor
//**************************************************************************************



} //end namespace femus



#endif



//********* CHOICE of RESTRICTOR TYPE *********
//The MG matrices are generated for ONE SCALAR VARIABLE,
//both the QUADRATIC and the LINEAR part,
//and referred to the GLOBAL MATRIX.
//They are available from the SolverBase class but
//they are read in the SolDA class
//I think they depend on the choice of the FE spaces
//for solving the considered equation

/*
 * This file is related to the configuration of the multigrid algorithm
 * Now the idea is the following. We cannot have so many files for configuration,
 * at least not in this manner.
 * If we want to have separate files, we may consider having one file per CLASS.
 * So we should do a class for multigrid management, which shouldnt be a bad idea 
 * actually. Now, the point is: we may want to configure classes at run-time or at compile time.
 * so we need a .in in one case, or a .h .
 * 
 * This MG_conf is basically a configuration of the EqnBASE class for now.
 * 
 * Notice that this conf file is used by both MAIN and GENCASE
 * 
 * Now, the point is: if you have to configure something,
 * you can decide some ways.
 *    At COMPILE TIME, you may either give the #define directives,
 * or you specify the parameters in the MAIN function,
 * so that they are passed to your class or function somehow.
 *    At RUN TIME, you decide a STANDARD way to READ either from file,
 * or from command line. In any case, the way you specify 
 * a parameter must be DOCUMENTED: you must know the exact word
 * for that parameter, for instance "dt" in a file or "--dt" at command line.
 * Well, this also holds for the #define variables.
 */

#ifndef __femus_utils_FemusDefault_hpp__
#define __femus_utils_FemusDefault_hpp__


// This file is meant to hold all the DEFAULT configuration of the library.
// The idea is that EACH PARAMETER HERE CAN EVENTUALLY BE OVERRIDDEN at RUNTIME,
// either from a file to be read at runtime, or from the main, or from command line, or something

#define DEFAULT_SOL_NCHARS   32    //GMV accepts at most 32 characters for a field

//*********************************************
//************** FILES ************************
#define DEFAULT_EXT_H5         ".h5"
#define DEFAULT_EXT_XDMF       ".xmf" 
#define DEFAULT_EXT_IN         ".in" 
#define DEFAULT_EXT_LOG        ".log"
#define DEFAULT_AUX_XDMF       "Xdmf.dtd" 
// // // // # ----  Mesh class ----
#define DEFAULT_BASEMESH_LIN  "mesh_linear"
#define DEFAULT_BASEMESH_BIQ  "_biquadratic"
#define DEFAULT_BASEMESH      "mesh_biquadratic" 
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
// // // // # ----  BCHandling  --------------
#define DEFAULT_BDRY_SUFFIX   "_bd"
// // // // # ----  log  --------------
#define DEFAULT_RUN_LOG       "run"
 

#define DEFAULT_BASEPATH     "./"
#define DEFAULT_INPUTDIR    "input/"
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
#endif
//*********************************************

//*********************************************
//********** BOUNDARY CONDITIONS **************
#define DEFAULT_BC_FLAG 1 //=0 if you put a function on the RHS //you also have to comment the bc_read
#define DEFAULT_BDRY_TOLL  0.0000001
//*********************************************


//**************************************************************************************
//************************EQNBASE - MULTIGRID ******************************************
//**************************************************************************************
#define DEFAULT_EPS_LSOLV  1.e-6 //1.e-20
#define DEFAULT_EPS_LSOLV_C  1.e-20//1.e-10
#define DEFAULT_EPS_PREPOST 1.e-20
#define DEFAULT_MAXITS_LSOLV  40
#define DEFAULT_MG_GAMMA 1
#define DEFAULT_NC_PRE    8 //16
#define DEFAULT_NC_COARSE 40
#define DEFAULT_NC_POST   8 //16 
#define DEFAULT_REST_SIMPLE  0   //no simple restrictor
//**************************************************************************************

//**************************************************************************************
#define DEFAULT_NDIGITS   4  //n of digits for the timestep printing

#endif
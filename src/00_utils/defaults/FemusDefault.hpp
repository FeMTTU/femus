#ifndef __femus_utils_FemusDefault_hpp__
#define __femus_utils_FemusDefault_hpp__


// This file is meant to hold all the DEFAULT configuration of the library.
// The idea is that EACH PARAMETER HERE CAN EVENTUALLY BE OVERRIDDEN at RUNTIME,
// either from a file to be read at runtime, or from the main, or from command line, or something


//*********************************************
//************** FILES - BEGIN ************************
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
#define DEFAULT_PACKAGE_DIR    "femus/"
#define DEFAULT_INPUTDIR    "input/"
#define DEFAULT_OUTPUTDIR    "output/"  //we have to hardcode it here otherwise we cannot perform restart without logical inconsistencies
#define DEFAULT_SOURCE_DIR    "src/"
#define DEFAULT_APPLICATIONS_DIR    "applications/"
#define DEFAULT_MESH_FILES_PATH    "src/06_mesh/00_single_level/01_input/00_mesh_files/"
#define DEFAULT_RUNTIMECONF  "femus_conf.in"
//************** FILES - END ************************
//*********************************************


//*********************** FILE PRINTING - BEGIN ***************************************************************
#define DEFAULT_SOL_NCHARS   32    //GMV accepts at most 32 characters for a field
#define DEFAULT_NDIGITS   4  //n of digits for the timestep printing
//*********************** FILE PRINTING - END ***************************************************************


//*********************************************
//********** PRINT INFO - BEGIN ***********************
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
//********** PRINT INFO - END ***********************
//*********************************************






#endif

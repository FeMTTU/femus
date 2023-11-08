#ifndef __femus_utils_FemusDefault_hpp__
#define __femus_utils_FemusDefault_hpp__


// This file is meant to hold all the DEFAULT configuration of the library.
// The idea is that EACH PARAMETER HERE CAN EVENTUALLY BE OVERRIDDEN at RUNTIME,
// either from a file to be read at runtime, or from the main, or from command line, or something


//*********************************************
//************** FILES - BEGIN ****************
 // # if a run reaches the end, then we write it as a "default restart" run
#define DEFAULT_LAST_RUN      "run_to_restart_from"
 

#define DEFAULT_BASEPATH     "./" //perhaps one day we should consider to treat these slashes appropriately in case we compile on Windows etc.

// # ----  Applications - BEGIN --------------
#define DEFAULT_RUNTIMECONF  "femus_conf.in"
// # ----  Applications - END --------------

// # ----  Library - BEGIN --------------
#define DEFAULT_SOURCE_DIR          "src/"
#define DEFAULT_MESH_FILES_PATH     "src/06_mesh/00_single_level/01_input/00_mesh_files/"
// # ----  Library - END --------------

//************** FILES - END ************************
//*********************************************




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

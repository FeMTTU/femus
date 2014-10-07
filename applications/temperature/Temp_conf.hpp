#ifndef __physuserconf_h__
#define __physuserconf_h__


//********* SPACE DIMENSION ************

#define DIMENSION    2
//   #define DIMENSION    3
// **************************************
  
//TODO it seems like this .h file is NOT in the .depend, how come?
 
//let us give a definition of configuration.
//What you pass to the constructor is configuration?
//What you use to fully define the function is configuration?

//What you pass at runtime is configuration or switch?


//does this go only in the main?



//====== PHYSICAL FRAMEWORKS and EQUATIONS ==============
// This file is only for ACTIVATING the QUANTITIES and EQUATIONS,
//There is no configuration in it


//=============================
//===================================
//===========PHYSICAL FRAMEWORKS FLAGS 
//====================================
//==============================


//**********************************
#define T_EQUATIONS     1


//********************************************//********************************************
//********************************************//********************************************


//********************************************
//******************Quantity configuration
//********************************************
// #define QTYFE_



#define FE_TEMPERATURE  QQ     //quadratic
#define FE_TEMPERATURE2  LL
#define FE_TEMPERATURE3  KK


//**************************************************************************************
//***************************** EQNT *********************************************************
//**************************************************************************************

#define TEMP_DIR_PENALTY 0


// #define SOLVERT VANKATM =======================
#define SOLVERT GMRES



#endif
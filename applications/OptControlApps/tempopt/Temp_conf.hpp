#ifndef __physuserconf_h__
#define __physuserconf_h__


//********* SPACE DIMENSION ************

#define DIMENSION    2
//   #define DIMENSION    3
// **************************************
  
 
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

//*********************************
#define VELOCITY_QTY 1
//Provided quantities: Velocity
//*********************************
#define PRESSURE_QTY 1
//Provided quantities: Pressure

//*********************************
#define TEMP_QTY 1
//Provided quantities: Temperature
//*********************************


//===============================
//===============================
//===========EQUATION FLAGS 
//==========================
//===============================
// Required Quantities, i.e. related Unknowns, i.e. Solves for, i.e. Provides (like Linux Packages!):


//**********************************
#define NS_EQUATIONS    1

#if (NS_EQUATIONS==1)  //Dependency of this Equation on the Quantities it computes for
    #define VELOCITY_QTY 1
    #define PRESSURE_QTY 1
#endif
//**********************************


//**********************************
#define T_EQUATIONS     1


//********************************************//********************************************
//********************************************//********************************************


//********************************************
//******************Quantity configuration
//********************************************
// #define QTYFE_

#define FE_VELOCITY   QQ     //quadratic
#define FE_PRESSURE   LL      // linear

#define FE_TEMPERATURE  QQ     //quadratic

#define  T4_ORD KK

#define FOURTH_ROW 1

//what happens with FE linear magnetic field
//what happens with quadratic pressure....
//what happens with LINEAR velocity...



//********************************************
//******************Equations configuration
//********************************************




//**************************************************************************************
//***************************** EQNNS *********************************************************
//**************************************************************************************

//use or not penalty for Dirichlet
//now penalty Dirichlet works for all domains,
// nodal Dirichlet only for straight domains
#define NS_DIR_PENALTY 0


// #define SOLVERNS VANKANSM
#define SOLVERNS GMRES


//  ===============================
// 3D NAVIER-STOKES ADVECTION term
// ========================================
// A) Stokes flow    ADVPIC_NS 0. ADVNEW_NS=0
// B)  Navier-Stokes ADVPIC_NS 1. ADVNEW_NS={0,1}
#define ADVPIC_NS 1.
//  Navier-Stokes  nonlinear iterations
// A) Picard iteration  ADVPIC_NS=1, ADVNEW_NS=0 
// B) Newton iteration  ADVPIC_NS=1, ADVNEW_NS=1
#define ADVNEW_NS 0.
//  Navier-Stokes  stab  +0.5*(div u,u.v)
#define STAB_NS 0.
#define  KOMP_NS 0. //1.e-6





//Pay attention, from now on the things were REDEFINED!

//**************************************************************************************
//***************************** EQNT *********************************************************
//**************************************************************************************

#define TEMP_DIR_PENALTY 0


// #define SOLVERT VANKATM =======================
#define SOLVERT GMRES



#endif
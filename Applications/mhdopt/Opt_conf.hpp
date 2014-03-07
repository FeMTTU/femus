#ifndef __physuserconf_h__
#define __physuserconf_h__



//********* SPACE DIMENSION, ONLY IN APPLICATION!!! ************
// #define DIMENSION    2
  #define DIMENSION    3

//These  will also be considered in another light
//The library should have a 1D 2D 3D version of everything
//Then here you only say what u want to USE.

// DIMENSION is STILL NEEDED for:
// ASSEMBLE FUNCTION,
// BOUNDARY CONDITIONS,
// JACOBIANS

//****************************************

//********* DOMAIN SHAPE ***************
 //box-shaped domn
//i still have to use DOM_BOX in tune with the MESH FILE
//so the problem is that i tune the mesh file at runtime
//but i have to define the domain at compile time...
//the point is that when you switch the domain
// this switch cannot be done ONLY at RUN-TIME,
// by changing things in the mesh file or in the femuscconf.in file
// some things must be done at compile time also, 
// like the instantiation of the related shape.


//****************************************
//****************************************
//****************************************
//****************************************


//TODO the Quantity #defines should be used to 
//INCLUDE or EXCLUDE a certain PHYSICAL PHENOMENON from the system.
// of course, this is not trivial, because not only you should 
// avoid defining quantities,
// but, if these quantities are unknowns of equations,
// you should avoid defining the equations.autoOr, if these quantities are not unknowns 
// of equations, but they participate to other equations as external quantities,
// you should remove them from the element matrices...
// So, that would mean removing the lines from the loop,
// and you can do that only with a define...
//So, it is not so easy to setup a so-called "PHYSICAL configuration..."


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
#define BMAG_QTY 1
//Provided quantities: MagnFieldHom
//                     MagnFieldHomLagMult
//                     MagnFieldExt
//                     MagnFieldExtLagMult
//*********************************
#define OPT_CONTROL 1
//Provided quantities: DesVelocity
//*********************************
#define TEMP_QTY 1
//Provided quantities: Temperature
//*********************************
// #define PHASE_QTY 1
//Provided quantities: VOFColor
//*********************************
// // // #define TEMP_DEPS 1
//Provided quantities: Density,
//                     Viscosity,
// 		       HeatConductivity,
// 		       SpecificHeatP
//*********************************


//===============================
//===============================
//===========EQUATION FLAGS 
//==========================
//===============================
// Required Quantities, i.e. related Unknowns, i.e. Solves for, i.e. Provides (like Linux Packages!):


//**********************************
 #define NS_EQUATIONS    1
//===== Provides: Velocity, Pressure
//===== C++Name: EqnNS
//===== MapName: "Eqn_NS"
#if (NS_EQUATIONS==1)  //Dependency of this Equation on the Quantities it computes for
    #define VELOCITY_QTY 1
    #define PRESSURE_QTY 1
#endif
//**********************************
 #define     MHD_EQUATIONS 1
//=== Provides: MagnFieldHom
//             (MagnFieldHomLagMult)
//===== C++Name: EqnMHD
//===== MapName: "Eqn_MHD"

//**********************************
   #define MHDCONT_EQUATIONS 1
//=== Provides: MagnFieldExt
//             (MagnFieldExtLagMult)
//===== C++Name: EqnMHDCONT
//===== MapName: "Eqn_MHDCONT"

//**********************************
   #define    NSAD_EQUATIONS 1
//=== Provides: AdjVelocity
//              AdjPressure
//===== C++Name: EqnNSAD
//===== MapName: "Eqn_NSAD"

//**********************************
   #define   MHDAD_EQUATIONS 1
//=== Provides: AdjMagnFieldHom
//             (AdjMagnFieldHomLagMult)
//===== C++Name: EqnMHDAD
//===== MapName: "Eqn_MHDAD"


//*********************************
//   #define PROP_EQUATIONS    1
//==== Provides: Pressure
//===== C++Name: EqnPROP
//===== MapName: "Proj_P"
//***********************************
//  #define PROQ_EQUATIONS    1
//==== Provides: CorrectedVelocity
//===== MapName: "Proj_Q"

//***********************************
// #define TWO_PHASE  1 
//==== Provides: VOFColor

//**********************************
// #define K_TURBULENCE   1
//=== Provides: Kappa

//***********************************
// #define EPS_EQS        1
//=== Provides: Epsilon

//**********************************
// #define OMEGA_EQS      1
//=== Provides: Omega

//**********************************
// #define T_EQUATIONS     1
//=== Provides: Temperature

//**********************************
//  #define TWO_EQUATIONS  1 //two-group neutron diffusion
//=== Provides: NeutrFluxTh,
//              NeutrFluxFast


//********************************************//********************************************
//********************************************//********************************************


//********************************************
//******************Quantity configuration
//********************************************
// #define QTYFE_

#define FE_MAGNFIELDHOM        0   //quadratic
#define FE_MAGNFIELDHOMLAGMULT 1    // linear
#define FE_MAGNFIELDEXT         0   //quadratic
#define FE_MAGNFIELDEXTLAGMULT  1      // linear
#define FE_VELOCITY   0     //quadratic
#define FE_PRESSURE  1      // linear
#define FE_DESVELOCITY  0

#define FE_TEMPERATURE  0    //quadratic

#endif


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
// // //   #define SOLVERNS CGM   //conjugate gradient, ONLY symmetric positive definite


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



//**************************************************************************************
//***************************** EQNMHD *********************************************************
//**************************************************************************************


#define MHD_DIR_PENALTY 0

#define SOLVERMHD GMRES

#define ADV_MHD  1.

#define LAP_MHD 1


//thin wall approximation
#define THIN_WALL 0   //0=no thin wall 
#define CONDRATIO 0.

#define JEXT 0.

 //the contribution of this BDRYelement to the global matrix
//Thin wall approximation in rectangular ducts: B = B_o + b
// b only in the LONGITUDINAL direction, Neumann bc => bc_y=1

//**************************************************************************************
//***************************** EQNNSAD *********************************************************
//**************************************************************************************

#define NSAD_DIR_PENALTY 0


#define SOLVERNSAD GMRES


//**************************************************************************************
//***************************** EQNMHDAD *********************************************************
//**************************************************************************************

#define MHDAD_DIR_PENALTY   0


#define LAPADJ 1.

#define SOLVERMHDAD GMRES


//**************************************************************************************
//***************************** EQNMHDCONT *********************************************************
//**************************************************************************************

#define MHDCONT_DIR_PENALTY 0

#define SOLVERMHDCONT GMRES



//Pay attention, from now on the things were REDEFINED!

//**************************************************************************************
//***************************** EQNPROP *********************************************************
//**************************************************************************************

// // // #define PROP_DIR_PENALTY 1
// // // 
// // // 



//**************************************************************************************
//***************************** EQNPROQ *********************************************************
//**************************************************************************************


// // // 
// // // //#define Newmark
// // // #define first_order
// // // 
// // // //  ===============================
// // // // MULTIGRID PARAMETERS
// // // // ========================================
// // // 
// // // // #define SOLVERNS VANKANSM
// // // #define SOLVERNS GMRES
// // // 
// // // //#define PROJ
// // // 
// // // // #ifdef PROJ
// // // // //the standard incremental pressure-correction scheme
// // // // //#define proj_std_incr
// // // // #endif
// // // 
// // // 
// // // //  ===============================
// // // // 3D NAVIER-STOKES ADVECTION term
// // // // ========================================
// // // // A) Stokes flow    ADVPIC_NS 0. ADVNEW_NS=0
// // // // B)  Navier-Stokes ADVPIC_NS 1. ADVNEW_NS={0,1}
// // // #define ADVPIC_NS 0.
// // // //  Navier-Stokes  nonlinear iterations
// // // // A) Picard iteration  ADVPIC_NS=1, ADVNEW_NS=0 
// // // // B) Newton iteration  ADVPIC_NS=1, ADVNEW_NS=1
// // // #define ADVNEW_NS 0.
// // // 
// // // 
// // // 
// // // //  Navier-Stokes  stab  +0.5*(div u,u.v)
// // // #define STAB_NS 0.
// // // //#define LAMBDA 1500
// // // //#define MU 1000
// // // #define KOMP_NS 1.e-20
// // // #define UP_WIND_NS (1.)
// // // // c -> antisymmetric (0.5) 
// // // #define ADV_ASYM 0.
// // // 
// // // // Crank-Nicolson first order 0. 2nd order 0.5
// // // #define CN_TIME 0.
// // // 
// // // #define SUPG_NS (0.)
// // // //  boundary integral
// // // #define P_0 (0.)
// // // 
// // // // turbulent 
// // // #define MU_LOW (1.e-12)
// // // #define MU_TOP (1.e+12)
// // // 
// // // 
// // // #define SQCMU (0.3)
// // // #define YPLUS (1.)
// // // #define ALPHA0 (1.)


//**************************************************************************************
//***************************** EQNT *********************************************************
//**************************************************************************************

// // // #define FE_T 0
// // // 
// // // // turbulence  ==========================
// // // #define ADVE 1.
// // // 
// // // // Turbulence Prandl number 
// // // #define PRT (0.85)
// // // 
// // // // #define SOLVERT VANKATM =======================
// // // #define SOLVERT GMRES
// // // 
// // // 
// // // // temperature lows ================================
// // //  #define  CONST 1
// // // // #define densityT(x) (1.)
// // // // #define cipiT(x) (1.)
// // // // #define kappa(x)  (1.)


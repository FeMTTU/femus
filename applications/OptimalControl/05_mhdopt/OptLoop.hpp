#ifndef __femus_mhdopt_OptLoop_hpp__
#define __femus_mhdopt_OptLoop_hpp__

#include "Typedefs.hpp"
#include "FemusInputParser.hpp"
#include "SystemTwo.hpp"
#include "TimeLoop.hpp"



namespace femus {


// Forward class
class Files;
class SystemTwo;


class OptLoop  : public TimeLoop {

public:

  OptLoop(Files& files_in, const FemusInputParser<double> & map_in);
 ~OptLoop();

  void optimization_loop(MultiLevelProblem& e_map_in);

  void init_equation_data(const SystemTwo* eqn);

  //====data
  NumericVector * _x_oldopt;  //old optimization step

};



//prototypes that can even stay outside of a class
  double ComputeIntegral (const uint Level, const MultiLevelMeshTwo* mesh, const SystemTwo* eqn);

  int ElFlagControl(const std::vector<double> el_xm, const MultiLevelMesh* mesh);

  void VelDesired(const MultiLevelProblem * ml_prob, CurrentQuantity& myvect, const CurrentElem<double> & currelem, const uint idim);

  double SetInitialCondition(const MultiLevelProblem * ml_prob, const std::vector <double> &xp, const char name[]);

  bool  SetBoundaryCondition(const MultiLevelProblem * ml_prob, const std::vector <double> &xp, const char * name, double &value, const int facename, const double time);



  //********* SPACE DIMENSION, ONLY IN APPLICATION!!! ************
// #define DIMENSION    2
  #define DIMENSION    3

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


//===================================
//===========PHYSICAL FRAMEWORKS FLAGS
//====================================

//**********************************
 #define NS_EQUATIONS    1

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

//**********************************
   #define   MHDAD_EQUATIONS 1
//=== Provides: AdjMagnFieldHom
//             (AdjMagnFieldHomLagMult)
//===== C++Name: EqnMHDAD

//**************************************************************************************
//***************************** EQNMHD *********************************************************
//**************************************************************************************

#define ADV_MHD  1.

#define LAP_MHD 1

//**************************************************************************************
//**************************************************************************************


#define LAPADJ 1.




} //end namespace femus



#endif

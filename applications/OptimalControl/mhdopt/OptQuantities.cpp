
//C++ includes
#include <cmath>
#include <iostream>

//library includes
#include "Typedefs.hpp"
#include "MultiLevelMeshTwo.hpp"

#include "Box.hpp"

//application
#include "OptLoop.hpp"
#include "OptQuantities.hpp"


//=================== BEGIN CONSTRUCTORS ================================
// ==================================================================
// ==================================================================

SpecificHeatP::SpecificHeatP(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in) 
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {
  //let us override any possible mistake in the instantiation
  //the rule is Qty_ + the class name.
  //putting this here i avoid using the default arguments which i dislike because they 
  //must be at the bottom of the list of arguments
  //we shouldnt do like this and give the user the possibility of 
  //defining MORE INSTANTIATIONS of the SpecificHeatP type 
  //TODO but, beware that more instantiations would mean for instance DIFFERENT FUNCTIONS of x,y,z,t,
  //so maybe you'll just want to define another specificHeatP class
  
}

HeatConductivity::HeatConductivity(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in) 
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {
  
}

Viscosity::Viscosity(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in) 
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

}

Density::Density(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in) 
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

}
//======end temp dependence


//===========================================================================
Temperature::Temperature(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in) 
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {

}

//===========================================================================
MagnFieldHom::MagnFieldHom(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

 for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in.GetInputParser()->get("Bref");
}

//===========================================================================
MagnFieldHomAdj::MagnFieldHomAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

  for (uint i=0;i<dim_in;i++) _refvalue[i]=1.;// qtymap_in.ml_prob.GetInputParser().get_par("Bref");
}

//==========================================================================
MagnFieldHomLagMult::MagnFieldHomLagMult(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 
  
  const double Bref = qtymap_in.GetInputParser()->get("Bref");
  const double Uref = qtymap_in.GetInputParser()->get("Uref");
  const double sigmaref = Uref*Bref;
  for (uint i=0;i<dim_in;i++) _refvalue[i]=sigmaref;
  
}

//==========================================================================
MagnFieldHomLagMultAdj::MagnFieldHomLagMultAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 
  
//   const double Bref = qtymap_in.ml_prob.GetInputParser().get_par("Bref");
//   const double Uref = qtymap_in.ml_prob.GetInputParser().get_par("Uref");
//   const double sigmaref = Uref*Bref;
  for (uint i=0;i<dim_in;i++) _refvalue[i]=1.;   //sigmaref;
  
}

//==========================================================================
MagnFieldExt::MagnFieldExt(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {

  for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in.GetInputParser()->get("Bref");
  
}

//===========================================================================
MagnFieldExtLagMult::MagnFieldExtLagMult(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {

  const double Bref = qtymap_in.GetInputParser()->get("Bref");
  const double Uref = qtymap_in.GetInputParser()->get("Uref");
  const double sigmaref = Uref*Bref;
  for (uint i=0;i<dim_in;i++) _refvalue[i]=sigmaref;

}

//===========================================================================
// an equation takes things from the equationsmap
//a quantity takes things from the quantity map

Pressure::Pressure(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

  
  for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in.GetInputParser()->get("pref");
}

//===========================================================================
PressureAdj::PressureAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

   for (uint i=0;i<dim_in;i++) _refvalue[i]= 1./*qtymap_in.ml_prob.GetInputParser()._pref*/;
}


//=========================================================================
Velocity::Velocity(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

    for (uint i=0;i<dim_in;i++) _refvalue[i] =  qtymap_in.GetInputParser()->get("Uref");
  
}

//=========================================================================
VelocityAdj::VelocityAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

    for (uint i=0;i<dim_in;i++) _refvalue[i] = 1.; /*qtymap_in.ml_prob.GetInputParser().get_par("Uref")*/ //TODO
                                                    //do i have to put the same reference value as 
                                                    // the corresponding direct variable?
                                                    //or should this be equal to the reference value for
                                                    //the test function of the corresponding state equation?
  
}

//==========================================================================
DesVelocity::DesVelocity(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 
  
   for (uint i=0;i<dim_in;i++) _refvalue[i] =  qtymap_in.GetInputParser()->get("Uref");

}

//=================== END CONSTRUCTORS ================================
// ==================================================================
// ==================================================================


//=============================================================
///analytical velocity for Hartmann flow
// difference between get_par and optsys:
// in both cases you are "dynamic" somehow

void Velocity::Function_txyz(const double t,const double* xp, double* func) const {

  Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  // we should do this static_cast in the QUANTITY or QUANTITY MAP constructor
  //if there is some domain shape, we see what type it is and we do the static cast
  //if there is no domain shape, we dont need the domain.
  
    //=====ROTATION of the Function
  const double thetaz = box->_domain_rtmap.get("thetaz");
  
  const double rhof   = _qtymap.GetInputParser()->get("rho0");
  const double Uref   = _qtymap.GetInputParser()->get("Uref");
  const double muvel  = _qtymap.GetInputParser()->get("mu0");
  const double MUMHD  = _qtymap.GetInputParser()->get("MUMHD");
  const double SIGMHD = _qtymap.GetInputParser()->get("SIGMHD");
  const double Bref   = _qtymap.GetInputParser()->get("Bref");
  const double Lref   = _qtymap.GetInputParser()->get("Lref");

  const double DpDz   = 1./*0.5*/;  //AAA: change it according to the pressure distribution!!!

  double DpDzad = DpDz*Lref/(rhof*Uref*Uref);

  double Re  = _qtymap.GetInputParser()->get("Re");
  double Rem = _qtymap.GetInputParser()->get("Rem");
  double Hm  = _qtymap.GetInputParser()->get("Hm");
  double S   = _qtymap.GetInputParser()->get("S");


  double Lhalf = 0.5*(box->_le[0] - box->_lb[0]);
  double Lmid  = 0.5*(box->_le[0] + box->_lb[0]);

  double xtr = xp[0] - Lmid/*/Lref*/;

  
 //constant for the real reference length in the Hartmann number
  const double LHm =2.;   //this is because the reference length for Hm is HALF THE WIDTH of the domain, which is Lref=1 now
  const double magnitude = DpDzad*Hm/LHm*(cosh(Hm/LHm) - cosh(Hm/LHm*xtr*Lref/Lhalf)) / (SIGMHD*Bref*Bref*sinh(Hm/LHm)*Uref);
 
  
  func[0] = -sin(thetaz)*magnitude/*/Uref*/;
  func[1] = cos(thetaz)*magnitude;
                                       //add a 4. to the denominator
				       //should check the difference between L and Lref
                                       //TODO check this nondimensionalization
#if (DIMENSION==3)
  func[2] = 0./*/Uref*/;
#endif

  return;

 
}

//============================================================= 
void Velocity::strain_txyz(const double t, const double* xyz,double strain[][DIMENSION]) const {

//here, tau is a tensor, so tau dot n is a vector which in general has a NORMAL and a TANGENTIAL component  
  
    const double Lref = _qtymap.GetInputParser()->get("Lref");
      double ILref = 1./Lref;
      const double lye = _qtymap.GetMeshTwo()->GetDomain()->_domain_rtmap.get("lye");
  const double x=xyz[0];
  const double y=xyz[1];
#if DIMENSION==3
  const double z=xyz[2];
#endif

  
  strain[0][0] = 0.;                     //ux,x
  strain[0][1] = strain[1][0] = 0. ;//0.5*(uy,x+ux,y) 
  strain[1][1] = 0.*(-(lye*ILref-y));                    //uy,y
#if (DIMENSION==3)
  strain[0][2] = strain[2][0] = 0. ;  //0.5*(uz,x+ux,z) 
  strain[1][2] = strain[2][1] = 0. ;  //0.5*(uy,z+uz,y)                                     
  strain[2][2] = 0. ;                    //uz,z
#endif

return;
}


//=============================================================
/// prescribed pressure at the boundary
//no initial condition for pressure is required, because it has no time derivative
//only the boundary condition in the Neumann part of the boundary has to be enforced
//so there is also no problem about consistency between IC and BC values,
// because there are NO IC values for pressure! 

void Pressure::Function_txyz(const double t, const double* xp,double* func) const {
  //the coordinates (x,y,z,t) of the VOLUME domain are NON-DIMENSIONAL
  //and the function value must be nondimensional as well
  //this function receives an ALREADY ROTATED COORDINATE!

  ///this function may receive the values at the DOFS or at the GAUSS POINTS also;
  ///in both cases, the coordinates must be in the REFBOX DOMAIN,
  ///so use xyzrefbox._val_g or xyzrefbox._val_dofs


//here you see that you perform a casting
//of BOTH the DOMAIN and the PHYSICS
//this is clear because the Base functions
//only yield the Base datatypes,
//so if you are in a CHILD CLASS you have to 
//do the CASTS towards the CHILD classes.
//Clearly, the goal would be to do the castings
//NOT INSIDE SMALL ROUTINES but as class members 
//or something, anyway at higher level,
//so that these small functions do not need to do it repeatedly

//the point is that while one of the two casts can be done "a priori"
//because one may establish "a priori" the set of specific Domain Shapes
//Box,Cylinder, whatever  (well, one could actually add it, and update the library...)
//the cast of the OptPhys cannot be done automatically
// so it should be done
//    in ALL the SPECIFIC Quantity constructors
//and in ALL the SPECIFIC Equation constructors
//since all these classes are application - specific,
//you do not perturb the library.

// also, we have to think how to do with the Domains, because we do not want 
// the user to need to update the library for every different domain
// we must think of an Application Specific domain 
// that must be cast after its introduction in the Mesh class.
// Everyone can reach the Domain through the mesh class as a FATHER domain;
//how can we convert it to a specific domain EXPLICITLY
// and STILL STAYING in the Mesh class WITHOUT PERTURBING it?
//the basic classes Mesh, MultiLevelProblemTwo, QuantityMap handle FATHER THINGS.
//the application-specific classes (Equation, Quantity, Physics) handle CHILD things.
//in this passage we must CONVERT at the highest possible level.
//we should INSTANTIATE the CHILDREN in the applications 
// and PASS THEM SEPARATELY to the basic classes and application-specific classes
// no we cant do like this, we pass only once to the basics and then
//convert in the app specific ones.

//yes, for the domain shapes we must find a way to avoid the switch Box Cylinder,
// because if one adds a new shape one SHOULD MODIFY the GENCASE also 
//and we'd want to keep it as small as possible.
// The point is that for now the Gencase CANNOT LIVE without SPECIFIC BOX information
// because we implemented the functions for BOX GENERATION.
//So if one has cylinder one should add the cylinder generation functions
// if one has a Xmas tree one should add the christmas tree generation functions...
// NO, we clearly cannot do that
// i would want the gencase to be AS GENERAL as POSSIBLE

Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  double         lb[DIMENSION];
  double         le[DIMENSION];
  lb[0] = box->_lb[0];//already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
#if DIMENSION == 3
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
#endif
  
 
  func[0] =  1./ _qtymap.GetInputParser()->get("pref")*( (le[1]-lb[1]) - /*x_rotshift*/xp[1] )*(cos(6.*0.*t));

//this equation is in the reference frame CENTERED AT (0,0,0)  
  
  return;
  }

// =================================================
void PressureAdj::Function_txyz(const double t, const double* xp,double* func) const {

  
  func[0] = 0.;
  
  return;
}


// =================================================
void Temperature::Function_txyz(const double t, const double* xp,double* temp) const {
  //the coordinates (x,y,z,t) of the VOLUME domain are NON-DIMENSIONAL
  //and the function value must be nondimensional as well

  //the compidx may be useful when you have to set a component other than zero... for instance,
  //you have three neutron fluxes, or pressure is u_value[4] instead of the first scalar 
  //We pass the pointer and the index, for scalar variables
  //for vector variables we always assume that they go from 0 to DIMENSION-1. Actually,
  //it might happen something different. So,maybe, to make things not very complicated,
  // it suffices to pass the pointer, and then externally one will think of shifting the indices and so on.. anyway.

  const double Tref = _qtymap.GetInputParser()->get("Tref");


  Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());  
  
  temp[0] = 100.*(xp[0])*( ( box->_le[0] - box->_lb[0]) - xp[0])/Tref;
 

  
  return;
  }
  
  
// =================================================
void Temperature::heatflux_txyz(const double t, const double* xyz, double* qflux) const {
  //the coordinates (x,y,z,t) of the VOLUME domain are NON-DIMENSIONAL
  //and the function value must be nondimensional as well
 //-----Nonhomogeneous Neumann-------
 // Qflux = - k grad(T) by definition
//  QfluxDOTn>0: energy flows outside (cooling)  QfluxDOTn<0: energy flows inside (heating)

std::cout << "Temperature: Heatflux, check which coordinates are passed in here" << std::endl;
  Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  const double thetaz = box->_domain_rtmap.get("thetaz");

     qflux[0]=+700.*cos(thetaz);
     qflux[1]=+700.*sin(thetaz);
 #if (DIMENSION==3)
      qflux[2]=0.;
 #endif

  return;
  }




//=============================================================
//this function receives the x values already in the Reference Box
//so no transformation occurs here
void MagnFieldExt::Function_txyz(const double t,const double* xp, double* func) const {

  
  //=====ROTATION of the Function
  Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  const double thetaz = box->_domain_rtmap.get("thetaz");

  //============== PICK THE REQUIRED REFERENCE VALUES for the FUNCTION
  const double Bref   = _qtymap.GetInputParser()->get("Bref");      //Uref*sqrt(rhof*MUMHD);   //in order to make S=1

//function
  func[0] = (cos(thetaz)*Bref
/*             + 0.*x*x
             + 0.*y
             + 0.*(1.-x)
             + 0.*(-2.*x*x +x +1.)
             + 0.*sin(pi*y)*cos(pi*y)*sin(pi*x)*sin(pi*x)
             + 0.*x*x*(x-1.)*(x-1.)*2.*y*(y-1.)*(2.*y-1.)    //div=0
             + 0.*x*(1.-x)*y*(1.-y)*/
             )/Bref; 
  func[1] = (sin(thetaz)*Bref
/*             + 0.*y
             - 0.*sin(pi*x)*cos(pi*x)*sin(pi*y)*sin(pi*y)
             - 0.*y*y*(y-1.)*(y-1.)*2.*x*(x-1.)*(2.*x-1.)    //div=0
             - 0.*x*(1.-x)*y*(1.-y)*/
              )/Bref;
#if (DIMENSION==3)
  func[2] = (0.)/Bref;
#endif
  

return;  
}




//This function receives the ABSOLUTE NON-DIMENSIONAL xp[]
//Then, you convert it
//Well, I'd better do the conversion OUTSIDE, in the DOF PART
//No, suppose you call this function ALONE
//All the things we need for the NONDIMENSIONAL EXPRESSION of the FUNCTION
//are called INSIDE HERE.
//The only things we pass are t,x,y,z of a point

//here, the box measures must be available, because e have to shift once more the reference frame
//now we'll pass it here, but later we'll put that as a class variable or something

void MagnFieldHom::Function_txyz(const double t, const double* xp, double* func) const {

//============== PICK THE REQUIRED REFERENCE VALUES
  const double Lref   = _qtymap.GetInputParser()->get("Lref");
  const double rhof   = _qtymap.GetInputParser()->get("rho0");
  const double Uref   = _qtymap.GetInputParser()->get("Uref");
  const double Bref   = _qtymap.GetInputParser()->get("Bref");      //Uref*sqrt(rhof*MUMHD);   //in order to make S=1

  const double DpDz   = 1.;  //AAA: change it according to the pressure distribution
  // TODO THIS IS DELICATE!!! Suppose you change this multiplicative coefficient, then you get DIFFERENT THINGS!!!!
  //it is just a multiplicative coefficient!
  
  double DpDzad = DpDz*Lref/(rhof*Uref*Uref);

  double Hm  = _qtymap.GetInputParser()->get("Hm");
  double S   = _qtymap.GetInputParser()->get("S");
//=========================================
  
//============= HERE, the analytical solution was given in a reference frame  [-LX,LX] 
//so, we must convert again
//but, NOW le and lb are NONDIMENSIONAL!

 //AAA now I must give a NON-DIMENSIONAL coordinate
//In Paraview I must give a dimensional function in dimensional coordinates, instead

//here the trick is: where you see Hartmann, put the half of it
//where you see S, leave it like that (even if it contains Hm...). This is because S does not contain the reference length, actually! Good.
 
  
    Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  
  double Lhalf = 0.5*(box->_le[0] - box->_lb[0]);
  double Lmid  = 0.5*(box->_le[0] + box->_lb[0]);

  double xtr = /*x_box*/xp[0] - Lmid /*/Lref*/;
  
   const double LHm =2.;  //this is because the reference length for Hm is HALF THE WIDTH of the domain, which is Lref=1 now

   const double thetaz = box->_domain_rtmap.get("thetaz");
 
   const double magnitude = 0.*DpDzad/S*Lhalf/Lref*(sinh(Hm/LHm*Lref/Lhalf*xtr) - xtr*Lref/Lhalf*sinh(Hm/LHm)) / sinh(Hm/LHm);

  func[0] = -sin(thetaz)*magnitude; //0./Bref;
  func[1] = cos(thetaz)*magnitude; //DpDzad/S*Lhalf/Lref*(sinh(Hm/LHm*Lref/Lhalf*xtr) - xtr*Lref/Lhalf*sinh(Hm/LHm)) / sinh(Hm/LHm) ;
#if (DIMENSION==3)
  func[2] = 0./*/Bref*/;
#endif
  
  
  
//  !!!!! REFERENCE TIME!!! if you change Lref, you have to change Uref so as to have TIMEref=1
//AAA p changes!!!
  //Wait a minute: has each equation a different reference time?!? No, the reference time always comes from
  //the Uref and Lref, which are representative of the advection term (u . Nabla)  !!
  
//pay attention! Also here you have to rotate
  // NOT ONLY THE DOMAIN
  // BUT ALSO THE VECTOR QUANTITIES
  //TODO this is DELICATE as well!
  // putting "sin,cos" instead of "-sin,cos" changes the solutions
  // also, it seems like there is a TWO factor, 
  // and the things do not look so perfect.
  //When things look correct but not perfect it is because 
  //the function we are interpolating is not QUADRATIC, but MORE THAN THAT,
  //like a trigonometric function, so it cannot be INTEGRATED perfectly with our quadrature rule
  //in those cases, REFINE THE MESH, which can be done by increasing the NUMBER OF LEVELS
  //where is the TWO FACTOR?
  //well, actually I think that  when the Hartmann number goes big the ratio 
  //between two maximum points is exactly TWO, so actually we have a different Bref... ?
return;  
}






//I'll do that after checking NS + MHD
///Desired velocity for optimal control
void DesVelocity::Function_txyz(const double t, const double* xp,double* func) const {
  
  
  const double Lref = _qtymap.GetInputParser()->get("Lref");
  const double Uref = _qtymap.GetInputParser()->get("Uref");
  double ILref = 1./Lref;
    
  const double rhof   = _qtymap.GetInputParser()->get("rho0");
  const double muvel  = _qtymap.GetInputParser()->get("mu0");
  const double MUMHD  = _qtymap.GetInputParser()->get("MUMHD");
  const double SIGMHD = _qtymap.GetInputParser()->get("SIGMHD");
  const double Bref   = _qtymap.GetInputParser()->get("Bref");

  const double DpDz   = 1./*0.5*/;  //AAA: change it according to the pressure distribution

  double DpDzad = DpDz*Lref/(rhof*Uref*Uref);

  double Re  = _qtymap.GetInputParser()->get("Re");
  double Rem = _qtymap.GetInputParser()->get("Rem");
  double Hm  = _qtymap.GetInputParser()->get("Hm");
  double S   = _qtymap.GetInputParser()->get("S");
 
  
  Box* box= static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());
  
  
  double Lhalf = 0.5*(box->_le[0] - box->_lb[0]);
  double Lmid  = 0.5*(box->_le[0] + box->_lb[0]);

  double xtr = xp[0] - Lmid;

  const double thetaz = box->_domain_rtmap.get("thetaz");

  //constant for the real reference length in the Hartmann number
  const double LHm =2.;   //this is because the reference length for Hm is HALF THE WIDTH of the domain, which is Lref=1 now

  const double magnitude = _qtymap.GetInputParser()->get("udes")*DpDzad*Hm/LHm*(cosh(Hm/LHm) - cosh(Hm/LHm*xtr*Lref/Lhalf)) / (SIGMHD*Bref*Bref*sinh(Hm/LHm)*Uref);
  
  func[0] = -sin(thetaz)*magnitude;
  func[1] = cos(thetaz)*magnitude;
                                       //add a 4 to the denominator
				       //should check the difference between L and Lref
                                       //TODO check this nondimensionalization
				       
//here, I give as target velocity the velocity that would be obtained WITHOUT CONTROL
//therefore, u starts very close to u_d, so I can put a very big alpha
//now, I'll just put get_par("udes") so that I choose to modify the "amplitude"
  
  // get_par("udes")*DpDz*Hm*(cosh(Hm) - cosh(Hm*xtr*Lref/Lhalf)) / (SIGMHD*Bref*Bref*sinh(Hm)*Uref);
//  get_par("udes")/**(x - lxb*ILref)*(lxe*ILref-x)*//Uref;
//  get_par("udes")/Uref;

#if (DIMENSION==3)
  func[2] = 0./*/Uref*/;
#endif

  
  return;

}
 

void VelocityAdj::Function_txyz(const double t, const double* xp,double* func) const{
  
  func[0] = 0./*/Uref*/;
  func[1] = 0./*/Uref*/;
#if (DIMENSION==3)
  func[2] = 0./*/Uref*/;
#endif
  
    return;
  } 

void MagnFieldHomAdj::Function_txyz(const double t, const double* xp,double* func) const{
  
  func[0] = 0./*/Uref*/;
  func[1] = 0./*/Uref*/;
#if (DIMENSION==3)
  func[2] = 0./*/Uref*/;
#endif
  
    return;
  } 



  void MagnFieldExtLagMult::Function_txyz(const double t, const double* xp,double* func) const{
  
  func[0] = 0./*/Uref*/;
  
    return;
  }
  
 void MagnFieldHomLagMult::Function_txyz(const double t, const double* xp,double* func) const{
  
  func[0] = 0./*/Uref*/;
  
    return;
  }
  
   void MagnFieldHomLagMultAdj::Function_txyz(const double t, const double* xp,double* func) const{
  
  func[0] = 0./*/Uref*/;
  
    return;
  }
  
  
  
// ========================================================  
// ========================================================  
// ========================================================  


void Velocity::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);


#if (DIMENSION==2)
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
     bc_flag[0]=0;
     bc_flag[1]=0;  
  }

 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
     bc_flag[0]=0;   //u dot t
//   bc_flag[1]=0;
   }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
    bc_flag[0]=0;
//     bc_flag[1]=0;
  }
  

#elif (DIMENSION==3)

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;    //u dot n
    bc_flag[1]=0;    //u x n
    bc_flag[2]=0;    //u x n 
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll ) {  //right of the RefBox
    bc_flag[0]=0;    //u dot n
    bc_flag[1]=0;   //u x n
    bc_flag[2]=0;   //u x n
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )  {  //bottom  of the RefBox
     bc_flag[0]=0;      //u x n
//      bc_flag[1]=0;   //u dot n   //leave this free for VELOCITY INLET
     bc_flag[2]=0;      //u x n
   }
  
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  {  //top of the  of the RefBox
     bc_flag[0]=0;     //u x n
//      bc_flag[1]=0;  //u dot n   //leave this free for outlet
     bc_flag[2]=0;     //u x n
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
       if (bc_flag[0] == 1)  bc_flag[0] = _qtymap.GetInputParser()->get("Fake3D");   //u x n  //check it for all equations
       if (bc_flag[1] == 1)  bc_flag[1] = _qtymap.GetInputParser()->get("Fake3D");   //u x n          //leave this free for 2D
      bc_flag[2]=0;                                               //u dot n  
  }
  
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  {
     if (bc_flag[0] == 1) bc_flag[0] = _qtymap.GetInputParser()->get("Fake3D");      //u x n
     if (bc_flag[1] == 1) bc_flag[1] = _qtymap.GetInputParser()->get("Fake3D");      //u x n      //leave this free for 2D
     bc_flag[2] = 0;                                                  //u dot n
  }
  
  #endif
  
  
  
  
  return;
 
}



void Pressure::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);


  #if (DIMENSION==2)
  
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
//   bc_flag[0]=0;
  }

 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
//  bc_flag[0]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
    bc_flag[0]=0;
                     //remember that this doesnt mean that the pressure at the boundary is fixed!
		     //its not a Dirichlet boundary condition for pressure!
		     //the initial guess of p at the boundary is changed by the equation after one step
		     //none of the pressure nodes are fixed, they are all COMPUTED after the first step
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
    bc_flag[0]=0;
  }
  

#elif (DIMENSION==3)

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
//      bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll ) {  //right of the RefBox
//      bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )  {  //bottom  of the RefBox
//      bc_flag[0]=0;
  }
  
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  {  //top of the RefBox
      bc_flag[0]=0;     //tau dot n  //PRESSURE OUTLET
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//      bc_flag[0]=0;
  }
  
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  {
//      bc_flag[0]=0;
  }

#endif
  
  return;
 
}



void MagnFieldHom::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);


 #if (DIMENSION==2)

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {
//     bc_flag[0]=0;  //b.n useless with curl curl
     bc_flag[1]=0;  //bxn
  }
  
 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){
//    bc_flag[0]=0;  //b.n useless with curl curl
    bc_flag[1]=0;  //bxn
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  {
     bc_flag[0]=0;  //bxn
//      bc_flag[1]=0;      //leave this free
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
    bc_flag[0]=0;  //bxn
//     bc_flag[1]=0;      //leave this free
 }
  
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) { //  INSULATING
    bc_flag[0]=0;   //b.n
    bc_flag[1]=0;   //bxn
    bc_flag[2]=0;   //bxn
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){ //   INSULATING
    bc_flag[0]=0;   //b.n
    bc_flag[1]=0;   //bxn
    bc_flag[2]=0;   //bxn
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {
    bc_flag[0]=0;       //bxn
//     bc_flag[1]=0;    //b.n     //leave this free for inlet
//     bc_flag[2]=0;    //bxn  //WHY ISNT THIS FIXED AS WELL?
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
    bc_flag[0]=0;         //bxn
//     bc_flag[1]=0;      //b.n    //leave this free for outlet
//     bc_flag[2]=0;     //bxn   //WHY ISNT THIS FIXED AS WELL?
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {   //  CONDUCTING, now insulating
   if (bc_flag[0] == 1) bc_flag[0] = _qtymap.GetInputParser()->get("Fake3D");     //bxn
   if (bc_flag[1] == 1) bc_flag[1] = _qtymap.GetInputParser()->get("Fake3D");     //bxn      //leave this free for 2D
//      bc_flag[2]=0;                                       //b.n 
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {   //  CONDUCTING, now insulating
   if (bc_flag[0] == 1) bc_flag[0] = _qtymap.GetInputParser()->get("Fake3D");   //bxn
   if (bc_flag[1] == 1) bc_flag[1] = _qtymap.GetInputParser()->get("Fake3D");   //bxn     //leave this free for 2D
//  bc_flag[2]=0;                           //b.n 
  }
  
#endif
  
  
  return;
 
}






void MagnFieldHomAdj::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);



 #if (DIMENSION==2)

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
     bc_flag[0]=0;
     bc_flag[1]=0;
  }
  
   if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
    bc_flag[0]=0; 
    bc_flag[1]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox  
     bc_flag[0]=0;
//      bc_flag[1]=0;      //u dot n leave this free
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox  
    bc_flag[0]=0;
//     bc_flag[1]=0;      //u dot n leave this free
  }
  
#else

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
  }
  
   if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
  }
  
     if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox  
     bc_flag[0]=0;
//      bc_flag[1]=0;    //u dot n 
     bc_flag[2]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox  
     bc_flag[0]=0;
//      bc_flag[1]=0;     //u dot n
     bc_flag[2]=0;
  }
  
if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //current  
     if (bc_flag[0] == 1) bc_flag[0] = _qtymap.GetInputParser()->get("Fake3D");
     if (bc_flag[1] == 1) bc_flag[1] = _qtymap.GetInputParser()->get("Fake3D");
//      bc_flag[2]=0;
  }
  
  if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
     if (bc_flag[0] == 1) bc_flag[0] = _qtymap.GetInputParser()->get("Fake3D");
     if (bc_flag[1] == 1) bc_flag[1] = _qtymap.GetInputParser()->get("Fake3D");
//    bc_flag[2]=0;
  }
#endif
 
  
  return;
 
}
 
 
 
 
 
 void MagnFieldHomLagMult::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);


#if (DIMENSION==2)

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {
//   bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){
//   bc_flag[0]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  {
//   bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
//   bc_flag[0]=0;
 }
  
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {
//  bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){
//   bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {
//   bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
//   bc_flag[0]=0;
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//     bc_flag[0]=0;
  }
  
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
//     bc_flag[0]=0;
  }
  
#endif
  
  return;
 
}




void MagnFieldHomLagMultAdj::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);



  #if (DIMENSION==2)

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
//   bc_flag[0]=0;
  }
  
   if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
//   bc_flag[0]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox  
//   bc_flag[0]=0;      //NO BOUNDARY ADJOINT MHD PRESSURE
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox  
//   bc_flag[0]=0;      //NO BOUNDARY ADJOINT MHD PRESSURE 
  }
#else

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
//  bc_flag[0]=0;
  }
  
   if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
//  bc_flag[0]=0;
  }
  
     if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox  
//  bc_flag[0]=0;   //NO BOUNDARY ADJOINT MHD PRESSURE
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox  
//  bc_flag[0]=0;   //NO BOUNDARY ADJOINT MHD PRESSURE
  }
  
if ( x_rotshift[2] > -bdry_toll &&  x_rotshift[2] < bdry_toll ) { //current  
//      bc_flag[0]=0;
  }
  
  if ((le[2]-lb[2]) - x_rotshift[2] > -bdry_toll &&  (le[2]-lb[2]) -x_rotshift[2] < bdry_toll)  {
//      bc_flag[0]=0;
  }
#endif
  
  return;
 
}

 
 
 

void VelocityAdj::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);



#if (DIMENSION==2)
  
 if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
     bc_flag[0]=0;
     bc_flag[1]=0;
  }
  
 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
    bc_flag[0]=0; 
    bc_flag[1]=0;
  }

 if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
     bc_flag[0]=0;
//   bc_flag[1]=0;   //comment it, because you leave u_y free
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
    bc_flag[0]=0;
//  bc_flag[1]=0;  //comment it, because you leave u_y free
  }
  
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
     bc_flag[0]=0;
//      bc_flag[1]=0;    //INSTEAD THIS MUST CORRESPOND TO THE DIRECT
                         //the SPACE of ADJOINT functions is the same as the SPACE for the DIRECT test functions
                         //if you fix this then you dont control...
     bc_flag[2]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
     bc_flag[0]=0;
//      bc_flag[1]=0;     //INSTEAD THIS MUST CORRESPOND TO THE DIRECT
     bc_flag[2]=0;
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
    if (bc_flag[0] == 1) bc_flag[0] = _qtymap.GetInputParser()->get("Fake3D");   //u x n
    if (bc_flag[1] == 1) bc_flag[1] = _qtymap.GetInputParser()->get("Fake3D");  //u x n             //leave this free for 2D
      bc_flag[2]=0;                                            //u dot n
  }
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
    if (bc_flag[0] == 1) bc_flag[0] = _qtymap.GetInputParser()->get("Fake3D");   //u x n
    if (bc_flag[1] == 1) bc_flag[1] = _qtymap.GetInputParser()->get("Fake3D");  //u x n             //leave this free for 2D
    bc_flag[2]=0;                                                //u dot n
  }
#endif
 
  return;
 
}  
 
 
 
 
 void PressureAdj::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);


#if (DIMENSION==2)
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
//   bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){ //right of the RefBox
//  bc_flag[0]=0;
  }

if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
//     bc_flag[0]=0;  //adjoint pressure: computing the integral is not correct, because the function p is prescribed
                    //at the boundary, so deltap = 0 at the boundary.
		    //You dont have to consider the symmetry with the direct equation!
                    //so THERE IS NO BOUNDARY INTEGRAL to be computed
		    //  COMMENT bc FOR THE ADJOINT PRESSURE
                    //The problem is that p_old is modified after one nonlinear step,
		    //so p was changing every time!
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
//    bc_flag[0]=0;   //  COMMENT bc FOR THE ADJOINT PRESSURE
  }
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
//  bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
//  bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
//   bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
//   bc_flag[0]=0;
  }
  
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//   bc_flag[0]=0;
  }
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
//   bc_flag[0]=0;
  }
#endif
  
  
  return;
 
}  



void MagnFieldExt::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);


#if (DIMENSION==2)
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) { //left
     bc_flag[0]=0;
     bc_flag[1]=0;
  }
  
  if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){  //right
    bc_flag[0]=0; 
    bc_flag[1]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
     bc_flag[0]=0;
     bc_flag[1]=0;  //here when you control you must fix it,even if the corresponding b is let free, otherwise the control may use this as well
                      //well,wait,it depends: if you are using the Laplacian for MHD state, then Becont bc's must be consistent with the Laplacian
		      // given there.
		      //if you are using curlxcurl, then Becont should be consistent with curl-curl
		      //on the other hand, Be has a gamma*Laplacian FOR HERSELF... !
		      //so the Becont BC's should be consistent BOTH WITH MHD EQUATION and with MHDCONT EQUATION!
		      
		      
// wait, if you fix all dirichlet, then things go into Becontp pressure
//instead if you fix one Dir and one Neu then nothing goes into pressure

  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {//top of the  of the RefBox
    
     bc_flag[0] = _qtymap.GetInputParser()->get("UseControl");
     bc_flag[1] = _qtymap.GetInputParser()->get("UseControl");
 }
  
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
    bc_flag[0] = _qtymap.GetInputParser()->get("UseControl");
    bc_flag[1] = _qtymap.GetInputParser()->get("UseControl");
    
    bc_flag[2]=0;
  }
   if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
  }
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
    bc_flag[0]=0;
    bc_flag[1]=0;
    bc_flag[2]=0;
  }
  
#endif //DIMENSION

  
  
  return;
 
} 



void MagnFieldExtLagMult::bc_flag_txyz(const double t, const double* xp, std::vector<int> & bc_flag) const  {
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
  Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);

  
  
#if (DIMENSION==2)
  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) { //left
//   bc_flag[0]=0;
  }
  
  if ( (le[0]-lb[0])  -(x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll){  //right
//  bc_flag[0]=0;
  }
  
   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
//   bc_flag[0]=0;     //comment it, don't do the integral
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {//top of the  of the RefBox
//   bc_flag[0]=0;       //comment it, don't do the integral

 }
  
#else

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll ) {  //left of the RefBox
//     bc_flag[0]=0;
  }
  
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll){  //right of the RefBox
//     bc_flag[0]=0;
  }
  
   if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)  {  //bottom  of the RefBox
//     bc_flag[0]=0;
  }
  
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
//     bc_flag[0]=0;
    
  }
   if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll ) {
//     bc_flag[0]=0;
  }
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  {
//     bc_flag[0]=0;
  }
  
#endif //DIMENSION
  
  
  return;
 
}

// =====================================================================
// ===================== INITIAL CONDITIONS ============================
// =====================================================================

void Velocity::initialize_xyz(const double* xp, std::vector< double >& value) const {

  const double Uref = _qtymap.GetInputParser()->get("Uref");
  const double pref = _qtymap.GetInputParser()->get("pref");
  const double udes = _qtymap.GetInputParser()->get("udes");
  
  const double bdry_toll = _qtymap.GetMeshTwo()->GetRuntimeMap().get("bdry_toll");
  
  Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);

//rotation of the function  
    double thetaz = box->_domain_rtmap.get("thetaz");

  
#if (DIMENSION==2)

const double magnitude = /*udes**/1.*(x_rotshift[0] - box->_lb[0])*(box->_le[0]-x_rotshift[0] )/( Uref); 
    value[0] = -sin(thetaz)*magnitude;
    value[1] = cos(thetaz)*magnitude; 

#elif (DIMENSION==3)

    value[0] = 0.;
    value[1] = 0.*/*udes**/(x_rotshift[0] - box->_lb[0])*(box->_le[0]-x_rotshift[0] )/( Uref);
    value[2] = 0.;

   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  {  value[1] = (x_rotshift[0] - lb[0])*(le[0] - x_rotshift[0])*(x_rotshift[2] - lb[2])*(le[2]-x_rotshift[2])/Uref;  }
  
#endif
  
  return;
}


void Pressure::initialize_xyz(const double* xp, std::vector< double >& value) const {

  Box* box = static_cast<Box*>(_qtymap.GetMeshTwo()->GetDomain());

  std::vector<double> lb(_qtymap.GetMeshTwo()->get_dim());
  std::vector<double> le(_qtymap.GetMeshTwo()->get_dim());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (_qtymap.GetMeshTwo()->get_dim() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }
  
  std::vector<double> x_rotshift(_qtymap.GetMeshTwo()->get_dim());
  _qtymap.GetMeshTwo()->_domain->TransformPointToRef(xp,&x_rotshift[0]);
  
#if (DIMENSION==2)
  
      Function_txyz(0.,&x_rotshift[0],&value[0]);
      
#elif (DIMENSION==3)
      
      value[0] = 0.;
  
#endif
  
  return;
}

void MagnFieldHom::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  const double Bref = _qtymap.GetInputParser()->get("Bref");
 
  value[0] = 0./Bref;
  value[1] = 0./Bref;
#if (DIMENSION==3)
  value[2] = 0./Bref;
#endif

  return;
}

void MagnFieldHomLagMult::initialize_xyz(const double* xp, std::vector< double >& value) const {

  const double Uref = _qtymap.GetInputParser()->get("Uref");
  const double Bref = _qtymap.GetInputParser()->get("Bref");
 
  value[0] = 0./(Uref*Bref);
  
  return;
}

void VelocityAdj::initialize_xyz(const double* xp, std::vector< double >& value) const {

  value[0] = 0.;
  value[1] = 0.;
#if (DIMENSION==3)
  value[2] = 0.;
#endif
  
  return;
}

void PressureAdj::initialize_xyz(const double* xp, std::vector< double >& value) const {

  value[0] = 0.;
  
  return;
}

void MagnFieldHomAdj::initialize_xyz(const double* xp, std::vector< double >& value) const {
  
  value[0] = 0.;
  value[1] = 0.;
#if (DIMENSION==3)
  value[2] = 0.;
#endif
  
  return;
}

void MagnFieldHomLagMultAdj::initialize_xyz(const double* xp, std::vector< double >& value) const {

  value[0] = 0.;

  return;
}

void MagnFieldExt::initialize_xyz(const double* xp, std::vector< double >& value) const {

  const double Bref = _qtymap.GetInputParser()->get("Bref");
 
  value[0] = Bref/Bref;
  value[1] = 0./Bref;
#if (DIMENSION==3)
  value[2] = 0./Bref;
#endif
  

  return;
}

void MagnFieldExtLagMult::initialize_xyz(const double* xp, std::vector< double >& value) const {

  const double Uref = _qtymap.GetInputParser()->get("Uref");
  const double Bref = _qtymap.GetInputParser()->get("Bref");
 
  value[0] = 0./(Uref*Bref);

  return;
}


//C++ includes
#include <cmath>
#include <iostream>

//library includes
#include "Typedefs.hpp"
#include "GeomEl.hpp"
#include "MultiLevelMeshTwo.hpp"

#include "Box.hpp"

//application
#include "Opt_conf.hpp"
#include "OptQuantities.hpp"


namespace femus {
  
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

 for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in._physmap->get("Bref");
}

//===========================================================================
MagnFieldHomAdj::MagnFieldHomAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

  for (uint i=0;i<dim_in;i++) _refvalue[i]=1.;// qtymap_in._phys.get_par("Bref");
}

//==========================================================================
MagnFieldHomLagMult::MagnFieldHomLagMult(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 
  
  const double Bref = qtymap_in._physmap->get("Bref");
  const double Uref = qtymap_in._physmap->get("Uref");
  const double sigmaref = Uref*Bref;
  for (uint i=0;i<dim_in;i++) _refvalue[i]=sigmaref;
  
}

//==========================================================================
MagnFieldHomLagMultAdj::MagnFieldHomLagMultAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 
  
//   const double Bref = qtymap_in._phys.get_par("Bref");
//   const double Uref = qtymap_in._phys.get_par("Uref");
//   const double sigmaref = Uref*Bref;
  for (uint i=0;i<dim_in;i++) _refvalue[i]=1.;   //sigmaref;
  
}

//==========================================================================
MagnFieldExt::MagnFieldExt(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {

  for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in._physmap->get("Bref");
  
}

//===========================================================================
MagnFieldExtLagMult::MagnFieldExtLagMult(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {

  const double Bref = qtymap_in._physmap->get("Bref");
  const double Uref = qtymap_in._physmap->get("Uref");
  const double sigmaref = Uref*Bref;
  for (uint i=0;i<dim_in;i++) _refvalue[i]=sigmaref;

}

//===========================================================================
// an equation takes things from the equationsmap
//a quantity takes things from the quantity map

Pressure::Pressure(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

  
  for (uint i=0;i<dim_in;i++) _refvalue[i]= qtymap_in._physmap->get("pref");
}

//===========================================================================
PressureAdj::PressureAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 

   for (uint i=0;i<dim_in;i++) _refvalue[i]= 1./*qtymap_in._phys._pref*/;
}


//=========================================================================
Velocity::Velocity(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

    for (uint i=0;i<dim_in;i++) _refvalue[i] =  qtymap_in._physmap->get("Uref");
  
}

//=========================================================================
VelocityAdj::VelocityAdj(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) {  

    for (uint i=0;i<dim_in;i++) _refvalue[i] = 1.; /*qtymap_in._phys.get_par("Uref")*/ //TODO
                                                    //do i have to put the same reference value as 
                                                    // the corresponding direct variable?
                                                    //or should this be equal to the reference value for
                                                    //the test function of the corresponding state equation?
  
}

//==========================================================================
DesVelocity::DesVelocity(std::string name_in, QuantityMap& qtymap_in, uint dim_in, uint FEord_in)
: Quantity(name_in,qtymap_in,dim_in,FEord_in) { 
  
   for (uint i=0;i<dim_in;i++) _refvalue[i] =  qtymap_in._physmap->get("Uref");

}

//=================== END CONSTRUCTORS ================================
// ==================================================================
// ==================================================================


//=============================================================
///analytical velocity for Hartmann flow
// difference between get_par and optsys:
// in both cases you are "dynamic" somehow

void Velocity::Function_txyz(const double t,const double* xp, double* func) const {

  Box* box = static_cast<Box*>(_qtymap._mesh.GetDomain());
  // we should do this static_cast in the QUANTITY or QUANTITY MAP constructor
  //if there is some domain shape, we see what type it is and we do the static cast
  //if there is no domain shape, we dont need the domain.
  
    //=====ROTATION of the Function
  const double thetaz = box->_domain_rtmap.get("thetaz");
  
  const double rhof   = _qtymap._physmap->get("rho0");
  const double Uref   = _qtymap._physmap->get("Uref");
  const double muvel  = _qtymap._physmap->get("mu0");
  const double MUMHD  = _qtymap._physmap->get("MUMHD");
  const double SIGMHD = _qtymap._physmap->get("SIGMHD");
  const double Bref   = _qtymap._physmap->get("Bref");
  const double Lref   = _qtymap._physmap->get("Lref");

  const double DpDz   = 1./*0.5*/;  //AAA: change it according to the pressure distribution!!!

  double DpDzad = DpDz*Lref/(rhof*Uref*Uref);

  double Re  = _qtymap._physmap->get("Re");
  double Rem = _qtymap._physmap->get("Rem");
  double Hm  = _qtymap._physmap->get("Hm");
  double S   = _qtymap._physmap->get("S");


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
  
    const double Lref = _qtymap._physmap->get("Lref");
      double ILref = 1./Lref;
      const double lye = _qtymap._mesh.GetDomain()->_domain_rtmap.get("lye");
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

Box* box= static_cast<Box*>(_qtymap._mesh.GetDomain());

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
  
 
  func[0] =  1./ _qtymap._physmap->get("pref")*( (le[1]-lb[1]) - /*x_rotshift*/xp[1] )*(cos(6.*0.*t));

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

  const double Tref = _qtymap._physmap->get("Tref");


  Box* box = static_cast<Box*>(_qtymap._mesh.GetDomain());  
  
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
  Box* box = static_cast<Box*>(_qtymap._mesh.GetDomain());
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
  Box* box = static_cast<Box*>(_qtymap._mesh.GetDomain());
  const double thetaz = box->_domain_rtmap.get("thetaz");

  //============== PICK THE REQUIRED REFERENCE VALUES for the FUNCTION
  const double Bref   = _qtymap._physmap->get("Bref");      //Uref*sqrt(rhof*MUMHD);   //in order to make S=1

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
  const double Lref   = _qtymap._physmap->get("Lref");
  const double rhof   = _qtymap._physmap->get("rho0");
  const double Uref   = _qtymap._physmap->get("Uref");
  const double Bref   = _qtymap._physmap->get("Bref");      //Uref*sqrt(rhof*MUMHD);   //in order to make S=1

  const double DpDz   = 1.;  //AAA: change it according to the pressure distribution
  // TODO THIS IS DELICATE!!! Suppose you change this multiplicative coefficient, then you get DIFFERENT THINGS!!!!
  //it is just a multiplicative coefficient!
  
  double DpDzad = DpDz*Lref/(rhof*Uref*Uref);

  double Hm  = _qtymap._physmap->get("Hm");
  double S   = _qtymap._physmap->get("S");
//=========================================
  
//============= HERE, the analytical solution was given in a reference frame  [-LX,LX] 
//so, we must convert again
//but, NOW le and lb are NONDIMENSIONAL!

 //AAA now I must give a NON-DIMENSIONAL coordinate
//In Paraview I must give a dimensional function in dimensional coordinates, instead

//here the trick is: where you see Hartmann, put the half of it
//where you see S, leave it like that (even if it contains Hm...). This is because S does not contain the reference length, actually! Good.
 
  
    Box* box= static_cast<Box*>(_qtymap._mesh.GetDomain());
  
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
  
  
  const double Lref = _qtymap._physmap->get("Lref");
  const double Uref = _qtymap._physmap->get("Uref");
  double ILref = 1./Lref;
    
  const double rhof   = _qtymap._physmap->get("rho0");
  const double muvel  = _qtymap._physmap->get("mu0");
  const double MUMHD  = _qtymap._physmap->get("MUMHD");
  const double SIGMHD = _qtymap._physmap->get("SIGMHD");
  const double Bref   = _qtymap._physmap->get("Bref");

  const double DpDz   = 1./*0.5*/;  //AAA: change it according to the pressure distribution

  double DpDzad = DpDz*Lref/(rhof*Uref*Uref);

  double Re  = _qtymap._physmap->get("Re");
  double Rem = _qtymap._physmap->get("Rem");
  double Hm  = _qtymap._physmap->get("Hm");
  double S   = _qtymap._physmap->get("S");
 
  
  Box* box= static_cast<Box*>(_qtymap._mesh.GetDomain());
  
  
  double Lhalf = 0.5*(box->_le[0] - box->_lb[0]);
  double Lmid  = 0.5*(box->_le[0] + box->_lb[0]);

  double xtr = xp[0] - Lmid;

  const double thetaz = box->_domain_rtmap.get("thetaz");

  //constant for the real reference length in the Hartmann number
  const double LHm =2.;   //this is because the reference length for Hm is HALF THE WIDTH of the domain, which is Lref=1 now

  const double magnitude = _qtymap._physmap->get("udes")*DpDzad*Hm/LHm*(cosh(Hm/LHm) - cosh(Hm/LHm*xtr*Lref/Lhalf)) / (SIGMHD*Bref*Bref*sinh(Hm/LHm)*Uref);
  
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
  


} //end namespace femus


  
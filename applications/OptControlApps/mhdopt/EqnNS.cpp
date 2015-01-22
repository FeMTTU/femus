#include "EqnNS.hpp"

#include "FemusDefault.hpp"

#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "DenseVector.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"

#include "Math.hpp"
#include "MeshTwo.hpp"
#include "GeomEl.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "FETypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "VBTypeEnum.hpp"
#include "QTYnumEnum.hpp"
#include "Domain.hpp"
#include "TimeLoop.hpp"
#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"

//application
#include "Opt_conf.hpp"
#include "EqnMHD.hpp"
#include "EqnMHDCONT.hpp"
#include "OptQuantities.hpp"


namespace femus {


///=============== Constructor
  EqnNS::EqnNS(    std::vector<Quantity*> int_map_in,  //no reference!
	           MultiLevelProblemTwo& equations_map_in,
                   std::string eqname_in,
                   std::string varname_in):
           SystemTwo(int_map_in,equations_map_in,eqname_in,varname_in),
     _AdvPic_fl(ADVPIC_NS),
     _AdvNew_fl(ADVNEW_NS),
     _Stab_fl(STAB_NS),
     _Komp_fac(KOMP_NS)   {

//=======  _var_names[]  ===========
    _var_names[0]="ux"; //variable names
    _var_names[1]="uy";
    _var_names[DIMENSION]="up";
#if (DIMENSION==3)
    _var_names[2]="uz";
#endif
    
//=======  _refvalue[] ==============   
     _refvalue[0] =  _QtyInternalVector[0]->_refvalue[0]; 
     _refvalue[1] =  _QtyInternalVector[0]->_refvalue[1]; 
#if (DIMENSION==3)
     _refvalue[2] =  _QtyInternalVector[0]->_refvalue[2]; 
#endif
     _refvalue[DIMENSION] = _QtyInternalVector[1]->_refvalue[0];

//========= MG solver ===================
  for(uint l=0;l<_NoLevels;l++)  _solver[l]->set_solver_type(SOLVERNS);

//============= DIR PENALTY===============
   _Dir_pen_fl = NS_DIR_PENALTY;
   
    }
//====== END CONSTRUCTOR    
    
    
//================ DESTRUCTOR    
      EqnNS::~EqnNS() {    }
    
    
  void EqnNS::GenMatRhs(const uint Level)  {

#if TEMP_DEPS==1
//for this one i decide not to use any Vect's
//i do not need the vects so badly if i do not have to do space interpolation
//here I need the pointer but not as a Quantity, but as the child class 
//also, since these quantities are not used for interpolation 
//i do not add them to the EXTERNAL MAP. The purpose of that map would be to
//loop automatically over it for getting the dofs.
//but actually it is not so necessary to have it. For instance you may need 
//to use a specific function, in which case you first should do the static cast
//for some Vect and nothing for others
//so, it could be better to do a Vect_LOCAL_EXTERNAL MAP inside here
Density* density_ptr     = static_cast<Density*>(_eqnmap._qtymap.get_qty("Qty_Density"));
Viscosity* viscosity_ptr = static_cast<Viscosity*>(_eqnmap._qtymap.get_qty("Qty_Viscosity"));
#endif  //temp deps
  //====== reference values ========================
//====== related to Quantities on which Operators act, and to the choice of the "LEADING" EQUATION Operator
  //====== Physics
  const double IRe = 1./_phys.get("Re");
  const double IFr = 1./_phys.get("Fr");
  const double   S = _phys.get("S");
//================================================  

//=============== electric current ===============
  //this is a particular Vect. It has already been defined in 3D
  //because we only need it for a CROSS product
  double Jext_g3D[3]={0.,0.,0.}; //Quantity
  //=======density and viscosity===================
  const double rhof = _phys.get("rho0");
  const double  muf = _phys.get("mu0");

//================================================  
//=======Operators @ gauss =======================
  double    dphijdx_g[DIMENSION]; // ShapeDer(): used for Laplacian,Divergence, ..., associated to an Unknown Quantity
  double    dphiidx_g[DIMENSION]; // Test(): used for Laplacian,Advection,Divergence... associated to an Unknown Quantity
  double     AdvRhs_g[DIMENSION]; //Operator: Adv(u,u,phi)
  double  curlBXB_g3D[3];         //Operator  CrossProd(curlB,B)
  double   JextXB_g3D[3];         //Operator: CrossProd(Jext,B)
//================================================  


   const double time =  0.;   //_eqnmap._timeloop._curr_time;
   
  
//========== PROCESSOR INDEX
//every routine we use here should depend directly on this one and not implicitly 
//through the class _iproc. This should be a sort of "function argument",
//like the Level
  const uint myproc= _iproc;

//==========FLAG FOR STATIONARITY OR NOT
//FLAG for the TIME DISCRETIZATION
//every Equation may have a TimeDiscretization
//If different equations have the same TimeDiscretization, 
//we can put them in the MultiLevelProblemTwo instead of the Equation
//TimeDiscretization can be considered as an OPERATOR, like LAPLACIAN or whatever
//What is an Operator characterized by?
//-- by the Quantities it acts on
//-- if it is bilinear or trilinear, two or three Quantities
//-- by a multiplicative coefficient in front of it, which is the REFERENCE VALUE
//-- well, actually the reference value is directly determined by the 
//-- QUANTITIES of the OPERATOR + the "LEADING OPERATOR" of the EQUATION
//-- Every Equation is a BUNCH of OPERATORS, whatever associated
//-- then, the Operator affects the way Ke and Fe are filled
//--- one thing to be careful about Operator is:
//--- do we have to define the list of Quantities the Operator acts on,
//--- or the list of Vects the Operator acts on?
//--- If we do the Quantities we are more general.
//--- on the other hand, if we act on the Vects, we may also have
//--- Operators WITHOUT needing the definition of the Quantities.
//--- for instance g DOT phi, where g is a Vect and phi is the test function of Velocity

//Concerning time, i think we have to consider it as a Quantity, and also space.
//space for sure.
//The different thing about time is that it doesnt have the FE discretization,
//because we use finite difference for it... unless we decide to do 
//finite elements for time as well...
//well, we might do Time as a Vect
const int NonStatNS = (int) _phys.get("NonStatNS");
  const double   dt = 1.; //_eqnmap._timeloop._timemap.get("dt");

//==========FLAG for NONLINEARITY, we might put here
//this flag in fact may involve BOTH MATRIX and RHS, so here we are on top of both
//this is another flag of an OPERATOR, a NONLINEAR Operator
//in case of nonlinear operator, you also have to specify the LINEARIZATION SCHEME.
//So, every Operator then will DROP in the Ke and Fe blocks
//The Choice of the desired nonlinearity should be similar to the choice 
//of the time scheme.
//basically, the point is deciding WHERE to put some operators, either in Ke or Fe,
//(explicit, implicit, semi-implicit, whatever...)  
//and how to loop on the equations and update the time/nonlinear iterations
//So, it may consist either of modifying the matrix in MatRhs, or modifying
//the TimeStep loop or the NonLinear loop
//So nonstationarity or nonlinearity may be characterized not only at the AssemblyMatRhs level,
//but also externally at the loop level.
//So they are a feature of an SystemTwo. Nevertheless they should carry their own stuff,
//for instance the vectors they need.
//for instance: _res,_x,_A, _Prl, _Rst are needed by the MULTIGRID LINEAR SOLVER
//So we should do a class that gathers those vectors which are only related to MG linear solver.
//Then, x_old or x_oold are related to the TimeStep, so we should make a class for them.
//if x_old or x_oold are then also used for NonLinear Steps, then it is an abuse
//to use them also for nonlinear steps, but it works ok.
  
  
  
//========== GEOMETRIC ELEMENT ========
  const uint           space_dim = _mesh.get_dim();
  const uint mesh_ord = (int) _mesh.GetRuntimeMap().get("mesh_ord");
  const uint meshql   = (int) _mesh.GetRuntimeMap().get("meshql"); //======== ELEMENT MAPPING =======
  
//========= BCHandling =========
  const double penalty_val  = _mesh.GetRuntimeMap().get("penalty_val");    

  
   {//BEGIN VOLUME
//===============================  
//===============================  
  
    const uint mesh_vb = VV;
  
    CurrentElem       currelem(VV,this,_mesh,_eqnmap._elem_type);    
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,_eqnmap, _mesh.get_dim());
  
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity VelOld(currgp);
    VelOld._qtyptr   = _QtyInternalVector[QTYZERO]; //an alternative cannot exist, because it is an Unknown of This Equation
    VelOld.VectWithQtyFillBasic();   //the internal quantities will eventually have *this as eqn pointer
    VelOld.Allocate();

   const uint   qtyzero_ord  = VelOld._FEord;
   const uint   qtyzero_ndof = VelOld._ndof; 
    Velocity*  vel_castqtyptr = static_cast<Velocity*>(VelOld._qtyptr); //casting for quantity-specific functions

//=========
    CurrentQuantity pressOld(currgp);
    pressOld._qtyptr   = _QtyInternalVector[QTYONE];
    pressOld.VectWithQtyFillBasic();
    pressOld.Allocate();

   const uint qtyone_ord  = pressOld._FEord;
   const uint qtyone_ndof = pressOld._ndof; 

   //order
   const uint  qtyZeroToOne_DofOffset = VelOld._ndof*VelOld._dim;
   
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();

    //==================Quadratic domain, auxiliary, must be QUADRATIC!!! ==========
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = _mesh.GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
  xyz_refbox.Allocate();
    
//============================ MAG WORLD =======================================
 #if BMAG_QTY==1  
    CurrentQuantity Bhom(currgp); //only to retrieve the dofs
    Bhom._qtyptr   = _eqnmap._qtymap.get_qty("Qty_MagnFieldHom");
    Bhom.VectWithQtyFillBasic();
    Bhom.Allocate();
 
//=========
    CurrentQuantity Bext(currgp);   //only to retrieve the dofs
    Bext._qtyptr   =  _eqnmap._qtymap.get_qty("Qty_MagnFieldExt");
    Bext.VectWithQtyFillBasic();
    Bext.Allocate();

//========= auxiliary, must be AFTER Bhom!
    CurrentQuantity Bmag(currgp); //total
    Bmag._dim        = Bhom._dim;
    Bmag._FEord      = Bhom._FEord;
    Bmag._ndof       = _eqnmap._elem_type[currelem.GetDim()-1][Bmag._FEord]->GetNDofs();
    Bmag.Allocate();
#endif
//======================== MAG WORLD ================================

//===================TEMPERATURE WORLD=============================
#if TEMP_QTY==1
    CurrentQuantity Temp(currgp);
    Temp._qtyptr   =  _eqnmap._qtymap.get_qty("Qty_Temperature");
    Temp.VectWithQtyFillBasic();
    Temp.Allocate();
#endif
//=================== TEMPERATURE WORLD============================

//=======gravity==================================
  CurrentQuantity gravity(currgp);
  gravity._dim=DIMENSION;
  gravity._val_g.resize(gravity._dim);
  gravity._val_g[0] = _phys.get("dirgx");
  gravity._val_g[1] = _phys.get("dirgy");
#if DIMENSION==3
  gravity._val_g[2] = _phys.get("dirgz");
#endif  
 
    const uint nel_e = _mesh._off_el[mesh_vb][_NoLevels*myproc+Level+1];
    const uint nel_b = _mesh._off_el[mesh_vb][_NoLevels*myproc+Level];

//=======================
//=======================    

    for (uint iel=0; iel < (nel_e - nel_b); iel++) {
 
    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.set_el_nod_conn_lev_subd(Level,myproc,iel);
    currelem.SetMidpoint();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);
    _mesh.TransformElemNodesToRef(currelem.GetDim(),currelem.GetNodeCoords(),&xyz_refbox._val_dofs[0]);    

//=======RETRIEVE the DOFS of the UNKNOWN QUANTITIES,i.e. MY EQUATION
    currelem.SetElDofsBc(Level);
    
      VelOld.GetElemDofs(Level);
    pressOld.GetElemDofs(Level);

    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),qtyzero_ord,currelem.GetBCDofFlag()); //only the Qtyzero Part is modified!

   
//=======RETRIEVE the DOFS of the COUPLED QUANTITIES    
 #if (BMAG_QTY==1)
  if ( Bext._eqnptr != NULL )  Bext.GetElemDofs(Level); 
  else                         Bext._qtyptr->FunctionDof(Bext,time,&xyz_refbox._val_dofs[0]);
  if ( Bhom._eqnptr != NULL )  Bhom.GetElemDofs(Level);   
  else                         Bhom._qtyptr->FunctionDof(Bhom,time,&xyz_refbox._val_dofs[0]);
#endif
#if (TEMP_QTY==1)
   if ( Temp._eqnptr != NULL ) Temp.GetElemDofs(Level);
     else                      Temp._qtyptr->FunctionDof(Temp,time,&xyz_refbox._val_dofs[0]);
#endif

//=== the connectivity is only related to the ELEMENT, so it is GEOMETRICAL
//===then, the DofMap is RELATED to the EQUATION the Vect comes from!
     //in fact, get_dof_val comes is related to the single equation
//===== After having all the dofs retrieved, we can start MANIPULATING them a bit!
     //Summing them    //Extending them,     //Whatever is performed on the Vect._val_dofs

  #if (BMAG_QTY==1)   
//======SUM Bhom and Bext  //from now on, you'll only use Bmag //Bmag,Bext and Bhom must have the same orders!
    Math::zeroN(&Bmag._val_dofs[0],Bmag._dim*Bmag._ndof);

    for (uint ivarq=0; ivarq < Bmag._dim; ivarq++)    { //ivarq is like idim
          for (uint d=0; d <  Bmag._ndof; d++)    {
          const uint     indxq  =         d + ivarq*Bmag._ndof;
          Bmag._val_dofs[indxq] = Bext._val_dofs[indxq] + Bhom._val_dofs[indxq];
	  }
    }
//=======after summing you EXTEND them to 3D
//for safety, i decide i'll put this INSIDE the CURLVECT, because it is the only routine 
//for which the _3D is used
//so you're sure you're doing it at a later time
//      ExtendDofs(Bmag); /*mgMHD->*//*mgMHDCONT->*/ // _utils.extend_nds(NDOF_FEM,B_qdofs,B_qdofs3D);//mgMHD2->_mixfe->_eldof[0][0]
//Do you need the DofExtension for other operators than CurlVect? Like for some vector products?
//if you have some multiplication v x phii , where phii is the test function?
//No, you already have the gauss values in that case

#endif

    double rho_nd =  1.;
    double  mu_nd =  1.;
    
//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
   const uint el_ngauss = _eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();
   
    for (uint qp = 0; qp < el_ngauss; qp++) {  

//=======here starts the "COMMON SHAPE PART"==================
// the call of these things should be related to the Operator
//also, the choice of the QuadratureRule should be dependent of the Involved Operators
//and the FE orders on which these Operators act
//after the COMMON SHAPE PART i dont want to have any dependence on qp anymore!
//these  phi and dphi are used for the different stages:
//BEFORE i: by the interpolation functions
//INSIDE i: for the tEST functions and derivatives
//INSIDE j: for the SHAPE functions and derivatives
//again, it should be the Operator to decide what functions to be called,
//for what FE ORDER and what DERIVATIVE ORDER
//if we decide that the PREPARATION of the tEST and SHAPE 
//of a certain Unknown are COMMON TO ALL,
//Then we must only concentrate on preparing the OTHER involved quantities in that Operator
for (uint fe = 0; fe < QL; fe++)     { 
    currgp.SetPhiElDofsFEVB_g (fe,qp);
    currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
}  
	  
const double      det = dt*currgp.JacVectVV_g(xyz);   //InvJac: is the same for both QQ and LL!
const double dtxJxW_g = det*_eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { 
  currgp.SetDPhiDxyzElDofsFEVB_g   (fe,qp);
  currgp.ExtendDphiDxyzElDofsFEVB_g(fe); 
}
//=======end of the "COMMON SHAPE PART"==================

//now we want to fill the element matrix and rhs with ONLY values AT GAUSS POINTS
//we can divide the values at Gauss points into three parts:
//1-constant values
//2-values which depend on (i)
//3-values which depend on (i,j)
//1- can be used for filling both Ke and Fe
//2- can be used to fill Fe and also Ke
//3- can be used only for filling Ke

//Here, before entering the (i,j) loop, you compute quantities that DO NOT depend on (i,j),
//therefore the ELEMENT SUM, GAUSS SUM And tEST/SHAPE sums have ALL been performed
//but, these quantities depend on idim and jdim, because they are  involved in multiplications with tEST and SHAPE functions

//Internal Quantities
      VelOld.val_g();   //fills _val_g, needs _val_dofs
      VelOld.grad_g();     //fills _grad_g, needs _val_dofs //TODO can we see the analogies between JacVectVV_g and grad_g?
//Advection all VelOld
      for (uint idim=0; idim<space_dim; idim++) { AdvRhs_g[idim]=0.;}
          for (uint idim=0; idim<space_dim; idim++)  {
	    for (uint b=0; b<space_dim; b++) {
	    AdvRhs_g[idim] += VelOld._val_g[b]*VelOld._grad_g[idim][b]; }   }   //TODO NONLIN // grad is [ivar][idim], i.e. [u v w][x y z]
//Divergence VelOld
	  double Div_g=0.;
          for (uint idim=0; idim<space_dim; idim++)    Div_g += VelOld._grad_g[idim][idim];     //TODO NONLIN

#if (BMAG_QTY==1)
//compute curlB
//what is needed for curlB, ie the dofs3D and the dphi3D, is done BEFORE for the  dphi3D and here for the dofs3D  
//      ExtendDphiDxyzElDofsFEVB_g(vb,Bmag._FEord/*FE_MAG*/);
     Bmag.curl_g(); 

//compute B
     Bmag.val_g();

//compute curlBxB
          Math::extend(&Bmag._val_g[0],&Bmag._val_g3D[0],space_dim);                    //fills _val_g3D
          Math::cross(&Bmag._curl_g3D[0],&Bmag._val_g3D[0],curlBXB_g3D);

//compute JxB
          Math::cross(Jext_g3D,&Bmag._val_g3D[0],JextXB_g3D);
#endif

#if (TEMP_QTY==1)
     Temp.val_g();
 #endif

#if (TEMP_DEPS==1)
  double rho_t=1.; density_ptr->Temp_dep(Temp._val_g[0],rho_t); rho_nd *= rho_t;
  double mu_t=1. ; density_ptr->Temp_dep(Temp._val_g[0],mu_t) ; mu_nd  *= mu_t;
#endif
       
#ifdef AXISYMX
      dtxJxW_g  *=yyg;  //what is the symmetry axis? I guess the x axis.
#endif

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
       for (uint i=0; i<qtyzero_ndof; i++)     {
//======="COMMON tEST PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
        //in this part we follow the Test Fncs of the FIRST QUANTITY of THIS Equation
        //since the Equation in this FE code is in WEAK FORM,
        //every Operator would have the test FUNCTIONS as ONE and ONLY ONE OF THEIR ARGUMENTS 
                                    const double phii_g       =      currgp._phi_ndsQLVB_g[qtyzero_ord][i];
        for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[qtyzero_ord][i+idim*qtyzero_ndof];
//======= END "COMMON tEST PART for QTYZERO" ==========
	
	  for (uint idim=0; idim<space_dim; idim++) {
            const uint irowq=i+idim*qtyzero_ndof;  //  (i):       dof of the tEST function
                                                  //(idim): component of the tEST function
           currelem.Rhs()(irowq) += 
         currelem.GetBCDofFlag()[irowq]*
           dtxJxW_g*(          NonStatNS*rho_nd*     VelOld._val_g[idim]*phii_g/dt  //time
                            + _AdvNew_fl*rho_nd*          AdvRhs_g[idim]*phii_g     //TODO NONLIN
                            +            rho_nd*IFr*gravity._val_g[idim]*phii_g     // gravity
#if (BMAG_QTY==1)
                            +                        S*curlBXB_g3D[idim]*phii_g   //TODO NONLIN external
                            +                           JextXB_g3D[idim]*phii_g  
#endif                            
                               )
            + (1-currelem.GetBCDofFlag()[irowq])*detb*VelOld._val_dofs[irowq] //Dirichlet bc    
	;
          }
           // end filling element rhs u

if (_Dir_pen_fl == 0)  { //faster than multiplying by _Dir_pen_fl
//actually this is not needed because you add only zeros
//it is only needed to avoid doing the operations on Ke
// Bc_ConvertToDirichletPenalty puts all ONE on the quadratic dofs
  for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
          const uint irowq = i+idim*qtyzero_ndof;
          currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
        }
}                                         // end filling diagonal for Dirichlet bc

//============ QTYZERO x QTYZERO dofs matrix (A matrix) ============
        for (uint j=0; j< qtyzero_ndof; j++) {
//======="COMMON SHAPE PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
//   (j):       dof of the SHAPE function
           double                                  phij_g       =     currgp._phi_ndsQLVB_g[qtyzero_ord][j];
           for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[qtyzero_ord][j+idim*qtyzero_ndof];
//======= END "COMMON SHAPE PART for QTYZERO" ==========
  
          double Lap_g=Math::dot(dphijdx_g,dphiidx_g,space_dim);
	  double Adv_g=Math::dot(&VelOld._val_g[0],dphijdx_g,space_dim);
          
          for (uint idim=0; idim<space_dim; idim++) { //filled in as 1-2-3 // 4-5-6 // 7-8-9
            int irowq = i+idim*qtyzero_ndof;      //(i) is still the dof of the tEST functions
                                                  //(idim): component of the tEST functions

            currelem.Mat()(irowq,j+idim*qtyzero_ndof)  // diagonal blocks [1-5-9] [idim(rows),idim(columns)]  //(idim): component of the SHAPE functions
               += 
            currelem.GetBCDofFlag()[irowq]*    
            dtxJxW_g*(
                   +NonStatNS*                           rho_nd*phij_g*phii_g/dt
                 + _AdvPic_fl*                           rho_nd* Adv_g*phii_g                //TODO NONLIN
                 + _AdvNew_fl*rho_nd*phij_g*VelOld._grad_g[idim][idim]*phii_g                //TODO NONLIN
                 + _AdvPic_fl*_Stab_fl*rho_nd*        0.5*Div_g*phij_g*phii_g                //TODO NONLIN
                 +                         mu_nd*IRe*(      dphijdx_g[idim]*dphiidx_g[idim] + Lap_g)
               );

            int idimp1=(idim+1)%space_dim;    // block +1 [2-6-7] [idim(rows),idim+1(columns)]  //(idimp1): component of the SHAPE functions
            currelem.Mat()(irowq,j+idimp1*qtyzero_ndof)
               +=
            currelem.GetBCDofFlag()[irowq]*
            dtxJxW_g*(
                   _AdvNew_fl*rho_nd*phij_g*VelOld._grad_g[idim][idimp1]*phii_g           //TODO NONLIN
                              +            mu_nd*IRe*(     dphijdx_g[idim]*dphiidx_g[idimp1])
               );
#if (DIMENSION==3)
	    
            int idimp2=(idim+2)%space_dim;   // block +2 [3-4-8] [idim(rows),idim+2(columns)]  //(idimp2): component of the SHAPE functions
            currelem.Mat()(irowq,j+idimp2*qtyzero_ndof)
               +=
           currelem.GetBCDofFlag()[irowq]*           
            dtxJxW_g*(
                   _AdvNew_fl*rho_nd*phij_g*VelOld._grad_g[idim][idimp2]*phii_g          //TODO NONLIN
                              +            mu_nd*IRe*(     dphijdx_g[idim]*dphiidx_g[idimp2])
               );
#endif
          }
 
        } 
//============ END QTYZERO x QTYZERO dofs matrix (A matrix) ============

//============ QTYZERO x QTYONE dofs matrix (B^T matrix) // ( p*div(v) ) (NS eq) ============
       for (uint j=0; j</*0*/qtyone_ndof; j++) {
//======="COMMON SHAPE PART for QTYONE" ==================
	 const double psij_g = currgp._phi_ndsQLVB_g[qtyone_ord][j];
//======="COMMON SHAPE PART for QTYONE" - END ============

          const int jclml = j + qtyZeroToOne_DofOffset;   /* space_dim*qtyzero_ndof */
          for (uint idim=0; idim<space_dim; idim++) {
            uint irowq = i+idim*qtyzero_ndof;
            currelem.Mat()(irowq,jclml) +=
           currelem.GetBCDofFlag()[irowq]*                    
            dtxJxW_g*(-psij_g*dphiidx_g[idim]);   /**   (-1.)*/
	    
           } 

        }
//============ END QTYZERO x QTYONE dofs matrix (B^T matrix) ============

          if (i</*0*/qtyone_ndof) {
//======="COMMON tEST PART for QTYONE" ============
          double psii_g = currgp._phi_ndsQLVB_g[qtyone_ord][i];
//======= "COMMON tEST PART for QTYONE" - END ============
	  const uint irowl = i + qtyZeroToOne_DofOffset;   /* space_dim*qtyzero_ndof */
          currelem.Rhs()(irowl)=0.;  // rhs
 //             currelem.Mat()(irowl,j+space_dim*qtyzero_ndof)  += (1./dt)*dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt;  //no bc here (KOMP dp/dt=rho*div)

          for (uint j=0; j<qtyzero_ndof; j++) { // B element matrix q*div(u)
//======="COMMON SHAPE PART for QTYZERO" ==================
            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[qtyzero_ord][j+idim*qtyzero_ndof];
//======="COMMON SHAPE PART for QTYZERO" - END ============
	    
            for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*qtyzero_ndof) += -/*(1./dt)**/dtxJxW_g*psii_g*dphijdx_g[idim]; 
                }

        }
                         // end pressure eq (cont)
      }
//===================================================================
//========= END FILLING ELEMENT MAT/RHS (i loop) =====================
//===================================================================
    }
//==============================================================
//================== END GAUSS LOOP (qp loop) ======================
//==============================================================
    
    ///  Add element matrix and rhs to the global ones.
                   _A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
                   _b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());

  } 
  // end of element loop
  

  }//END VOLUME
  
    // *****************************************************************
    // *****************************************************************
    // *****************************************************************

  {//BEGIN BOUNDARY  // *****************************************************************

    const uint mesh_vb = BB;
  
    CurrentElem       currelem(BB,this,_mesh,_eqnmap._elem_type);    
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,_eqnmap, _mesh.get_dim());
  
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity VelOld(currgp);
    VelOld._qtyptr   = _QtyInternalVector[QTYZERO]; //an alternative cannot exist, because it is an Unknown of This Equation
    VelOld.VectWithQtyFillBasic();   //the internal quantities will eventually have *this as eqn pointer
    VelOld.Allocate();

   const uint   qtyzero_ord  = VelOld._FEord;
   const uint   qtyzero_ndof = VelOld._ndof; 
    Velocity*  vel_castqtyptr = static_cast<Velocity*>(VelOld._qtyptr); //casting for quantity-specific functions

//=========
    CurrentQuantity pressOld(currgp);
    pressOld._qtyptr   = _QtyInternalVector[QTYONE];
    pressOld.VectWithQtyFillBasic();
    pressOld.Allocate();

   const uint qtyone_ord  = pressOld._FEord;
   const uint qtyone_ndof = pressOld._ndof; 

   //order
   const uint  qtyZeroToOne_DofOffset = VelOld._ndof*VelOld._dim;
   
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();

    //==================Quadratic domain, auxiliary, must be QUADRATIC!!! ==========
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = _mesh.GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
  xyz_refbox.Allocate();

  
    const uint nel_e = _mesh._off_el[mesh_vb][_NoLevels*myproc+Level+1];
    const uint nel_b = _mesh._off_el[mesh_vb][_NoLevels*myproc+Level];

//=======================
//=======================    
    
//=== auxiliary Operators at the boundary
  double strainU_g[DIMENSION][DIMENSION];
  double strainUtrDn_g[DIMENSION];


  for (uint iel=0; iel < (nel_e - nel_b) ; iel++) {

     currelem.Mat().zero();
     currelem.Rhs().zero();

     currelem.set_el_nod_conn_lev_subd(Level,myproc,iel);
     currelem.SetMidpoint();
     
     currelem.ConvertElemCoordsToMappingOrd(xyz);
     _mesh.TransformElemNodesToRef(currelem.GetDim(),currelem.GetNodeCoords(),&xyz_refbox._val_dofs[0]);    

     currelem.SetElDofsBc(Level);
     
     VelOld.GetElemDofs(Level);
     pressOld.GetElemDofs(Level);

    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),qtyzero_ord,currelem.GetBCDofFlag()); //only the Quadratic Part is modified! /*OK DIR_PEN*/
       

//============ BC =======
       int     el_flag[NT] = {0,0}; //normal and tangential flag
#if DIMENSION==2
       double el_value[N1T1] = {0.,0.};  //1 normal and 1 tangential
#elif DIMENSION==3
       double el_value[N1T3] = {0.,0.,0.,0.}; //1 normal and 3 tangential
#endif
       double  dbl_pen[NT] = {0.,0.};  //normal and tangential penalty value
   
    
       Bc_GetElFlagValLevSubd(Level,myproc,iel,el_flag,el_value);

if (_Dir_pen_fl == 1)  { 
       if (el_flag[NN] == 1) {   dbl_pen[NN]=penalty_val; } //normal dirichlet
       if (el_flag[TT] == 1) {   dbl_pen[TT]=penalty_val; } //tangential dirichlet
   }

//checks
//TODO here i should check that the nodal bc dirichlet i put correspond to the element NT flags

       uint press_fl=0;
       Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(currelem.GetBCDofFlag(),VelOld,pressOld,press_fl); //compute the PRESSURE FLAG with the PRESSURE nodal bc flags
 //only the LINEAR PART is USED!!
       
// // TODO  if ( (1-el_flag[NN]) != press_fl)  {std::cout << "Sthg wrong with press elflags" << std::endl;abort();}

//========END BC============

//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
   const uint el_ngauss = _eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();
   
    for (uint qp=0; qp< el_ngauss; qp++) {
            
//======= "COMMON SHAPE PART"============================
for (uint fe = 0; fe < QL; fe++)     {    
  currgp.SetPhiElDofsFEVB_g (fe,qp);
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp); 
}

        const double det   = dt*currgp.JacVectBB_g(xyz);
	const double dtxJxW_g = det * _eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
//=======end "COMMON SHAPE PART"===================================

//-------- pressure==============
    //"predof" VS "post gauss method"
    //pay attention to this fact: even if we are in the equation to which pressure is associated
    //(pressure is an unknown here!), we might prefer using the FUNCTION instead,
    //in order to get some values that DO NOT CHANGE with the solution
    //so i'm using the FUNCTION of an UNKNOWN, it is not an external quantity  
	  //clearly, this is ALTERNATIVE to the command
	  
// // //    val_g(vb,pressOld);
   
   xyz_refbox.val_g(); // val_g(vb,xyz);   //CHECK the QUADRATICS!!!!!!!!!
      pressOld._qtyptr->Function_txyz(time,&xyz_refbox._val_g[0]/*xyz._val_g*/,&pressOld._val_g[0]);  //i prefer using the function instead of the p_old vector
       
//--- strain, derivative of velocity ============== 
      
// // // TODO    /*VelOld._qtyptr*/vel_castqtyptr->strain_txyz_box(time,xyz_refbox._val_g/*xyz._val_g*/,strainU_g);
// // //          for (uint idim=0; idim< /*space_dim*/; idim++)  strainUtrDn_g[idim]=0.;
// // // 	    for (uint idim=0; idim< /*space_dim*/; idim++)  {   //beware that this space_dim is not IntDim
// // // 	      for (uint jdim=0; jdim< /*space_dim*/; jdim++)    {
// // //                 strainUtrDn_g[idim] += strainU_g[jdim][idim]*get_normal_ptr()[jdim];
// // // 	        }
// // // 	    }
//=================================================

      //-----old velocity=============
	  VelOld.val_g();  //substitute el_value...


//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
      for (uint i=0; i<qtyzero_ndof; i++) { // rhs (NS eq) //interpolation of the velocity test functions on the boundary

	const double phii_g = currgp._phi_ndsQLVB_g[qtyzero_ord][i];

         for (uint idim=0; idim< space_dim; idim++)    {
             uint irowq=i+idim*qtyzero_ndof;
            currelem.Rhs()(irowq)  += 
          currelem.GetBCDofFlag()[irowq]*           
           dtxJxW_g*(   -1.*/*press_fl*/(1-el_flag[NN])*pressOld._val_g[0]*currgp.get_normal_ptr()[idim]*phii_g  //  //OLD VALUES //AAA multiplying int times uint!!!

// // //             TODO STRAIN AT THE BOUNDARY            + /*stress_fl*/el_flag[1]*IRe*strainUtrDn_g[idim]*phii_g 

	  )
                                //projection over the physical (x,y,z)
      + _Dir_pen_fl *dtxJxW_g*phii_g*(dbl_pen[NN]*el_value[0]*currgp.get_normal_ptr()[idim] 
                                    + dbl_pen[TT]*el_value[1]*currgp.get_tangent_ptr()[0][idim]  // VelOld._val_g[idim] instead of el_value...
               #if DIMENSION==3
		                    + dbl_pen[TT]*el_value[1]*currgp.get_tangent_ptr()[1][idim]    
               #endif   
          )
	   ;   
	   
//====================
if (_Dir_pen_fl == 1) {  //much faster than multiplying by _Dir_pen_fl=0 , and much better than removing the code with the #ifdef //  #if (NS_DIR_PENALTY==1)  
	   for (uint jdim=0; jdim< space_dim; jdim++)    {

	   for (uint j=0; j<qtyzero_ndof; j++) {
          const double phij_g = currgp._phi_ndsQLVB_g[qtyzero_ord][j];

  currelem.Mat()(irowq,j+jdim*qtyzero_ndof) +=                //projection over the physical (x,y,z) 
      + /*_Dir_pen_fl**/dtxJxW_g*phii_g*phij_g*(dbl_pen[NN]*currgp.get_normal_ptr()[jdim]*currgp.get_normal_ptr()[idim]   //the PENALTY is BY ELEMENT, but the (n,t) is BY GAUSS because we cannot compute now a nodal normal
                                              + dbl_pen[TT]*currgp.get_tangent_ptr()[0][jdim]*currgp.get_tangent_ptr()[0][idim]
                 #if DIMENSION==3
                                              + dbl_pen[TT]*currgp.get_tangent_ptr()[1][jdim]*currgp.get_tangent_ptr()[1][idim]
                #endif   
                   );
	         } //end j
      
               } //end jdim
  
             }  //end penalty if
//====================

	 }
           //end of idim loop

     }
//==============================================================
//========= END FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
    } 
//==================================================================
//================== END GAUSS LOOP (qp loop) ======================
//==================================================================
    
    _A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
    _b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());

    
  }
  // end of BDRYelement loop

  
    }
  // END BOUNDARY ******************************
  

#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << _eqname << ": assembled  Level " << Level
           << " with " << _A[Level]->m() << " dofs" << std::endl;
#endif

    return;
}


} //end namespace femus



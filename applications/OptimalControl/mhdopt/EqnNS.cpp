
#include "FemusDefault.hpp"

#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "DenseVector.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"

#include "Math.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "MultiLevelProblem.hpp"
#include "FETypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "VBTypeEnum.hpp"
#include "GeomElTypeEnum.hpp"
#include "Domain.hpp"
#include "TimeLoop.hpp"
#include "CurrentElem.hpp"

//application
#include "OptLoop.hpp"
#include "OptQuantities.hpp"
#include "Box.hpp"


    
  void GenMatRhsNS(MultiLevelProblem &ml_prob)  {

  SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_NS"); 
  const unsigned Level = my_system.GetLevelToAssemble();
    
    
  //====== reference values ========================
//====== related to Quantities on which Operators act, and to the choice of the "LEADING" EQUATION Operator
  //====== Physics
  const double IRe = 1./ml_prob.GetInputParser().get("Re");
  const double IFr = 1./ml_prob.GetInputParser().get("Fr");
  const double   S = ml_prob.GetInputParser().get("S");
//================================================  

//=============== electric current ===============
  //this is a particular Vect. It has already been defined in 3D
  //because we only need it for a CROSS product
  double Jext_g3D[3]={0.,0.,0.}; //Quantity
  //=======density and viscosity===================
  const double rhof = ml_prob.GetInputParser().get("rho0");
  const double  muf = ml_prob.GetInputParser().get("mu0");

//================================================  
//=======Operators @ gauss =======================
  double    dphijdx_g[DIMENSION]; // ShapeDer(): used for Laplacian,Divergence, ..., associated to an Unknown Quantity
  double    dphiidx_g[DIMENSION]; // Test(): used for Laplacian,Advection,Divergence... associated to an Unknown Quantity
  double     AdvRhs_g[DIMENSION]; //Operator: Adv(u,u,phi)
  double  curlBXB_g3D[3];         //Operator  CrossProd(curlB,B)
  double   JextXB_g3D[3];         //Operator: CrossProd(Jext,B)
//================================================  


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
const int NonStatNS = (int) ml_prob.GetInputParser().get("NonStatNS");
  const double   dt = 1.; //ml_prob._timeloop._timemap.get("dt");

//==========FLAG for NONLINEARITY, we might put here
   const uint   _AdvPic_fl = 1;
   const uint   _AdvNew_fl = 0;
   const uint   _Stab_fl = 0;
   const double _Komp_fac = 0.;
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
  const uint           space_dim = ml_prob._ml_msh->GetDimension();
  
        my_system._LinSolver[Level]->_KK->zero();
        my_system._LinSolver[Level]->_RESC->zero();

// ==========================================  
  Mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(Level);
  elem		*myel		=  mymsh->el;
  const unsigned myproc  = mymsh->processor_id();
	
// ==========================================  
// ==========================================  

   {//BEGIN VOLUME
//===============================  
//===============================  
  
    const uint mesh_vb = VV;
 
    const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level+1];
    const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level];
    
    for (uint iel=0; iel < (nel_e - nel_b); iel++) {
      
    CurrentElem<double>       currelem(iel,myproc,Level,VV,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType(),mymsh);    
    //   CurrentGaussPointBase & currgp = //   CurrentGaussPointBase::build(currelem,ml_prob.GetQuadratureRule(currelem.GetDim()));
  
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
      CurrentQuantity VelOldX(currelem);
    VelOldX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity0"); 
    VelOldX.VectWithQtyFillBasic();
    VelOldX.Allocate();
    
    CurrentQuantity VelOldY(currelem);
    VelOldY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity1"); 
    VelOldY.VectWithQtyFillBasic();
    VelOldY.Allocate();
    
    CurrentQuantity VelOldZ(currelem);
    VelOldZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity2"); 
    VelOldZ.VectWithQtyFillBasic();
    VelOldZ.Allocate();
    
    std::vector<CurrentQuantity*> VelOld_vec;   
    VelOld_vec.push_back(&VelOldX);
    VelOld_vec.push_back(&VelOldY);
    VelOld_vec.push_back(&VelOldZ);
    
   const uint   qtyzero_ord  = VelOldX._FEord;
   const uint   qtyzero_ndof = VelOldX._ndof; 

//=========
    CurrentQuantity pressOld(currelem);
    pressOld._qtyptr = ml_prob.GetQtyMap().GetQuantity("Qty_Pressure");
    pressOld.VectWithQtyFillBasic();
    pressOld.Allocate();

   const uint qtyone_ord  = pressOld._FEord;
   const uint qtyone_ndof = pressOld._ndof; 

   //order
   const uint  qtyZeroToOne_DofOffset = VelOldX._ndof*space_dim;
   
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currelem);
    xyz._dim      = DIMENSION;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

//============================ MAG WORLD =======================================
    CurrentQuantity BhomX(currelem);
    BhomX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHom0"); 
    BhomX.VectWithQtyFillBasic();
    BhomX.Allocate();
    
    CurrentQuantity BhomY(currelem);
    BhomY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHom1"); 
    BhomY.VectWithQtyFillBasic();
    BhomY.Allocate();
    
    CurrentQuantity BhomZ(currelem);
    BhomZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHom2"); 
    BhomZ.VectWithQtyFillBasic();
    BhomZ.Allocate();
    
    std::vector<CurrentQuantity*> Bhom_vec;   
    Bhom_vec.push_back(&BhomX);
    Bhom_vec.push_back(&BhomY);
    Bhom_vec.push_back(&BhomZ);
 
    //Bhom only to retrieve the dofs

    CurrentQuantity BextX(currelem);
    BextX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt0"); 
    BextX.VectWithQtyFillBasic();
    BextX.Allocate();
    
    CurrentQuantity BextY(currelem);
    BextY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt1"); 
    BextY.VectWithQtyFillBasic();
    BextY.Allocate();
    
    CurrentQuantity BextZ(currelem);
    BextZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt2"); 
    BextZ.VectWithQtyFillBasic();
    BextZ.Allocate();
    
    std::vector<CurrentQuantity*> Bext_vec;   
    Bext_vec.push_back(&BextX);
    Bext_vec.push_back(&BextY);
    Bext_vec.push_back(&BextZ);
 
    //Bext only to retrieve the dofs

//========= auxiliary, must be AFTER Bhom!
    CurrentQuantity Bmag(currelem); //total
    Bmag._dim        = Bhom_vec.size();
    Bmag._FEord      = BhomX._FEord;
    Bmag._ndof       = BhomX._ndof;
    Bmag.Allocate();
//======================== MAG WORLD ================================

//=======gravity==================================
  CurrentQuantity gravity(currelem);
  gravity._dim=DIMENSION;
  gravity._val_g.resize(gravity._dim);
  gravity._val_g[0] = ml_prob.GetInputParser().get("dirgx");
  gravity._val_g[1] = ml_prob.GetInputParser().get("dirgy");
#if DIMENSION==3
  gravity._val_g[2] = ml_prob.GetInputParser().get("dirgz");
#endif  

//=======================
//=======================    

 
    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.SetDofobjConnCoords();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);

//=======RETRIEVE the DOFS of the UNKNOWN QUANTITIES,i.e. MY EQUATION
    currelem.SetElDofsBc();
    
    for (uint idim=0; idim < space_dim; idim++) VelOld_vec[idim]->GetElemDofs();
     
    pressOld.GetElemDofs();

//=======RETRIEVE the DOFS of the COUPLED QUANTITIES    
  
  for (uint idim=0; idim < space_dim; idim++) { 
         Bhom_vec[idim]->GetElemDofs();   
         Bext_vec[idim]->GetElemDofs(); 
     }

//=== the connectivity is only related to the ELEMENT, so it is GEOMETRICAL
//===then, the DofMap is RELATED to the EQUATION the Vect comes from!
     //in fact, get_dof_val comes is related to the single equation
//===== After having all the dofs retrieved, we can start MANIPULATING them a bit!
     //Summing them    //Extending them,     //Whatever is performed on the Vect._val_dofs

//======SUM Bhom and Bext  //from now on, you'll only use Bmag //Bmag,Bext and Bhom must have the same orders!
    Math::zeroN(&Bmag._val_dofs[0],Bmag._dim*Bmag._ndof);

    for (uint ivarq=0; ivarq < Bmag._dim; ivarq++)    { //ivarq is like idim
          for (uint d=0; d <  Bmag._ndof; d++)    {
          const uint     indxq  =         d + ivarq*Bmag._ndof;
          Bmag._val_dofs[indxq] = Bext_vec[ivarq]->_val_dofs[d] + Bhom_vec[ivarq]->_val_dofs[d];
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


    double rho_nd =  1.;
    double  mu_nd =  1.;
    
//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
   const uint el_ngauss = ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussPointsNumber();
   
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
//     currgp.SetPhiElDofsFEVB_g (fe,qp);
//     currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
}  
	  
const double      det = dt*1.; //currgp.JacVectVV_g(xyz);
const double dtxJxW_g = det*ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussWeight(qp);
const double     detb = det/el_ngauss;
	  
for (uint fe = 0; fe < QL; fe++)     { 
//   currgp.SetDPhiDxyzElDofsFEVB_g   (fe,qp);
//   currgp.ExtendDphiDxyzElDofsFEVB_g(fe); 
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
      for (uint idim=0; idim < space_dim; idim++) { 
	VelOld_vec[idim]->val_g(); 
	VelOld_vec[idim]->grad_g();
      }

//Advection all VelOld
      for (uint idim=0; idim<space_dim; idim++) { AdvRhs_g[idim]=0.;}
          for (uint idim=0; idim<space_dim; idim++)  {
	    for (uint b=0; b<space_dim; b++) {
	    AdvRhs_g[idim] += VelOld_vec[b]->_val_g[0]*VelOld_vec[idim]->_grad_g[0][b]; }   }   // grad is [ivar][idim], i.e. [u v w][x y z]
//Divergence VelOld
	  double Div_g=0.;
          for (uint idim=0; idim<space_dim; idim++)    Div_g += VelOld_vec[idim]->_grad_g[0][idim];


	  
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

	  
	  
//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
// // //        for (uint i=0; i<qtyzero_ndof; i++)     {
// // // //======="COMMON tEST PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
// // //         //in this part we follow the Test Fncs of the FIRST QUANTITY of THIS Equation
// // //         //since the Equation in this FE code is in WEAK FORM,
// // //         //every Operator would have the test FUNCTIONS as ONE and ONLY ONE OF THEIR ARGUMENTS 
// // //                                     const double phii_g       =      currgp._phi_ndsQLVB_g[qtyzero_ord][i];
// // //         for (uint idim=0; idim<space_dim; idim++)  dphiidx_g[idim] = currgp._dphidxyz_ndsQLVB_g[qtyzero_ord][i+idim*qtyzero_ndof];
// // // //======= END "COMMON tEST PART for QTYZERO" ==========
// // // 	
// // // 	  for (uint idim=0; idim<space_dim; idim++) {
// // //             const uint irowq=i+idim*qtyzero_ndof;  //  (i):       dof of the tEST function
// // //                                                   //(idim): component of the tEST function
// // //            currelem.Rhs()(irowq) += 
// // //          currelem.GetBCDofFlag()[irowq]*
// // //            dtxJxW_g*(          NonStatNS*rho_nd*     VelOld_vec[idim]->_val_g[0]*phii_g/dt  //time
// // //                             + _AdvNew_fl*rho_nd*          AdvRhs_g[idim]*phii_g     //TODO NONLIN
// // //                             +            rho_nd*IFr*gravity._val_g[idim]*phii_g     // gravity
// // //                             
// // //                             +                        S*curlBXB_g3D[idim]*phii_g   //TODO NONLIN external
// // //                             +                           JextXB_g3D[idim]*phii_g  
// // //                                )
// // //             + (1-currelem.GetBCDofFlag()[irowq])*detb*VelOld_vec[idim]->_val_dofs[i] //Dirichlet bc    
// // // 	;
// // //           }
// // // 
// // //   for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
// // //           const uint irowq = i+idim*qtyzero_ndof;
// // //           currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
// // //         }
// // // 
// // // //============ QTYZERO x QTYZERO dofs matrix (A matrix) ============
// // //         for (uint j=0; j< qtyzero_ndof; j++) {
// // // //======="COMMON SHAPE PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
// // // //   (j):       dof of the SHAPE function
// // //            double                                  phij_g       =     currgp._phi_ndsQLVB_g[qtyzero_ord][j];
// // //            for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[qtyzero_ord][j+idim*qtyzero_ndof];
// // // //======= END "COMMON SHAPE PART for QTYZERO" ==========
// // //   
// // //           double Lap_g=Math::dot(dphijdx_g,dphiidx_g,space_dim);
// // // 	  double Adv_g=0.;
// // // 	  for (uint idim=0; idim<space_dim; idim++) Adv_g += VelOld_vec[idim]->_val_g[0] * dphijdx_g[idim]; // =Math::dot(&VelOld._val_g[0],dphijdx_g,space_dim);
// // //           
// // //           for (uint idim=0; idim<space_dim; idim++) { //filled in as 1-2-3 // 4-5-6 // 7-8-9
// // //             int irowq = i+idim*qtyzero_ndof;      //(i) is still the dof of the tEST functions
// // //                                                   //(idim): component of the tEST functions
// // // 
// // //             currelem.Mat()(irowq,j+idim*qtyzero_ndof)  // diagonal blocks [1-5-9] [idim(rows),idim(columns)]  //(idim): component of the SHAPE functions
// // //                += 
// // //             currelem.GetBCDofFlag()[irowq]*    
// // //             dtxJxW_g*(
// // //                    +NonStatNS*                           rho_nd*phij_g*phii_g/dt
// // //                  + _AdvPic_fl*                           rho_nd* Adv_g*phii_g                //TODO NONLIN
// // //                  + _AdvNew_fl*rho_nd*phij_g*VelOld_vec[idim]->_grad_g[0][idim]*phii_g                //TODO NONLIN
// // //                  + _AdvPic_fl*_Stab_fl*rho_nd*        0.5*Div_g*phij_g*phii_g                //TODO NONLIN
// // //                  +                         mu_nd*IRe*(      dphijdx_g[idim]*dphiidx_g[idim] + Lap_g)
// // //                );
// // // 
// // //             int idimp1=(idim+1)%space_dim;    // block +1 [2-6-7] [idim(rows),idim+1(columns)]  //(idimp1): component of the SHAPE functions
// // //             currelem.Mat()(irowq,j+idimp1*qtyzero_ndof)
// // //                +=
// // //             currelem.GetBCDofFlag()[irowq]*
// // //             dtxJxW_g*(
// // //                    _AdvNew_fl*rho_nd*phij_g*VelOld_vec[idim]->_grad_g[0][idimp1]*phii_g           //TODO NONLIN
// // //                               +            mu_nd*IRe*(     dphijdx_g[idim]*dphiidx_g[idimp1])
// // //                );
// // // #if (DIMENSION==3)
// // // 	    
// // //             int idimp2=(idim+2)%space_dim;   // block +2 [3-4-8] [idim(rows),idim+2(columns)]  //(idimp2): component of the SHAPE functions
// // //             currelem.Mat()(irowq,j+idimp2*qtyzero_ndof)
// // //                +=
// // //            currelem.GetBCDofFlag()[irowq]*           
// // //             dtxJxW_g*(
// // //                    _AdvNew_fl*rho_nd*phij_g*VelOld_vec[idim]->_grad_g[0][idimp2]*phii_g          //TODO NONLIN
// // //                               +            mu_nd*IRe*(     dphijdx_g[idim]*dphiidx_g[idimp2])
// // //                );
// // // #endif
// // //           }
// // //  
// // //         } 
// // // //============ END QTYZERO x QTYZERO dofs matrix (A matrix) ============
// // // 
// // // //============ QTYZERO x QTYONE dofs matrix (B^T matrix) // ( p*div(v) ) (NS eq) ============
// // //        for (uint j=0; j</*0*/qtyone_ndof; j++) {
// // // //======="COMMON SHAPE PART for QTYONE" ==================
// // // 	 const double psij_g = currgp._phi_ndsQLVB_g[qtyone_ord][j];
// // // //======="COMMON SHAPE PART for QTYONE" - END ============
// // // 
// // //           const int jclml = j + qtyZeroToOne_DofOffset;   /* space_dim*qtyzero_ndof */
// // //           for (uint idim=0; idim<space_dim; idim++) {
// // //             uint irowq = i+idim*qtyzero_ndof;
// // //             currelem.Mat()(irowq,jclml) +=
// // //            currelem.GetBCDofFlag()[irowq]*                    
// // //             dtxJxW_g*(-psij_g*dphiidx_g[idim]);   /**   (-1.)*/
// // // 	    
// // //            } 
// // // 
// // //         }
// // // //============ END QTYZERO x QTYONE dofs matrix (B^T matrix) ============
// // // 
// // //           if (i</*0*/qtyone_ndof) {
// // // //======="COMMON tEST PART for QTYONE" ============
// // //           double psii_g = currgp._phi_ndsQLVB_g[qtyone_ord][i];
// // // //======= "COMMON tEST PART for QTYONE" - END ============
// // // 	  const uint irowl = i + qtyZeroToOne_DofOffset;   /* space_dim*qtyzero_ndof */
// // //           currelem.Rhs()(irowl)=0.;  // rhs
// // //  //             currelem.Mat()(irowl,j+space_dim*qtyzero_ndof)  += (1./dt)*dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt;  //no bc here (KOMP dp/dt=rho*div)
// // // 
// // //           for (uint j=0; j<qtyzero_ndof; j++) { // B element matrix q*div(u)
// // // //======="COMMON SHAPE PART for QTYZERO" ==================
// // //             for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[qtyzero_ord][j+idim*qtyzero_ndof];
// // // //======="COMMON SHAPE PART for QTYZERO" - END ============
// // // 	    
// // //             for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*qtyzero_ndof) += -/*(1./dt)**/dtxJxW_g*psii_g*dphijdx_g[idim]; 
// // //                 }
// // // 
// // //         }
// // //                          // end pressure eq (cont)
// // //       }
//===================================================================
//========= END FILLING ELEMENT MAT/RHS (i loop) =====================
//===================================================================
    }
//==============================================================
//================== END GAUSS LOOP (qp loop) ======================
//==============================================================
    
    ///  Add element matrix and rhs to the global ones.
                   my_system._LinSolver[Level]->_KK->add_matrix(currelem.Mat(),currelem.GetDofIndices());
                   my_system._LinSolver[Level]->_RESC->add_vector(currelem.Rhs(),currelem.GetDofIndices());

  } 
  // end of element loop
  

  }//END VOLUME
  
    // *****************************************************************
    // *****************************************************************
    // *****************************************************************

  {//BEGIN BOUNDARY  // *****************************************************************

    const uint mesh_vb = BB;
  
    const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level+1];
    const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level];

  for (uint iel=0; iel < (nel_e - nel_b) ; iel++) {
  
    CurrentElem<double>       currelem(iel,myproc,Level,BB,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType(),mymsh);    
    //   CurrentGaussPointBase & currgp = //   CurrentGaussPointBase::build(currelem,ml_prob.GetQuadratureRule(currelem.GetDim()));
  
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity VelOldX(currelem);
    VelOldX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity0"); 
    VelOldX.VectWithQtyFillBasic();
    VelOldX.Allocate();

   const uint   qtyzero_ord  = VelOldX._FEord;
   const uint   qtyzero_ndof = VelOldX._ndof; 

    CurrentQuantity pressOld(currelem);
    pressOld._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_Pressure");
    pressOld.VectWithQtyFillBasic();
    pressOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currelem);
    xyz._dim      = DIMENSION;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

//=======================    
    

     currelem.Mat().zero();
     currelem.Rhs().zero();

     currelem.SetDofobjConnCoords();
     currelem.SetElDofsBc();
     
     currelem.ConvertElemCoordsToMappingOrd(xyz);

     
     pressOld.GetElemDofs();

//============ BC =======
       int press_fl = currelem.Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(VelOldX._ndof,space_dim,pressOld); 
//========END BC============

//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
   const uint el_ngauss = ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussPointsNumber();
   
    for (uint qp=0; qp< el_ngauss; qp++) {
            
//======= "COMMON SHAPE PART"============================
for (uint fe = 0; fe < QL; fe++)     {    
//   currgp.SetPhiElDofsFEVB_g (fe,qp);
//   currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp); 
}

        const double det   = 1.;//dt*currgp.JacVectBB_g(xyz);
	const double dtxJxW_g = det * ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussWeight(qp);
//=======end "COMMON SHAPE PART"===================================

   xyz.val_g(); 
   
//      std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
//      std::cout << xyz._val_g[1] << std::endl;

   
   pressOld._qtyptr->Function_txyz(0.,&xyz._val_g[0],&pressOld._val_g[0]);  //i prefer using the function instead of the p_old vector
//        pressOld.val_g();  //this is the alternative
   if (press_fl == 1) {
     std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
     std::cout << pressOld._val_g[0] << std::endl;
   }
   
//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
// // //       for (uint i=0; i<qtyzero_ndof; i++) { // rhs (NS eq) //interpolation of the velocity test functions on the boundary
// // // 
// // // 	const double phii_g = currgp._phi_ndsQLVB_g[qtyzero_ord][i];
// // // 
// // //          for (uint idim=0; idim< space_dim; idim++)    {
// // //              uint irowq=i+idim*qtyzero_ndof;
// // //             currelem.Rhs()(irowq)  += 
// // //           currelem.GetBCDofFlag()[irowq]*           
// // //            dtxJxW_g*(   -1.*pressOld._val_g[0]*currgp.get_normal_ptr()[idim]*phii_g/**press_fl*/  //TODO if you uncomment this press_fl, which I think you should, it gives a different result...
// // // 
// // // 	  ) ;   
// // // 	   
// // // 	 }
// // //            //end of idim loop
// // //      }
//==============================================================
//========= END FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
    } 
//==================================================================
//================== END GAUSS LOOP (qp loop) ======================
//==================================================================
//     std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
//     currelem.Rhs().print_scientific(std::cout);
//     std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;

    my_system._LinSolver[Level]->_KK->add_matrix(currelem.Mat(),currelem.GetDofIndices());
    my_system._LinSolver[Level]->_RESC->add_vector(currelem.Rhs(),currelem.GetDofIndices());

    
  }
  // end of BDRYelement loop

  
    }
  // END BOUNDARY ******************************
  
        my_system._LinSolver[Level]->_KK->close();
        my_system._LinSolver[Level]->_RESC->close();

#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << my_system.name() << ": assembled  Level " << Level
           << " with " << my_system._LinSolver[Level]->_KK->m() << " dofs" << std::endl;
#endif

    return;
}






/// This function generates the initial conditions for the NS system:
//xp[] is the NON-DIMENSIONAL node coordinate
//when this function is called,
//the domain has been NONDIMENSIONALIZED
//but NOT ROTATED YET
//on the other hand, the functions of the type _txyz
//only accept a ROTATED xp, so let us not forget about rotating

 //this law is given in the NON-DIMENSIONAL domain, with NON-DIMENSIONAL values
 //NONDIMENSIONAL pressure distribution, fundamental!!
 

//=======================
///COMMENTS on the BC for MHD ====================

/// This function defines the boundary conditions for the MHD system:
    
 //Bxn = (b+Be)xn = 0 would mean perfectly insulating   
 //our present situation is such that:
//on the walls:
//bxn = Bexn = 0 => Bxn = 0
//on the inlet and outlet we cannot speak of WALLS
//how do we set BOUNDARY/INTERFACE conditions on the INLET/OUTLET?
//right out of the boundary there will still be the same material
//so we have to put INFINITY|SYMMETRY conditions
//because we are studying an infinite channel:
//   dB/dy=0
//that would also come out from a formulation with the Laplacian
//So you dont have to fix the test function, because the integrand 
//is specified (equal to zero... but specified)
// So for symmetry you LET FREE all the quantities
// and that's ok. That means that the possible
//BOUNDARY INTEGRALS containing FIRST DERIVATIVES 
//of as a form of NORMAL DERIVATIVE, CURL DERIVATIVE or whatever
//are set TO ZERO, 
//not because the projection onto the normal is ZERO
//but because you assume that ALL THE NORMAL DERIVATIVES ARE ZERO.

//however we have
 //bxn = 0
 //Bexn = either FIXED or CONTROLLED value


  



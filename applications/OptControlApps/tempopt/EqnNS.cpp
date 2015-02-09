
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
#include "QTYnumEnum.hpp"
#include "Domain.hpp"
#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"
#include "Box.hpp"
  

//application
#include "OptLoop.hpp"
#include "TempQuantities.hpp"


//===================================================
/// This function assembles the matrix and the rhs:
 void GenMatRhsNS(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assemble_matrix)  {

   SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_NS");
    
   const uint   _AdvPic_fl = 1;
   const uint   _AdvNew_fl = 0;
   const uint   _Stab_fl = 0;
   const double _Komp_fac = 0.;

   const double time =  0.;//ml_prob._timeloop._curr_time;
   
//========== PROCESSOR INDEX
  const uint myproc = ml_prob.GetMeshTwo()._iproc;

//==========FLAG FOR STATIONARITY OR NOT
  const int NonStatNS = (int) ml_prob.GetInputParser().get("NonStatNS");
  const double     dt = 1.;

//========= BCHandling =========
  const double penalty_val = ml_prob.GetMeshTwo().GetRuntimeMap().get("penalty_val");    

  //========== GEOMETRIC ELEMENT ========
  const uint           space_dim = ml_prob.GetMeshTwo().get_dim();
  const uint mesh_ord = (int) ml_prob.GetMeshTwo().GetRuntimeMap().get("mesh_ord");
  const uint meshql   = (int) ml_prob.GetMeshTwo().GetRuntimeMap().get("meshql"); //======== ELEMENT MAPPING =======


  //=======density and viscosity===================
  const double rhof = ml_prob.GetInputParser().get("rho0");
  const double  muf = ml_prob.GetInputParser().get("mu0");
  //====== reference values ========================
//====== related to Quantities on which Operators act, and to the choice of the "LEADING" EQUATION Operator
  const double IRe = 1./ml_prob.GetInputParser().get("Re");
  const double IFr = 1./ml_prob.GetInputParser().get("Fr");
//================================================  

//================================================  
//=======Operators @ gauss =======================
  std::vector<double>    dphijdx_g(space_dim); // ShapeDer(): used for Laplacian,Divergence, ..., associated to an Unknown Quantity
  std::vector<double>    dphiidx_g(space_dim); // Test(): used for Laplacian,Advection,Divergence... associated to an Unknown Quantity
  std::vector<double>     AdvRhs_g(space_dim); //Operator: Adv(u,u,phi)
//================================================  

        my_system._A[Level]->zero();
        my_system._b[Level]->zero();


  {//BEGIN VOLUME
//========================
//========================
    const uint mesh_vb = VV;
  
    const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level+1];
    const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level];

  for (int iel=0; iel < (nel_e - nel_b); iel++) {
    
    CurrentElem       currelem(Level,VV,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType());
    
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,ml_prob.GetQrule(currelem.GetDim()));
 
  
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity VelOld(currgp);
    VelOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYZERO]; //an alternative cannot exist, because it is an Unknown of This Equation
    VelOld.VectWithQtyFillBasic();
    VelOld.Allocate();

   const uint   qtyzero_ord  = VelOld._FEord;
   const uint   qtyzero_ndof = VelOld._ndof; 
//     Velocity*  vel_castqtyptr = static_cast<Velocity*>(VelOld._qtyptr); //casting for quantity-specific functions

//=========
    CurrentQuantity pressOld(currgp);
    pressOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYONE];
    pressOld.VectWithQtyFillBasic();
    pressOld.Allocate();

   const uint qtyone_ord  = pressOld._FEord;
   const uint qtyone_ndof = pressOld._ndof; 

   //order
   const uint  qtyZeroToOne_DofOffset = VelOld._ndof*VelOld._dim;
   
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp); //domain
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();
    
    //==================Quadratic domain, auxiliary, must be QUADRATIC!!! ==========
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = NVE[ ml_prob.GetMeshTwo()._geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];
  xyz_refbox.Allocate();
    
//other Physical constant Quantities
//=======gravity==================================
  CurrentQuantity gravity(currgp);
  gravity._dim = space_dim;
//   gravity.Allocate(); CANNOT DO THIS NOW BECAUSE NOT ALL THE DATA FOR THE ALLOCATION ARE FILLED
  gravity._val_g.resize(gravity._dim);
  gravity._val_g[0] = ml_prob.GetInputParser().get("dirgx");
  gravity._val_g[1] = ml_prob.GetInputParser().get("dirgy");
  if ( space_dim == 3 )   gravity._val_g[2] = ml_prob.GetInputParser().get("dirgz"); 

//========================
//========================
 
    currelem.Mat().zero();
    currelem.Rhs().zero(); 

    currelem.SetDofobjConnCoords(myproc,iel);
    currelem.SetMidpoint();

    currelem.ConvertElemCoordsToMappingOrd(xyz);
    currelem.TransformElemNodesToRef(ml_prob.GetMeshTwo().GetDomain(),&xyz_refbox._val_dofs[0]);    

//=======RETRIEVE the DOFS of the UNKNOWN QUANTITIES,i.e. MY EQUATION
    currelem.SetElDofsBc();
      VelOld.GetElemDofs();
    pressOld.GetElemDofs();

//======== TWO PHASE WORLD
    double rho_nd =  1.;
    double  mu_nd =  1.;
//======== TWO PHASE WORLD

    
//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
   const uint el_ngauss = ml_prob.GetQrule(currelem.GetDim()).GetGaussPointsNumber();
   
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
const double dtxJxW_g = det*ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
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
      VelOld.val_g();     //fills _val_g, needs _val_dofs
      VelOld.grad_g();     //fills _grad_g, needs _val_dofs //TODO can we see the analogies between JacVectVV_g and grad_g?
//Advection all VelOld
      for (uint idim=0; idim<space_dim; idim++) { AdvRhs_g[idim]=0.;}
          for (uint idim=0; idim<space_dim; idim++)  {
	    for (uint b=0; b<space_dim; b++) {
	    AdvRhs_g[idim] += VelOld._val_g[b]*VelOld._grad_g[idim][b]; }   }   //TODO NONLIN // grad is [ivar][idim], i.e. [u v w][x y z]
//Divergence VelOld
	  double Div_g=0.;
          for (uint idim=0; idim<space_dim; idim++)    Div_g += VelOld._grad_g[idim][idim];     //TODO NONLIN

       
#ifdef AXISYMX
      dtxJxW_g  *=yyg;  //what is the symmetry axis? I guess the x axis.
#endif

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
// TODO according to the order we should switch DIM loop and DOF loop
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
                               )
            + (1-currelem.GetBCDofFlag()[irowq])*detb*VelOld._val_dofs[irowq] //Dirichlet bc    
	;
          }
           // end filling element rhs u

  for (uint idim=0; idim<space_dim; idim++) { // filling diagonal for Dirichlet bc
          const uint irowq = i+idim*qtyzero_ndof;
          currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
        }
                                       // end filling diagonal for Dirichlet bc

//============ QTYZERO x QTYZERO dofs matrix (A matrix) ============
        for (uint j=0; j< qtyzero_ndof; j++) {
//======="COMMON SHAPE PART for QTYZERO": func and derivative, of the QTYZERO FE ORD ==========
//   (j):       dof of the SHAPE function
           double                                  phij_g       =      currgp._phi_ndsQLVB_g[qtyzero_ord][j];
           for (uint idim=0; idim<space_dim; idim++) dphijdx_g[idim] = currgp._dphidxyz_ndsQLVB_g[qtyzero_ord][j+idim*qtyzero_ndof];
//======= END "COMMON SHAPE PART for QTYZERO" ==========
  
          double Lap_g=Math::dot(&dphijdx_g[0],&dphiidx_g[0],space_dim);
	  double Adv_g=Math::dot(&VelOld._val_g[0],&dphijdx_g[0],space_dim);
          
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
       for (uint j=0; j<qtyone_ndof; j++) {
//======="COMMON SHAPE PART for QTYONE" ==================
	 const double psij_g = currgp._phi_ndsQLVB_g[qtyone_ord][j];
//======="COMMON SHAPE PART for QTYONE" - END ============

          const int jclml = j + qtyZeroToOne_DofOffset;
          for (uint idim=0; idim<space_dim; idim++) {
            uint irowq = i+idim*qtyzero_ndof;
            currelem.Mat()(irowq,jclml) +=
               currelem.GetBCDofFlag()[irowq]*                    
               dtxJxW_g*(-psij_g*dphiidx_g[idim]);   /**   (-1.)*/
	    
           } 

        }
//============ END QTYZERO x QTYONE dofs matrix (B^T matrix) ============

          if (i < qtyone_ndof) {
//======="COMMON tEST PART for QTYONE" ============
          double psii_g = currgp._phi_ndsQLVB_g[qtyone_ord][i];
//======= "COMMON tEST PART for QTYONE" - END ============
	  const uint irowl = i + qtyZeroToOne_DofOffset;
          currelem.Rhs()(irowl)=0.;  // rhs
 //             Mat()(irowl,j+space_dim*qtyzero_ndof)  += (1./dt)*dtxJxW_g*(psii_g*psij_g)*_Komp_fac/dt;  //no bc here (KOMP dp/dt=rho*div)

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
                   my_system._A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
                   my_system._b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());

    
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
  
    CurrentElem       currelem(Level,BB,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType());
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,ml_prob.GetQrule(currelem.GetDim()));
 
  
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity VelOld(currgp);
    VelOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYZERO]; //an alternative cannot exist, because it is an Unknown of This Equation
    VelOld.VectWithQtyFillBasic();
    VelOld.Allocate();

   const uint   qtyzero_ord  = VelOld._FEord;
   const uint   qtyzero_ndof = VelOld._ndof; 
//     Velocity*  vel_castqtyptr = static_cast<Velocity*>(VelOld._qtyptr); //casting for quantity-specific functions

//=========
    CurrentQuantity pressOld(currgp);
    pressOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYONE];
    pressOld.VectWithQtyFillBasic();
    pressOld.Allocate();

   const uint qtyone_ord  = pressOld._FEord;
   const uint qtyone_ndof = pressOld._ndof; 

   //order
   const uint  qtyZeroToOne_DofOffset = VelOld._ndof*VelOld._dim;
   
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp); //domain
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();
    
    //==================Quadratic domain, auxiliary, must be QUADRATIC!!! ==========
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = NVE[ ml_prob.GetMeshTwo()._geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];
  xyz_refbox.Allocate();

//=== auxiliary Operators at the boundary
//   double strainU_g['dimension']['dimension'];
//   double strainUtrDn_g['dimension'];

     currelem.Mat().zero();  
     currelem.Rhs().zero();

     currelem.SetDofobjConnCoords(myproc,iel);
     currelem.SetMidpoint();
     
     currelem.ConvertElemCoordsToMappingOrd(xyz);
     currelem.TransformElemNodesToRef(ml_prob.GetMeshTwo().GetDomain(),&xyz_refbox._val_dofs[0]);    

     currelem.SetElDofsBc();
     
       VelOld.GetElemDofs();
     pressOld.GetElemDofs();

//============ BC =======
       int press_fl = currelem.Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(VelOld,pressOld); 
//========END BC============

//==============================================================
//================== GAUSS LOOP (qp loop) ======================
//==============================================================
   const uint el_ngauss = ml_prob.GetQrule(currelem.GetDim()).GetGaussPointsNumber();
   
    for (uint qp=0; qp< el_ngauss; qp++) {
            
//======= "COMMON SHAPE PART"============================
 for (uint fe = 0; fe < QL; fe++)     {
   currgp.SetPhiElDofsFEVB_g (fe,qp);
   currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp); 
}

        const double det   = dt*currgp.JacVectBB_g(xyz);
	const double dtxJxW_g = det*ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
//=======end "COMMON SHAPE PART"===================================

//-------- pressure==============
    //"predof" VS "post gauss method"
    //pay attention to this fact: even if we are in the equation to which pressure is associated
    //(pressure is an unknown here!), we might prefer using the FUNCTION instead,
    //in order to get some values that DO NOT CHANGE with the solution
    //so i'm using the FUNCTION of an UNKNOWN, it is not an external quantity  
	  //clearly, this is ALTERNATIVE to the command
	  
   
   xyz_refbox.val_g();
      pressOld._qtyptr->Function_txyz(time,&xyz_refbox._val_g[0],&pressOld._val_g[0]);
       
	  VelOld.val_g();

//==============================================================
//========= FILLING ELEMENT MAT/RHS (i loop) ====================
//==============================================================
      for (uint i=0; i<qtyzero_ndof; i++) {

	const double phii_g = currgp._phi_ndsQLVB_g[qtyzero_ord][i];

         for (uint idim=0; idim< space_dim; idim++)    {
	   
             uint irowq=i+idim*qtyzero_ndof;
            currelem.Rhs()(irowq)  += 
         currelem.GetBCDofFlag()[irowq]*           
           dtxJxW_g*( -1.*0.*press_fl*pressOld._val_g[0]*currgp.get_normal_ptr()[idim]*phii_g );
	   
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
    
    my_system._A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
    my_system._b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());

    
  }
  // end of BDRYelement loop

    
    }
  // END BOUNDARY ******************************
  
        my_system._A[Level]->close();
        my_system._b[Level]->close();


#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << my_system.name() << ": assembled  Level " << Level
           << " with " << my_system._A[Level]->m() << " dofs" << std::endl;
#endif

    return;
}

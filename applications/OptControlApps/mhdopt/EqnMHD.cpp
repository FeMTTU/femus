
//library includes
#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "DenseVector.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"

#include "Math.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "Quantity.hpp"
#include "MultiLevelProblem.hpp"
#include "FETypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "QTYnumEnum.hpp"
#include "TimeLoop.hpp"


#include "OptLoop.hpp"
#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"

#include "FemusDefault.hpp"

using namespace femus;


void GenMatRhsMHD(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assemble_matrix) {

  SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_MHD");
  
  //========= reference values =========
  const double IRem = 1./ml_prob.GetInputParser().get("Rem");
//============================

//==== Operators @ gauss ======== 
  double         vXBe_g3D[3];
  double   vXBeXdphii_g3D[3];  
  double curlBeXdphii_g3D[3];
  
  double LapBe_g[DIMENSION];
//============================
  
   const double time =  0.; //ml_prob._timeloop._curr_time;    
    
//======= TIME - STATIONARY OR NOT =======
 const int NonStatMHD = (int) ml_prob.GetInputParser().get("NonStatMHD");
    const double dt   = 1.; //ml_prob._timeloop._timemap.get("dt");

//========= BCHandling =========
  const double penalty_val = ml_prob.GetMeshTwo().GetRuntimeMap().get("penalty_val");    
  
//======== GEOMETRICAL ELEMENT =======
  const uint space_dim =       ml_prob.GetMeshTwo().get_dim();
  const uint  mesh_ord = (int) ml_prob.GetMeshTwo().GetRuntimeMap().get("mesh_ord");
  const uint  meshql   = (int) ml_prob.GetMeshTwo().GetRuntimeMap().get("meshql");   //======== ELEMENT MAPPING =======
  
        my_system._A[Level]->zero();
        my_system._b[Level]->zero();

  {//BEGIN VOLUME
//======================
//======================
   const uint mesh_vb = VV;
   
    const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*ml_prob.GetMeshTwo()._iproc+Level+1];
    const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*ml_prob.GetMeshTwo()._iproc+Level];
 
    for (uint iel=0; iel < (nel_e - nel_b); iel++) {
   
    CurrentElem       currelem(Level,VV,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType());
   
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,ml_prob.GetQrule(currelem.GetDim()));
    
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity bhomOld(currgp);
    bhomOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYZERO];
    bhomOld.VectWithQtyFillBasic();
    bhomOld.Allocate();

    CurrentQuantity LagMultOld(currgp);
    LagMultOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYONE];
    LagMultOld.VectWithQtyFillBasic();
    LagMultOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

    
//========= tEST AND SHAPE FOR QTYZERO AND QTYONE =================
//QTYZERO SHAPE: shape of the first Unknown
    CurrentQuantity Phij(currgp); //TODO this is another Vect that doesnt have an associated quantity still
    Phij._dim      = 1;                                                         //scalar!
    Phij._FEord    = bhomOld._FEord;
    Phij._ndof     = currelem.GetElemType(Phij._FEord)->GetNDofs(); 
    Phij.Allocate();
        
//QTYZERO tEST:  test of the first Unknown    
    CurrentQuantity Phii(currgp);
    Phii._dim      = 1;
    Phii._FEord    = bhomOld._FEord;
    Phii._ndof     = currelem.GetElemType(Phii._FEord)->GetNDofs();
    Phii.Allocate();
    
//QTYONE SHAPE: Shape of the second Unknown
    CurrentQuantity Psij(currgp);
    Psij._dim      = 1;
    Psij._FEord    = LagMultOld._FEord;
    Psij._ndof     = currelem.GetElemType(Psij._FEord)->GetNDofs(); 
    Psij.Allocate();
    
//QTYONE tEST: test of the second Unknown
    CurrentQuantity Psii(currgp);
    Psii._dim      = 1;
    Psii._FEord    = LagMultOld._FEord;
    Psii._ndof     = currelem.GetElemType(Psii._FEord)->GetNDofs();
    Psii.Allocate();

//========= END tEST AND SHAPE FOR QTYZERO AND QTYONE =================

    
//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();
    
    //================== Quadratic domain, auxiliary
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = NVE[ ml_prob.GetMeshTwo()._geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];
  xyz_refbox.Allocate();

    CurrentQuantity Bext(currgp);
    Bext._qtyptr     = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt");
    Bext.VectWithQtyFillBasic();
    Bext.Allocate();

#if VELOCITY_QTY==1
    CurrentQuantity Vel(currgp);
    Vel._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity"); 
    Vel.VectWithQtyFillBasic();
    Vel.Allocate();
#endif  

//=========END EXTERNAL QUANTITIES (couplings) =====
//======================
//======================
    
    
    currelem.Mat().zero();
    currelem.Rhs().zero();

     currelem.SetDofobjConnCoords(ml_prob.GetMeshTwo()._iproc,iel);
     currelem.SetMidpoint();
     
    currelem.ConvertElemCoordsToMappingOrd(xyz);
    currelem.TransformElemNodesToRef(ml_prob.GetMeshTwo().GetDomain(),&xyz_refbox._val_dofs[0]);    

    currelem.SetElDofsBc();
    
       bhomOld.GetElemDofs();
    LagMultOld.GetElemDofs();
  
#if BMAG_QTY==1
    if ( Bext._eqnptr != NULL )  Bext.GetElemDofs();
    else                         Bext._qtyptr->FunctionDof(Bext,time,&xyz_refbox._val_dofs[0]);
#endif
#if VELOCITY_QTY==1
    if ( Vel._eqnptr != NULL )  Vel.GetElemDofs();      //----- for Advection MAT & RHS
    else                        Vel._qtyptr->FunctionDof(Vel,time,&xyz_refbox._val_dofs[0]);
#endif
    
   const uint el_ngauss = ml_prob.GetQrule(currelem.GetDim()).GetGaussPointsNumber();

    for (uint qp = 0; qp < el_ngauss; qp++) {

//======= here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     { 
  currgp.SetPhiElDofsFEVB_g (fe,qp);
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
}

 const double      det = dt*currgp.JacVectVV_g(xyz);   //InvJac: is unique!
 const double dtxJxW_g = det*ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
 const double     detb = det/el_ngauss;

for (uint fe = 0; fe < QL; fe++)     {    
  currgp.SetDPhiDxyzElDofsFEVB_g (fe,qp);
  currgp.ExtendDphiDxyzElDofsFEVB_g (fe); 
}
//======= end of the "COMMON SHAPE PART"==================

      bhomOld.val_g();         //---for Time

#if VELOCITY_QTY==1
      Vel.val_g();            //---- for Advection MAT & RHS
#endif
#if BMAG_QTY==1
      Bext.val_g();          //----- for Advection RHS
      Bext.grad_g();          //----- for Laplacian RHS
      Bext.curl_g();          //----- for Curl Curl RHS         //THE EXTENSION of the DOFS to 3D is done INSIDE!!
#endif

#if (VELOCITY_QTY==1) && (BMAG_QTY==1) //in this case we have two couplings with external quantities
       Math::extend(&Vel._val_g[0],&Vel._val_g3D[0],space_dim);                    //----- for Advection RHS
       Math::extend(&Bext._val_g[0],&Bext._val_g3D[0],space_dim);                    //----- for Advection RHS
       Math::cross(&Vel._val_g3D[0],&Bext._val_g3D[0],vXBe_g3D);          //----- for Advection RHS
#endif

//================================
//========= FILLING ELEMENT MAT/RHS
//=================================
       for (uint i=0; i < bhomOld._ndof; i++)     {
//============ preparation for (i) (Phii) ============
           Phii._val_g[0]       =      currgp._phi_ndsQLVB_g[Phii._FEord][i];          /*const double  phii_g*/ 
        for (uint idim=0; idim<space_dim; idim++)
	  Phii._grad_g[0][idim] = currgp._dphidxyz_ndsQLVB_g[Phii._FEord][i+idim*Phii._ndof];   /*dphiidx_g*/
//=== from this point on, the dependency on /*(i)*/, the ROW INDEX, of the quantities is given by Phii  
//===  /*(i)*/ is ONLY used for the POSITIONS, later
	  
    //--------- LAPLACIAN: Operator, RHS: grad Be . grad phi ---------------------
	 for (uint idim=0; idim<space_dim/*Bext._dim*/; idim++) LapBe_g[idim]=0.;

          for (uint idim=0; idim<space_dim/*Bext._dim*/; idim++) {
 	    for (uint jdim=0; jdim<space_dim; jdim++) {
              LapBe_g[idim] += Bext._grad_g[idim][jdim]*Phii._grad_g[0][jdim];
	    }
 	  }

             Math::extend(&Phii._grad_g[0][0],&Phii._grad_g3D[0][0],space_dim);

    //--------- CURL CURL: Operator, RHS: curl Be . curl phi -------------------
             Math::cross(&Bext._curl_g3D[0],&Phii._grad_g3D[0][0],curlBeXdphii_g3D);

    //--------- ADVECTION: Operator, RHS: v x Be . curl phi -------------------
             Math::cross(      vXBe_g3D,&Phii._grad_g3D[0][0],  vXBeXdphii_g3D);       // _utils.cross(vXBe_g3D,dphiidx_g3D,vXBeXdphii_g3D);
//============end preparation for (i) ============
	     
         for (uint idim=0; idim<space_dim/*bhomOld._dim*/; idim++) {
            const uint irowq = i + idim*bhomOld._ndof;
            
            currelem.Rhs()(irowq) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                           NonStatMHD*bhomOld._val_g[idim]*Phii._val_g[0]/dt //time
                           - LAP_MHD*IRem*LapBe_g[idim]
                            - (1-LAP_MHD)*IRem*curlBeXdphii_g3D[idim]              //phii of bhomOld /*CurlCurl(RHS,vb,Phij,Phii,idim,idimp1);*/
                           + ADV_MHD*           vXBeXdphii_g3D[idim]               //phii of bhomOld
                         )
                         + (1-currelem.GetBCDofFlag()[irowq])*detb*bhomOld._val_dofs[irowq]; //Dirichlet bc
	   }


        for (uint idim=0; idim<space_dim/*bhomOld._dim*/; idim++) { // filling diagonal block for Dirichlet bc
          const uint irowq = i + idim*bhomOld._ndof;
          currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
        }
                                           // end filling diagonal for Dirichlet bc
	 
        for (uint j=0; j<bhomOld._ndof; j++) {  // A element matrix
//============ preparation for (j) (Phij) ============
            Phij._val_g[0]        =      currgp._phi_ndsQLVB_g[Phij._FEord][j]; 
          for (uint idim=0; idim<space_dim; idim++)
	    Phij._grad_g[0][idim] = currgp._dphidxyz_ndsQLVB_g[Phij._FEord][j+idim*Phij._ndof];   //real shape  /*dphijdx_g[idim]*/
//=== from this point on, the dependency on /*(j)*/, the COLUMN INDEX, of the quantities is given by Phij  
//===  /*(j)*/ is ONLY used for the POSITIONS, later

    //--------- LAPLACIAN: Operator, MAT: grad b . grad phi ---------------------
          double     Lap_g = Math::dot(&Phij._grad_g[0][0],&Phii._grad_g[0][0],space_dim);   /*(i,j)*/  //part independent of idim
    //--------- ADVECTION: Operator, MAT: v x b . curl phi -------------------
          double Advphii_g = Math::dot( &Vel._val_g[0],&Phii._grad_g[0][0],space_dim);  /*(i)*/ //part independent of idim //TODO what about putting it OUTSIDE?
     
          for (uint idim=0; idim<space_dim/*bhomOld._dim*/; idim++) { //filled in as 1-2-3 // 5-6-4 // 9-7-8
            int irowq=i+idim*bhomOld._ndof;   //idim gives the row index  //test of bhomOld
            // diagonal blocks [1-5-9]  idim = row index = column index  //shape of bhomOld
            currelem.Mat()(irowq,j+idim*bhomOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                NonStatMHD*Phij._val_g[0]*Phii._val_g[0]/dt 
                 + LAP_MHD*IRem*(Lap_g) 
                 + (1-LAP_MHD)*IRem*(           Lap_g -  Phij._grad_g[0][idim]*Phii._grad_g[0][idim] )
                 - ADV_MHD * Phij._val_g[0]* (Advphii_g -   Vel._val_g[idim]    *Phii._grad_g[0][idim] )   //TODO Phij here does not depend on idim, but it depends on j
               );
            // block +1 [2-6-7]
            int idimp1=(idim+1)%space_dim;
	    //idim: component of the SHAPE
	    //idimp1: component of the tEST
            currelem.Mat()(irowq,j+idimp1*bhomOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                  + (1-LAP_MHD)*IRem*(        - Phij._grad_g[0][idim]*Phii._grad_g[0][idimp1] )       /*(i,j)*/    /*CurlCurl(MAT,vb,Phij,Phii,idim,idimp1);*/
                  - ADV_MHD * Phij._val_g[0]* ( -  Vel._val_g[idim]    *Phii._grad_g[0][idimp1] )       /*(i,j)*/    /*AdvCurl(MAT,Vel,Phij,Phii,idim,idimp1)*/
               );

#if (DIMENSION==3)
            // block +2 [3-4-8] 
	    int idimp2=(idim+2)%space_dim;//idimp2 column index
            currelem.Mat()(irowq,j+idimp2*bhomOld._ndof)
            += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                  + (1-LAP_MHD)*IRem*(       - Phij._grad_g[0][idim]*Phii._grad_g[0][idimp2] )               /*(i,j)*/
                  - ADV_MHD * Phij._val_g[0]* ( -  Vel._val_g[idim]    *Phii._grad_g[0][idimp2] )               /*(i,j)*/
               );
#endif
          }
                    
        } 
                                      // end A element matrix
				      
       for (uint j = 0; j < LagMultOld._ndof; j++) {// B^T element matrix ( p*div(v) )
          Psij._val_g[0] = currgp._phi_ndsQLVB_g[Psij._FEord][j];   /*const double psij_g */
          const int jclml= j+/*bhomOld._dim*/space_dim*bhomOld._ndof;
          for (uint idim=0; idim<space_dim; idim++) {  //bhomOld._dim==spacedimension
            uint irowq=i+idim*bhomOld._ndof;
            currelem.Mat()(irowq,jclml) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*( - Psij._val_g[0]*Phii._grad_g[0][idim] );   /*(i,j)*/ /*PDiv(MAT,Psij,Phii,idim)*/
           }
        }
                                     // end B^T element matrix
 //============ begin linear rows
          if (i < LagMultOld._ndof) {
           Psii._val_g[0]  = currgp._phi_ndsQLVB_g[Psii._FEord][i];   /*double psii_g*/
	  const uint irowl=i+space_dim*bhomOld._ndof;
          currelem.Rhs()(irowl)=0.;

          for (uint j=0; j<bhomOld._ndof; j++) { // B element matrix q*div(u)
            for (uint idim=0; idim<space_dim; idim++) Phij._grad_g[0][idim] =  currgp._dphidxyz_ndsQLVB_g[Phij._FEord][j+idim*Phij._ndof];
            for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*bhomOld._ndof) += - dtxJxW_g*Psii._val_g[0]*Phij._grad_g[0][idim];  /*(i,j)*/  /*PDiv(MAT,Psii,Phij,idim)*/
                }

        }
//============ end linear rows
      }
    //end filling element mat/rhs (i loop)
    }
    // end element gaussian integration loop
    
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
  
//======================
//======================
    const uint mesh_vb = BB;
    
    const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*ml_prob.GetMeshTwo()._iproc+Level+1];
    const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*ml_prob.GetMeshTwo()._iproc+Level];
  
  for (uint iel=0;iel < (nel_e - nel_b) ; iel++) {
   
    CurrentElem       currelem(Level,BB,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType());
   
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,ml_prob.GetQrule(currelem.GetDim()));
    
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity bhomOld(currgp);
    bhomOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYZERO];
    bhomOld.VectWithQtyFillBasic();
    bhomOld.Allocate();

    CurrentQuantity LagMultOld(currgp);
    LagMultOld._qtyptr   = my_system.GetUnknownQuantitiesVector()[QTYONE];
    LagMultOld.VectWithQtyFillBasic();
    LagMultOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

    
//========= tEST AND SHAPE FOR QTYZERO AND QTYONE =================
//QTYZERO SHAPE: shape of the first Unknown
    CurrentQuantity Phij(currgp); //TODO this is another Vect that doesnt have an associated quantity still
    Phij._dim      = 1;                                                         //scalar!
    Phij._FEord    = bhomOld._FEord;
    Phij._ndof     = currelem.GetElemType(Phij._FEord)->GetNDofs(); 
    Phij.Allocate();
        
//QTYZERO tEST:  test of the first Unknown    
    CurrentQuantity Phii(currgp);
    Phii._dim      = 1;
    Phii._FEord    = bhomOld._FEord;
    Phii._ndof     = currelem.GetElemType(Phii._FEord)->GetNDofs();
    Phii.Allocate();
    
//QTYONE SHAPE: Shape of the second Unknown
    CurrentQuantity Psij(currgp);
    Psij._dim      = 1;
    Psij._FEord    = LagMultOld._FEord;
    Psij._ndof     = currelem.GetElemType(Psij._FEord)->GetNDofs(); 
    Psij.Allocate();
    
//QTYONE tEST: test of the second Unknown
    CurrentQuantity Psii(currgp);
    Psii._dim      = 1;
    Psii._FEord    = LagMultOld._FEord;
    Psii._ndof     = currelem.GetElemType(Psii._FEord)->GetNDofs();
    Psii.Allocate();

//========= END tEST AND SHAPE FOR QTYZERO AND QTYONE =================

    
//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();
    
    //================== Quadratic domain, auxiliary
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = NVE[ ml_prob.GetMeshTwo()._geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];
  xyz_refbox.Allocate();

    CurrentQuantity Bext(currgp);
    Bext._qtyptr     = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldExt");
    Bext.VectWithQtyFillBasic();
    Bext.Allocate();

#if VELOCITY_QTY==1
    CurrentQuantity Vel(currgp);
    Vel._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity"); 
    Vel.VectWithQtyFillBasic();
    Vel.Allocate();
#endif  

//=========END EXTERNAL QUANTITIES (couplings) =====
  
  //=======Auxiliary Operators for the Boundary Integrals
    double velXBext_g3D[3];  //Operations at the boundary
      double normal_g3D[3];  //Operations at the boundary
  double velXBextXn_g3D[3];  //Operations at the boundary
    
//======================
//======================
 

     currelem.Mat().zero();
     currelem.Rhs().zero(); 
     
     currelem.SetDofobjConnCoords(ml_prob.GetMeshTwo()._iproc,iel);
     currelem.SetMidpoint();

     currelem.ConvertElemCoordsToMappingOrd(xyz);
     currelem.TransformElemNodesToRef(ml_prob.GetMeshTwo().GetDomain(),&xyz_refbox._val_dofs[0]);    
   
     currelem.SetElDofsBc();
     
        bhomOld.GetElemDofs();
     LagMultOld.GetElemDofs();

//============ BC =======
       int press_fl = currelem.Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(bhomOld,LagMultOld); 
//========END BC=========    
    
//========== EXTERNAL DOFS ===   
#if BMAG_QTY==1
    if ( Bext._eqnptr != NULL )   Bext.GetElemDofs();
    else                          Bext._qtyptr->FunctionDof(Bext,time,&xyz_refbox._val_dofs[0]);
#endif
#if VELOCITY_QTY==1
    if ( Vel._eqnptr != NULL )  Vel.GetElemDofs();
    else                        Vel._qtyptr->FunctionDof(Vel,time,&xyz_refbox._val_dofs[0]);
#endif
    

   const uint el_ngauss = ml_prob.GetQrule(currelem.GetDim()).GetGaussPointsNumber();
   
    for (uint qp=0; qp< el_ngauss; qp++) {
      
//=======here starts the "COMMON SHAPE PART"============================
      for (uint fe = 0; fe < QL; fe++)     {
	currgp.SetPhiElDofsFEVB_g (fe,qp);  
        currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp); 
      }
      
        const double det      = dt * currgp.JacVectBB_g(xyz);
	const double dtxJxW_g = det * ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
//=======end of the "COMMON SHAPE PART"===================================

//---------lagmult
	    //"post gauss method"
   xyz_refbox.val_g(); // val_g(vb,xyz);   //CHECK the QUADRATICS!!!!!!!!!
      LagMultOld._qtyptr->Function_txyz(time,&xyz_refbox._val_g[0],&LagMultOld._val_g[0]);  //check that you have ZERO here
      
#if VELOCITY_QTY==1
     Vel.val_g();
#endif
#if BMAG_QTY==1
     Bext.val_g();
#endif
#if (VELOCITY_QTY==1) && (BMAG_QTY==1)
          Math::extend( &Vel._val_g[0],&Vel._val_g3D[0],space_dim);
	  Math::extend(&Bext._val_g[0],&Bext._val_g3D[0],space_dim);
	  Math::cross(&Vel._val_g3D[0],&Bext._val_g3D[0],velXBext_g3D);
  	  Math::extend(currgp.get_normal_ptr(),normal_g3D,space_dim);
	  Math::cross(velXBext_g3D,normal_g3D,velXBextXn_g3D);
#endif


      // BDRYelement rhs filling
      for (uint i=0; i< bhomOld._ndof; i++) {  //loop over the BDRYelement dofs

	 Phii._val_g[0] = currgp._phi_ndsQLVB_g[Phii._FEord][i];

          // rhs
         for (uint idim=0; idim< space_dim; idim++)    {
             uint irowq=i+idim*bhomOld._ndof;
            currelem.Rhs()(irowq)  +=
            currelem.GetBCDofFlag()[irowq]*
            dtxJxW_g*(    -1.*press_fl*LagMultOld._val_g[0]*currgp.get_normal_ptr()[idim]*Phii._val_g[0]
                           + velXBextXn_g3D[idim]*Phii._val_g[0]  //advection  //this type of advection has also a boundary integral, if integration by parts is performed in similar manner
                             );

          }
          //end of idim loop

     }
      // end of i loop

    } 
    // end BDRYelement gaussian integration loop

   my_system._A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
   my_system._b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());

 
  }
  // end of BDRYelement loop  
  
    }  
    
// END BOUNDARY  // *****************************************************************

        my_system._A[Level]->close();
        my_system._b[Level]->close();
    
#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << my_system.name() << ": assembled  Level " << Level
           << " with " << my_system._A[Level]->m() << " dofs " << std::endl;
#endif     


  return;
}


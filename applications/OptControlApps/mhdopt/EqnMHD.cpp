#include "EqnMHD.hpp"

//library includes
#include "DenseMatrix.hpp"
#include "SparseMatrix.hpp"
#include "DenseVector.hpp"
#include "NumericVector.hpp"
#include "LinearEquationSolver.hpp"

#include "Math.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "GeomEl.hpp"
#include "Quantity.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "FETypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "QTYnumEnum.hpp"
#include "TimeLoop.hpp"


#include "Opt_conf.hpp"
#include "EqnNS.hpp"
#include "EqnMHDCONT.hpp"

#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"

#include "FemusDefault.hpp"


namespace femus {

///Constructor
  EqnMHD::EqnMHD(  std::vector<Quantity*> int_map_in,
	           MultiLevelProblemTwo& mg_equations_map_in,
                   std::string eqname_in,
                   std::string varname_in):
      SystemTwo(int_map_in,mg_equations_map_in,eqname_in,varname_in)
                  {

//====== VARNAMES of the equation
   _var_names[0]="Bx";
   _var_names[1]="By";
   _var_names[DIMENSION]="Bp";
#if (DIMENSION==3)
    _var_names[2]="Bz";
#endif

//=========  REFVALUES of the Unknown quantities of the equation
   _refvalue[0] = _QtyInternalVector[0]->_refvalue[0];/* Bref;*/
   _refvalue[1] = _QtyInternalVector[0]->_refvalue[1]; /*Bref*/
#if (DIMENSION==3)
    _refvalue[2]= _QtyInternalVector[0]->_refvalue[2];/*Bref*/
#endif
   _refvalue[DIMENSION]= _QtyInternalVector[1]->_refvalue[0];/*sigmaref = Uref*Bref*/
   

  //===== solver  
  for(uint l=0;l<_NoLevels;l++)  _solver[l]->set_solver_type(SOLVERMHD);

    //////DIR PENALTY ///////////
   _Dir_pen_fl = MHD_DIR_PENALTY;

   
  }

  
///============ DESTRUCTOR    
   EqnMHD::~EqnMHD() {     }
    
  


/// This function assembles the matrix and the rhs:

  void EqnMHD::GenMatRhs(const uint Level)  {

  //========= reference values =========
  const double IRem = 1./_phys.get("Rem");
//============================

//==== Operators @ gauss ======== 
  double         vXBe_g3D[3];
  double   vXBeXdphii_g3D[3];  
  double curlBeXdphii_g3D[3];
  
  double LapBe_g[DIMENSION];
//============================
  
   const double time =  0.; //_eqnmap._timeloop._curr_time;    
    
//======= TIME - STATIONARY OR NOT =======
 const int NonStatMHD = (int) _phys.get("NonStatMHD");
    const double dt   = 1.; //_eqnmap._timeloop._timemap.get("dt");

//========= BCHandling =========
  const double penalty_val = _mesh.GetRuntimeMap().get("penalty_val");    
  
//======== GEOMETRICAL ELEMENT =======
  const uint space_dim =       _mesh.get_dim();
  const uint  mesh_ord = (int) _mesh.GetRuntimeMap().get("mesh_ord");
  const uint  meshql   = (int) _mesh.GetRuntimeMap().get("meshql");   //======== ELEMENT MAPPING =======
     
  {//BEGIN VOLUME
//======================
//======================
   const uint mesh_vb = VV;
   
   CurrentElem       currelem(VV,this,_mesh,_eqnmap._elem_type);
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,_eqnmap, _mesh.get_dim());
    
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity bhomOld(currgp);
    bhomOld._qtyptr   = _QtyInternalVector[QTYZERO];
    bhomOld.VectWithQtyFillBasic();
    bhomOld.Allocate();

    CurrentQuantity LagMultOld(currgp);
    LagMultOld._qtyptr   = _QtyInternalVector[QTYONE];
    LagMultOld.VectWithQtyFillBasic();
    LagMultOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

    
//========= tEST AND SHAPE FOR QTYZERO AND QTYONE =================
//QTYZERO SHAPE: shape of the first Unknown
    CurrentQuantity Phij(currgp); //TODO this is another Vect that doesnt have an associated quantity still
    Phij._dim      = 1;                                                         //scalar!
    Phij._FEord    = bhomOld._FEord;
    Phij._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][Phij._FEord]->GetNDofs(); 
    Phij.Allocate();
        
//QTYZERO tEST:  test of the first Unknown    
    CurrentQuantity Phii(currgp);
    Phii._dim      = 1;
    Phii._FEord    = bhomOld._FEord;
    Phii._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][Phii._FEord]->GetNDofs();
    Phii.Allocate();
    
//QTYONE SHAPE: Shape of the second Unknown
    CurrentQuantity Psij(currgp);
    Psij._dim      = 1;
    Psij._FEord    = LagMultOld._FEord;
    Psij._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][Psij._FEord]->GetNDofs(); 
    Psij.Allocate();
    
//QTYONE tEST: test of the second Unknown
    CurrentQuantity Psii(currgp);
    Psii._dim      = 1;
    Psii._FEord    = LagMultOld._FEord;
    Psii._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][Psii._FEord]->GetNDofs();
    Psii.Allocate();

//========= END tEST AND SHAPE FOR QTYZERO AND QTYONE =================

    
//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();
    
    //================== Quadratic domain, auxiliary
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = _mesh.GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
  xyz_refbox.Allocate();

    CurrentQuantity Bext(currgp);
    Bext._qtyptr     = _eqnmap._qtymap.get_qty("Qty_MagnFieldExt");
    Bext.VectWithQtyFillBasic();
    Bext.Allocate();

#if VELOCITY_QTY==1
    CurrentQuantity Vel(currgp);
    Vel._qtyptr      = _eqnmap._qtymap.get_qty("Qty_Velocity"); 
    Vel.VectWithQtyFillBasic();
    Vel.Allocate();
#endif  

//=========END EXTERNAL QUANTITIES (couplings) =====
    
    const uint nel_e = _mesh._off_el[mesh_vb][_NoLevels*_iproc+Level+1];
    const uint nel_b = _mesh._off_el[mesh_vb][_NoLevels*_iproc+Level];

//======================
//======================
    
 
    for (uint iel=0; iel < (nel_e - nel_b); iel++) {

    currelem.Mat().zero();
    currelem.Rhs().zero();

     currelem.set_el_nod_conn_lev_subd(Level,_iproc,iel);
     currelem.SetMidpoint();
     
     currelem.ConvertElemCoordsToMappingOrd(xyz);
     _mesh.TransformElemNodesToRef(currelem.GetDim(),currelem.GetNodeCoords(),&xyz_refbox._val_dofs[0]);    

    currelem.SetElDofsBc(Level);
    
       bhomOld.GetElemDofs(Level);
    LagMultOld.GetElemDofs(Level);
  
      if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),bhomOld._FEord,currelem.GetBCDofFlag());  //only the Quadratic Part is modified!



#if BMAG_QTY==1
    if ( Bext._eqnptr != NULL )  Bext.GetElemDofs(Level);
    else                         Bext._qtyptr->FunctionDof(Bext,time,&xyz_refbox._val_dofs[0]);
#endif
#if VELOCITY_QTY==1
    if ( Vel._eqnptr != NULL )  Vel.GetElemDofs(Level);      //----- for Advection MAT & RHS
    else                        Vel._qtyptr->FunctionDof(Vel,time,&xyz_refbox._val_dofs[0]);
#endif
    
   const uint el_ngauss = _eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();

    for (uint qp = 0; qp < el_ngauss; qp++) {

//======= here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     { 
  currgp.SetPhiElDofsFEVB_g (fe,qp);
  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
}

 const double      det = dt*currgp.JacVectVV_g(xyz);   //InvJac: is unique!
 const double dtxJxW_g = det*_eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
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
                      _A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
                      _b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices()); 

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
   
   CurrentElem       currelem(BB,this,_mesh,_eqnmap._elem_type);
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,_eqnmap, _mesh.get_dim());
    
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity bhomOld(currgp);
    bhomOld._qtyptr   = _QtyInternalVector[QTYZERO];
    bhomOld.VectWithQtyFillBasic();
    bhomOld.Allocate();

    CurrentQuantity LagMultOld(currgp);
    LagMultOld._qtyptr   = _QtyInternalVector[QTYONE];
    LagMultOld.VectWithQtyFillBasic();
    LagMultOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

    
//========= tEST AND SHAPE FOR QTYZERO AND QTYONE =================
//QTYZERO SHAPE: shape of the first Unknown
    CurrentQuantity Phij(currgp); //TODO this is another Vect that doesnt have an associated quantity still
    Phij._dim      = 1;                                                         //scalar!
    Phij._FEord    = bhomOld._FEord;
    Phij._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][Phij._FEord]->GetNDofs(); 
    Phij.Allocate();
        
//QTYZERO tEST:  test of the first Unknown    
    CurrentQuantity Phii(currgp);
    Phii._dim      = 1;
    Phii._FEord    = bhomOld._FEord;
    Phii._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][Phii._FEord]->GetNDofs();
    Phii.Allocate();
    
//QTYONE SHAPE: Shape of the second Unknown
    CurrentQuantity Psij(currgp);
    Psij._dim      = 1;
    Psij._FEord    = LagMultOld._FEord;
    Psij._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][Psij._FEord]->GetNDofs(); 
    Psij.Allocate();
    
//QTYONE tEST: test of the second Unknown
    CurrentQuantity Psii(currgp);
    Psii._dim      = 1;
    Psii._FEord    = LagMultOld._FEord;
    Psii._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][Psii._FEord]->GetNDofs();
    Psii.Allocate();

//========= END tEST AND SHAPE FOR QTYZERO AND QTYONE =================

    
//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = _eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();
    
    //================== Quadratic domain, auxiliary
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = _mesh.GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
  xyz_refbox.Allocate();

    CurrentQuantity Bext(currgp);
    Bext._qtyptr     = _eqnmap._qtymap.get_qty("Qty_MagnFieldExt");
    Bext.VectWithQtyFillBasic();
    Bext.Allocate();

#if VELOCITY_QTY==1
    CurrentQuantity Vel(currgp);
    Vel._qtyptr      = _eqnmap._qtymap.get_qty("Qty_Velocity"); 
    Vel.VectWithQtyFillBasic();
    Vel.Allocate();
#endif  

//=========END EXTERNAL QUANTITIES (couplings) =====
    
    const uint nel_e = _mesh._off_el[mesh_vb][_NoLevels*_iproc+Level+1];
    const uint nel_b = _mesh._off_el[mesh_vb][_NoLevels*_iproc+Level];

//======================
//======================
  
  //=======Auxiliary Operators for the Boundary Integrals
    double velXBext_g3D[3];  //Operations at the boundary
      double normal_g3D[3];  //Operations at the boundary
  double velXBextXn_g3D[3];  //Operations at the boundary
 
  
  for (uint iel=0;iel < (nel_e - nel_b) ; iel++) {

     currelem.Mat().zero();
     currelem.Rhs().zero(); 
     
     currelem.set_el_nod_conn_lev_subd(Level,_iproc,iel);
     currelem.SetMidpoint();

     currelem.ConvertElemCoordsToMappingOrd(xyz);
    _mesh.TransformElemNodesToRef(currelem.GetDim(),currelem.GetNodeCoords(),&xyz_refbox._val_dofs[0]);    
   
     currelem.SetElDofsBc(Level);
     
        bhomOld.GetElemDofs(Level);
     LagMultOld.GetElemDofs(Level);

    if (_Dir_pen_fl == 1) Bc_ConvertToDirichletPenalty(currelem.GetDim(),bhomOld._FEord,currelem.GetBCDofFlag()); //only the Quadratic Part is modified! /*OK DIR_PEN*/
       

//============ BC =======
       int     el_flag[NT] = {0,0}; //normal and tangential flag
#if DIMENSION==2
       double el_value[N1T1] = {0.,0.};  //1 normal and 1 tangential
#elif DIMENSION==3
       double el_value[N1T3] = {0.,0.,0.,0.}; //1 normal and 3 tangential
#endif
       double  dbl_pen[NT] = {0.,0.};  //normal and tangential penalty value
   
       Bc_GetElFlagValLevSubd(Level,_iproc,iel,el_flag,el_value);

if (_Dir_pen_fl == 1)  { 
       if (el_flag[NN] == 1) {   dbl_pen[NN]=penalty_val; } //normal dirichlet
       if (el_flag[TT] == 1) {   dbl_pen[TT]=penalty_val; } //tangential dirichlet
   }

//checks
//TODO here i should check that the nodal bc dirichlet i put correspond to the element NT flags

       uint press_fl=0;
       Bc_ComputeElementBoundaryFlagsFromNodalFlagsForPressure(currelem.GetBCDofFlag(),bhomOld,LagMultOld,press_fl); //compute the PRESSURE FLAG with the PRESSURE nodal bc flags
 //only the LINEAR PART is USED!!
       
// // //   if ( (1-el_flag[NN]) != press_fl)  {std::cout << _eqname << ": sthg. wrong with press elflags" << std::endl;abort();}
//========END BC=========    
    
//========== EXTERNAL DOFS ===   
#if BMAG_QTY==1
    if ( Bext._eqnptr != NULL )   Bext.GetElemDofs(Level);
    else                          Bext._qtyptr->FunctionDof(Bext,time,&xyz_refbox._val_dofs[0]);
#endif
#if VELOCITY_QTY==1
    if ( Vel._eqnptr != NULL )  Vel.GetElemDofs(Level);
    else                        Vel._qtyptr->FunctionDof(Vel,time,&xyz_refbox._val_dofs[0]);
#endif
    

   const uint el_ngauss = _eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();
   
    for (uint qp=0; qp< el_ngauss; qp++) {
      
//=======here starts the "COMMON SHAPE PART"============================
      for (uint fe = 0; fe < QL; fe++)     {
	currgp.SetPhiElDofsFEVB_g (fe,qp);  
        currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);   }
      
        const double det      = dt * currgp.JacVectBB_g(xyz);
	const double dtxJxW_g = det * _eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);
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
            dtxJxW_g*(    -1.*(1-el_flag[NN])*LagMultOld._val_g[0]*currgp.get_normal_ptr()[idim]*Phii._val_g[0]   //AAA multiplying int times uint!!!   /*PDiv(RHS,(LagMultOld(Psiii),Phii,idim)*/
// // //                            + Neum_fl*IRe*strainUtrDn_g[idim]*Phii._val_g[0]          // TODO: \tau \vect{n} on the boundary for the Laplacian
                           + velXBextXn_g3D[idim]*Phii._val_g[0]  //advection  //this type of advection has also a boundary integral, if integration by parts is performed in similar manner
                             )

                                //projection over the physical (x,y,z)
      + _Dir_pen_fl *dtxJxW_g*Phii._val_g[0]*(dbl_pen[NN]*el_value[0]*currgp.get_normal_ptr()[idim] 
                                    + dbl_pen[TT]*el_value[1]*currgp.get_tangent_ptr()[0][idim]  // VelOld._val_g[idim] instead of el_value...
               #if DIMENSION==3
		                    + dbl_pen[TT]*el_value[1]*currgp.get_tangent_ptr()[1][idim]    
               #endif   
          )    
 

	    ;   
	    
//====================
if (_Dir_pen_fl == 1) {  //much faster than multiplying by _Dir_pen_fl=0 , and much better than removing the code with the #ifdef //  #if (NS_DIR_PENALTY==1)  
	   for (uint jdim=0; jdim< space_dim; jdim++)    {

	   for (uint j=0; j<bhomOld._ndof; j++) {
         Phij._val_g[0] = currgp._phi_ndsQLVB_g[Phij._FEord][j];      /*const double phij_g*/

     currelem.Mat()(irowq,j+jdim*bhomOld._ndof) +=                //projection over the physical (x,y,z) 
         + /*_Dir_pen_fl**/dtxJxW_g*Phii._val_g[0]*Phij._val_g[0]*(dbl_pen[NN]*currgp.get_normal_ptr()[jdim]*currgp.get_normal_ptr()[idim]   //the PENALTY is BY ELEMENT, but the (n,t) is BY GAUSS because we cannot compute now a nodal normal
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

#if THIN_WALL==1
    for (uint j=0; j<el_ndof_q; j++) {
	const double phij_g=_phi_nds_gVB_QQ[vb][j];

//b is in the y direction
           uint irowq=i+1*el_ndof_q;
           uint jrowq=j+1*el_ndof_q;

         currelem.Mat()(irowq,jrowq) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
                                  -CONDRATIO*phij_g*Phii._val_g[0]/*phii_g*/ //-
                          );

        }
#endif

      }
      // end of i loop

    } 
    // end BDRYelement gaussian integration loop

   _A[Level]->add_matrix(currelem.Mat(),currelem.GetDofIndices());
   _b[Level]->add_vector(currelem.Rhs(),currelem.GetDofIndices());

 
  }
  // end of BDRYelement loop  
  
    }  
    
// END BOUNDARY  // *****************************************************************
     
    
#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << _eqname << ": assembled  Level " << Level
           << " with " << _A[Level]->m() << " dofs " << std::endl;
#endif     


  return;
}



} //end namespace femus



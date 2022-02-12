
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
#include "GeomElTypeEnum.hpp"
#include "NormTangEnum.hpp"
#include "TimeLoop.hpp"


#include "OptLoop.hpp"
#include "CurrentElem.hpp"

#include "FemusDefault.hpp"

using namespace femus;


void GenMatRhsMHD(MultiLevelProblem &ml_prob) {

  SystemTwo & my_system = ml_prob.get_system<SystemTwo>("Eqn_MHD");
  const unsigned Level = my_system.GetLevelToAssemble();
  
  //========= reference values =========
  const double IRem = 1./ml_prob.GetInputParser().get("Rem");
//============================
  
//======== GEOMETRICAL ELEMENT =======
  const uint space_dim =       ml_prob._ml_msh->GetDimension();

//==== Operators @ gauss ======== 
  std::vector<double> Vel_vec_val_g(space_dim);
  std::vector<double> Vel_vec_val_g3D(3);
  std::vector<double> Bext_vec_val_g(space_dim);
  std::vector<double> Bext_vec_val_g3D(3);
  double         vXBe_g3D[3];
  double   vXBeXdphii_g3D[3];  
  double curlBeXdphii_g3D[3];
  
  double LapBe_g[DIMENSION];
//============================
  
        my_system._LinSolver[Level]->_KK->zero();
        my_system._LinSolver[Level]->_RESC->zero();

// ==========================================  
  Mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(Level);
  elem		*myel		=  mymsh->el;
  const unsigned myproc  = mymsh->processor_id();
	
// ==========================================  
// ==========================================
  
  {//BEGIN VOLUME

   const uint mesh_vb = VV;
   
    const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level+1];
    const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level];
 
    for (uint iel=0; iel < (nel_e - nel_b); iel++) {
   
    CurrentElem<double>       currelem(iel,myproc,Level,VV,&my_system,ml_prob.GetMeshTwo(),ml_prob.GetElemType(),mymsh);
    //   CurrentGaussPointBase & currgp = //   CurrentGaussPointBase::build(currelem,ml_prob.GetQuadratureRule(currelem.GetDim()));
    
//=========INTERNAL QUANTITIES (unknowns of the equation) ==================
    CurrentQuantity BhomOldX(currelem);
    BhomOldX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHom0"); 
    BhomOldX.VectWithQtyFillBasic();
    BhomOldX.Allocate();
    
    CurrentQuantity BhomOldY(currelem);
    BhomOldY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHom1"); 
    BhomOldY.VectWithQtyFillBasic();
    BhomOldY.Allocate();
    
    CurrentQuantity BhomOldZ(currelem);
    BhomOldZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHom2"); 
    BhomOldZ.VectWithQtyFillBasic();
    BhomOldZ.Allocate();
    
    std::vector<CurrentQuantity*> BhomOld_vec;   
    BhomOld_vec.push_back(&BhomOldX);
    BhomOld_vec.push_back(&BhomOldY);
    BhomOld_vec.push_back(&BhomOldZ);

    CurrentQuantity LagMultOld(currelem);
    LagMultOld._qtyptr   = ml_prob.GetQtyMap().GetQuantity("Qty_MagnFieldHomLagMult");
    LagMultOld.VectWithQtyFillBasic();
    LagMultOld.Allocate();
//========= END INTERNAL QUANTITIES (unknowns of the equation) =================

    
//========= tEST AND SHAPE FOR QTYZERO AND QTYONE =================
//QTYZERO SHAPE: shape of the first Unknown
    CurrentQuantity Phij(currelem); //TODO this is another Vect that doesnt have an associated quantity still
    Phij._dim      = 1;                                                         //scalar!
    Phij._FEord    = BhomOldX._FEord;
    Phij._ndof     = currelem.GetElemType(Phij._FEord)->GetNDofs(); 
    Phij.Allocate();
        
//QTYZERO tEST:  test of the first Unknown    
    CurrentQuantity Phii(currelem);
    Phii._dim      = 1;
    Phii._FEord    = BhomOldX._FEord;
    Phii._ndof     = currelem.GetElemType(Phii._FEord)->GetNDofs();
    Phii.Allocate();
    
//QTYONE SHAPE: Shape of the second Unknown
    CurrentQuantity Psij(currelem);
    Psij._dim      = 1;
    Psij._FEord    = LagMultOld._FEord;
    Psij._ndof     = currelem.GetElemType(Psij._FEord)->GetNDofs(); 
    Psij.Allocate();
    
//QTYONE tEST: test of the second Unknown
    CurrentQuantity Psii(currelem);
    Psii._dim      = 1;
    Psii._FEord    = LagMultOld._FEord;
    Psii._ndof     = currelem.GetElemType(Psii._FEord)->GetNDofs();
    Psii.Allocate();

//========= END tEST AND SHAPE FOR QTYZERO AND QTYONE =================

    
//=========EXTERNAL QUANTITIES (couplings) =====
  //========= //DOMAIN MAPPING
    CurrentQuantity xyz(currelem);
    xyz._dim      = DIMENSION;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();
    

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
  
    CurrentQuantity Bext_vecQuant(currelem);   //without quantity, nor equation
    Bext_vecQuant._dim     = Bext_vec.size();
    Bext_vecQuant._FEord   = Bext_vec[0]->_FEord;
    Bext_vecQuant._ndof    = Bext_vec[0]->_ndof;
    Bext_vecQuant.Allocate();
   
    CurrentQuantity VelX(currelem);
    VelX._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity0"); 
    VelX.VectWithQtyFillBasic();
    VelX.Allocate();
    
    CurrentQuantity VelY(currelem);
    VelY._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity1"); 
    VelY.VectWithQtyFillBasic();
    VelY.Allocate();
    
    CurrentQuantity VelZ(currelem);
    VelZ._qtyptr      = ml_prob.GetQtyMap().GetQuantity("Qty_Velocity2"); 
    VelZ.VectWithQtyFillBasic();
    VelZ.Allocate();
    
    std::vector<CurrentQuantity*> Vel_vec;   
    Vel_vec.push_back(&VelX);
    Vel_vec.push_back(&VelY);
    Vel_vec.push_back(&VelZ);

//=========END EXTERNAL QUANTITIES (couplings) =====
//======================
//======================
    
    
    currelem.Mat().zero();
    currelem.Rhs().zero();

    currelem.SetDofobjConnCoords();
     currelem.SetElDofsBc();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);

    
    for (uint idim=0; idim < space_dim; idim++)    {    
        BhomOld_vec[idim]->GetElemDofs();
            Vel_vec[idim]->GetElemDofs();
           Bext_vec[idim]->GetElemDofs();
    }
    
   Bext_vecQuant.GetElemDofs(Bext_vec); 
   
    
   const uint el_ngauss = ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussPointsNumber();

    for (uint qp = 0; qp < el_ngauss; qp++) {

//======= here starts the "COMMON SHAPE PART"==================
for (uint fe = 0; fe < QL; fe++)     { 
//   currgp.SetPhiElDofsFEVB_g (fe,qp);
//   currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
}

 const double      det = 1.;  //currgp.JacVectVV_g(xyz);
 const double dtxJxW_g = det*ml_prob.GetQuadratureRule(currelem.GetDim()).GetGaussWeight(qp);
 const double     detb = det/el_ngauss;

for (uint fe = 0; fe < QL; fe++)     {    
//   currgp.SetDPhiDxyzElDofsFEVB_g (fe,qp);
//   currgp.ExtendDphiDxyzElDofsFEVB_g (fe); 
}
//======= end of the "COMMON SHAPE PART"==================


    for (uint idim=0; idim < space_dim; idim++)  {
        BhomOld_vec[idim]->val_g();
            Vel_vec[idim]->val_g();            //---- for Advection MAT & RHS
      Vel_vec_val_g[idim] = Vel_vec[idim]->_val_g[0];
           Bext_vec[idim]->val_g();          //----- for Advection RHS
     Bext_vec_val_g[idim] = Bext_vec[idim]->_val_g[0];
           Bext_vec[idim]->grad_g();          //----- for Laplacian RHS
    }
    
      Bext_vecQuant.curl_g();          //----- for Curl Curl RHS         //THE EXTENSION of the DOFS to 3D is done INSIDE!!

       Math::extend(&Vel_vec_val_g[0],&Vel_vec_val_g3D[0],space_dim);                    //----- for Advection RHS
       Math::extend(&Bext_vec_val_g[0],&Bext_vec_val_g3D[0],space_dim);                    //----- for Advection RHS
       Math::cross(&Vel_vec_val_g3D[0],&Bext_vec_val_g3D[0],vXBe_g3D);          //----- for Advection RHS

//================================
//========= FILLING ELEMENT MAT/RHS
//=================================
// // //        for (uint i=0; i < BhomOldX._ndof; i++)     {
// // // //============ preparation for (i) (Phii) ============
// // //            Phii._val_g[0]       =      currgp._phi_ndsQLVB_g[Phii._FEord][i];          /*const double  phii_g*/ 
// // //         for (uint idim=0; idim<space_dim; idim++)
// // // 	  Phii._grad_g[0][idim] = currgp._dphidxyz_ndsQLVB_g[Phii._FEord][i+idim*Phii._ndof];   /*dphiidx_g*/
// // // //=== from this point on, the dependency on /*(i)*/, the ROW INDEX, of the quantities is given by Phii  
// // // //===  /*(i)*/ is ONLY used for the POSITIONS, later
// // // 	  
// // //     //--------- LAPLACIAN: Operator, RHS: grad Be . grad phi ---------------------
// // // 	 for (uint idim=0; idim<space_dim/*Bext._dim*/; idim++) LapBe_g[idim]=0.;
// // // 
// // //           for (uint idim=0; idim<space_dim/*Bext._dim*/; idim++) {
// // //  	    for (uint jdim=0; jdim<space_dim; jdim++) {
// // //               LapBe_g[idim] += Bext_vec[idim]->_grad_g[0][jdim]*Phii._grad_g[0][jdim];
// // // 	    }
// // //  	  }
// // // 
// // //              Math::extend(&Phii._grad_g[0][0],&Phii._grad_g3D[0][0],space_dim);
// // // 
// // //     //--------- CURL CURL: Operator, RHS: curl Be . curl phi -------------------
// // //              Math::cross(&Bext_vecQuant._curl_g3D[0],&Phii._grad_g3D[0][0],curlBeXdphii_g3D);
// // // 
// // //     //--------- ADVECTION: Operator, RHS: v x Be . curl phi -------------------
// // //              Math::cross(      vXBe_g3D,&Phii._grad_g3D[0][0],  vXBeXdphii_g3D);       // _utils.cross(vXBe_g3D,dphiidx_g3D,vXBeXdphii_g3D);
// // // //============end preparation for (i) ============
// // // 	     
// // //          for (uint idim=0; idim<space_dim/*bhomOld._dim*/; idim++) {
// // //             const uint irowq = i + idim*BhomOldX._ndof;
// // //             
// // //             currelem.Rhs()(irowq) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                            - LAP_MHD*IRem*LapBe_g[idim]
// // //                             - (1-LAP_MHD)*IRem*curlBeXdphii_g3D[idim]              //phii of bhomOld /*CurlCurl(RHS,vb,Phij,Phii,idim,idimp1);*/
// // //                            + ADV_MHD*           vXBeXdphii_g3D[idim]               //phii of bhomOld
// // //                          )
// // //                          + (1-currelem.GetBCDofFlag()[irowq])*detb*BhomOld_vec[idim]->_val_dofs[i]; //Dirichlet bc
// // // 	   }
// // // 
// // // 
// // //         for (uint idim=0; idim<space_dim/*bhomOld._dim*/; idim++) { // filling diagonal block for Dirichlet bc
// // //           const uint irowq = i + idim*BhomOldX._ndof;
// // //           currelem.Mat()(irowq,irowq) += (1-currelem.GetBCDofFlag()[irowq])*detb;
// // //         }
// // //                                            // end filling diagonal for Dirichlet bc
// // // 	 
// // //         for (uint j=0; j<BhomOldX._ndof; j++) {  // A element matrix
// // // //============ preparation for (j) (Phij) ============
// // //             Phij._val_g[0]        =      currgp._phi_ndsQLVB_g[Phij._FEord][j]; 
// // //           for (uint idim=0; idim<space_dim; idim++)
// // // 	    Phij._grad_g[0][idim] = currgp._dphidxyz_ndsQLVB_g[Phij._FEord][j+idim*Phij._ndof];   //real shape  /*dphijdx_g[idim]*/
// // // //=== from this point on, the dependency on /*(j)*/, the COLUMN INDEX, of the quantities is given by Phij  
// // // //===  /*(j)*/ is ONLY used for the POSITIONS, later
// // // 
// // //     //--------- LAPLACIAN: Operator, MAT: grad b . grad phi ---------------------
// // //           double     Lap_g = Math::dot(&Phij._grad_g[0][0],&Phii._grad_g[0][0],space_dim);   /*(i,j)*/  //part independent of idim
// // //     //--------- ADVECTION: Operator, MAT: v x b . curl phi -------------------
// // //           double Advphii_g = Math::dot( &Vel_vec_val_g[0],&Phii._grad_g[0][0],space_dim);  /*(i)*/ //part independent of idim //TODO what about putting it OUTSIDE?
// // //      
// // //           for (uint idim=0; idim<space_dim/*bhomOld._dim*/; idim++) { //filled in as 1-2-3 // 5-6-4 // 9-7-8
// // //             int irowq=i+idim*BhomOldX._ndof;   //idim gives the row index  //test of bhomOld
// // //             // diagonal blocks [1-5-9]  idim = row index = column index  //shape of bhomOld
// // //             currelem.Mat()(irowq,j+idim*BhomOldX._ndof)
// // //             += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                  + LAP_MHD*IRem*(Lap_g) 
// // //                  + (1-LAP_MHD)*IRem*(           Lap_g -  Phij._grad_g[0][idim]*Phii._grad_g[0][idim] )
// // //                  - ADV_MHD * Phij._val_g[0]* (Advphii_g -   Vel_vec[idim]->_val_g[0]    *Phii._grad_g[0][idim] )   //TODO Phij here does not depend on idim, but it depends on j
// // //                );
// // //             // block +1 [2-6-7]
// // //             int idimp1=(idim+1)%space_dim;
// // // 	    //idim: component of the SHAPE
// // // 	    //idimp1: component of the tEST
// // //             currelem.Mat()(irowq,j+idimp1*BhomOldX._ndof)
// // //             += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                   + (1-LAP_MHD)*IRem*(        - Phij._grad_g[0][idim]*Phii._grad_g[0][idimp1] )       /*(i,j)*/    /*CurlCurl(MAT,vb,Phij,Phii,idim,idimp1);*/
// // //                   - ADV_MHD * Phij._val_g[0]* ( -  Vel_vec[idim]->_val_g[0] * Phii._grad_g[0][idimp1] )       /*(i,j)*/    /*AdvCurl(MAT,Vel,Phij,Phii,idim,idimp1)*/
// // //                );
// // // 
// // // #if (DIMENSION==3)
// // //             // block +2 [3-4-8] 
// // // 	    int idimp2=(idim+2)%space_dim;//idimp2 column index
// // //             currelem.Mat()(irowq,j+idimp2*BhomOldX._ndof)
// // //             += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*(
// // //                   + (1-LAP_MHD)*IRem*(       - Phij._grad_g[0][idim]*Phii._grad_g[0][idimp2] )               /*(i,j)*/
// // //                   - ADV_MHD * Phij._val_g[0]* ( -  Vel_vec[idim]->_val_g[0]    *Phii._grad_g[0][idimp2] )               /*(i,j)*/
// // //                );
// // // #endif
// // //           }
// // //                     
// // //         } 
// // //                                       // end A element matrix
// // // 				      
// // //        for (uint j = 0; j < LagMultOld._ndof; j++) {// B^T element matrix ( p*div(v) )
// // //           Psij._val_g[0] = currgp._phi_ndsQLVB_g[Psij._FEord][j];   /*const double psij_g */
// // //           const int jclml= j+/*bhomOld._dim*/space_dim*BhomOldX._ndof;
// // //           for (uint idim=0; idim<space_dim; idim++) {  //bhomOld._dim==spacedimension
// // //             uint irowq=i+idim*BhomOldX._ndof;
// // //             currelem.Mat()(irowq,jclml) += currelem.GetBCDofFlag()[irowq]*dtxJxW_g*( - Psij._val_g[0]*Phii._grad_g[0][idim] );   /*(i,j)*/ /*PDiv(MAT,Psij,Phii,idim)*/
// // //            }
// // //         }
// // //                                      // end B^T element matrix
// // //  //============ begin linear rows
// // //           if (i < LagMultOld._ndof) {
// // //            Psii._val_g[0]  = currgp._phi_ndsQLVB_g[Psii._FEord][i];   /*double psii_g*/
// // // 	  const uint irowl=i+space_dim*BhomOldX._ndof;
// // //           currelem.Rhs()(irowl)=0.;
// // // 
// // //           for (uint j=0; j<BhomOldX._ndof; j++) { // B element matrix q*div(u)
// // //             for (uint idim=0; idim<space_dim; idim++) Phij._grad_g[0][idim] =  currgp._dphidxyz_ndsQLVB_g[Phij._FEord][j+idim*Phij._ndof];
// // //             for (uint idim=0; idim<space_dim; idim++) currelem.Mat()(irowl,j+idim*BhomOldX._ndof) += - dtxJxW_g*Psii._val_g[0]*Phij._grad_g[0][idim];  /*(i,j)*/  /*PDiv(MAT,Psii,Phij,idim)*/
// // //                 }
// // // 
// // //         }
// // // //============ end linear rows
// // //       }
// // //     //end filling element mat/rhs (i loop)
    }
    // end element gaussian integration loop
    
    ///  Add element matrix and rhs to the global ones.
                      my_system._LinSolver[Level]->_KK->add_matrix(currelem.Mat(),currelem.GetDofIndices());
                      my_system._LinSolver[Level]->_RESC->add_vector(currelem.Rhs(),currelem.GetDofIndices()); 

  } 
  // end of element loop
  
 
  }//END VOLUME
  

        my_system._LinSolver[Level]->_KK->close();
        my_system._LinSolver[Level]->_RESC->close();
    
#ifdef DEFAULT_PRINT_INFO
 std::cout << " GenMatRhs " << my_system.name() << ": assembled  Level " << Level
           << " with " << my_system._LinSolver[Level]->_KK->m() << " dofs " << std::endl;
#endif     


  return;
}


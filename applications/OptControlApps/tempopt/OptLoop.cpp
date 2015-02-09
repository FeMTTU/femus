#include "OptLoop.hpp"

#include "MultiLevelProblem.hpp"
#include "CurrentGaussPoint.hpp"
#include "CurrentElem.hpp"
#include "CurrentQuantity.hpp"
#include "Quantity.hpp"
#include "ElemType.hpp"
#include "Box.hpp"
#include "paral.hpp"
#include "Files.hpp"
#include "XDMFWriter.hpp"


namespace femus {
  
  
 OptLoop::OptLoop(Files& files_in, const FemusInputParser<double> & map_in): TimeLoop(files_in, map_in) { }

 //=================
    void OptLoop::optimization_loop(MultiLevelProblem & eqmap_in)  {

    //  parameters
    int    print_step = _timemap.get("printstep");

    double curr_time = _time_in;  //initialize current time
    
 for (uint curr_step = _t_idx_in + 1; curr_step <= _t_idx_final; curr_step++) {

   curr_time += 1.;

#if DEFAULT_PRINT_TIME==1
      std::clock_t  start_time=std::clock();
#endif

       _curr_t_idx = curr_step;
       _curr_time  = curr_time;

      std::cout << "\n  ** Solving time step " << _curr_t_idx
                << ", time = "                 << _curr_time   << " ***" << std::endl;

       
	const uint delta_t_step = curr_step - _t_idx_in;

      //  time step for each system, without printing (good)
      OneTimestepEqnLoop(delta_t_step,eqmap_in);

#if DEFAULT_PRINT_TIME==1
      std::clock_t    end_time=std::clock();
#endif 

      // print solution
      if (delta_t_step%print_step == 0) XDMFWriter::PrintSolLinear(_files.GetOutputPath(),curr_step,curr_time,eqmap_in);   //print sol.N.h5 and sol.N.xmf
    

#if DEFAULT_PRINT_TIME==1
      std::clock_t    end_time2=std::clock();
      std::cout << " Time solver ----->= "   << double(end_time- start_time)/ CLOCKS_PER_SEC
                << " Time printing ----->= " << double(end_time2- end_time) / CLOCKS_PER_SEC <<
                std::endl;
#endif 

//=====functional evaluations=======
   SystemTwo & eqnT = eqmap_in.get_system<SystemTwo>("Eqn_T");

		
     double J = 0.;
J = ComputeIntegral    ( eqmap_in.GetMeshTwo()._NoLevels - 1,&eqmap_in.GetMeshTwo(),&eqnT,_files.GetOutputTime());
J = ComputeNormControl ( eqmap_in.GetMeshTwo()._NoLevels - 1,&eqmap_in.GetMeshTwo(),&eqnT,0 );
J = ComputeNormControl ( eqmap_in.GetMeshTwo()._NoLevels - 1,&eqmap_in.GetMeshTwo(),&eqnT,1 );
//=====functional evaluations =======


    }   // end time loop

    return;
  }
  

  
// This function computes the integral only for the current processor

double ComputeIntegral (const uint Level, const MultiLevelMeshTwo* mesh, const SystemTwo* eqn, const std::string output_time) {

  //====== processor index
  const uint myproc = mesh->_iproc;
  const uint space_dim =      mesh->get_dim();
  const uint mesh_ord = (int) mesh->GetRuntimeMap().get("mesh_ord");  
  const uint meshql   = (int) mesh->GetRuntimeMap().get("meshql");   //======== ELEMENT MAPPING =======
 
  const uint mesh_vb = VV;

    CurrentElem       currelem(Level,VV,eqn,*mesh,eqn->GetMLProb().GetElemType());
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,eqn->GetMLProb().GetQrule(currelem.GetDim()));

  //========== 
    CurrentQuantity Tempold(currgp);
    Tempold._qtyptr   =  eqn->GetMLProb().GetQtyMap().GetQuantity("Qty_Temperature"); 
    Tempold.VectWithQtyFillBasic();
    Tempold.Allocate();

  //========== 
    CurrentQuantity Tlift(currgp);
    Tlift._qtyptr   =  eqn->GetMLProb().GetQtyMap().GetQuantity("Qty_TempLift"); 
    Tlift.VectWithQtyFillBasic();
    Tlift.Allocate();
    
 //===========
    CurrentQuantity Tdes(currgp);
    Tdes._qtyptr   = eqn->GetMLProb().GetQtyMap().GetQuantity("Qty_TempDes"); 
    Tdes.VectWithQtyFillBasic();
    Tdes.Allocate();
  
//========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

//========== Quadratic domain, auxiliary  
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = NVE[ mesh->_geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];
  xyz_refbox.Allocate();

  
   double integral = 0.;

      const uint el_ngauss = eqn->GetMLProb().GetQrule(currelem.GetDim()).GetGaussPointsNumber();

//parallel sum
    const uint nel_e = mesh->_off_el[mesh_vb][mesh->_NoLevels*myproc+Level+1];
    const uint nel_b = mesh->_off_el[mesh_vb][mesh->_NoLevels*myproc+Level];
  
    for (uint iel=0; iel < (nel_e - nel_b); iel++) {

      currelem.SetDofobjConnCoords(myproc,iel);
      currelem.SetMidpoint();
      
     currelem.ConvertElemCoordsToMappingOrd(xyz);
     currelem.TransformElemNodesToRef(eqn->GetMLProb().GetMeshTwo().GetDomain(),&xyz_refbox._val_dofs[0]);    

// =============== 
      xyz_refbox.SetElemAverage();
      int el_flagdom = ElFlagControl(xyz_refbox._el_average,mesh);
//====================     
 
    if ( Tempold._eqnptr != NULL )   Tempold.GetElemDofs();
    else                             Tempold._qtyptr->FunctionDof(Tempold,0.,&xyz_refbox._val_dofs[0]);
    if ( Tlift._eqnptr != NULL )       Tlift.GetElemDofs();
    else                               Tlift._qtyptr->FunctionDof(Tlift,0.,&xyz_refbox._val_dofs[0]);
    if ( Tdes._eqnptr != NULL )         Tdes.GetElemDofs();
    else                                Tdes._qtyptr->FunctionDof(Tdes,0.,&xyz_refbox._val_dofs[0]);    



    for (uint qp = 0; qp < el_ngauss; qp++) {

     for (uint fe = 0; fe < QL; fe++)     {  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  }  
     
   const double  Jac_g = currgp.JacVectVV_g(xyz);  //not xyz_refbox!      
   const double  wgt_g = eqn->GetMLProb().GetQrule(currelem.GetDim()).GetGaussWeight(qp);

     for (uint fe = 0; fe < QL; fe++)     {          currgp.SetPhiElDofsFEVB_g (fe,qp);  }

 Tempold.val_g();
   Tlift.val_g();
    Tdes.val_g();

  double deltau_squarenorm_g = 0.;
   deltau_squarenorm_g += (Tempold._val_g[0] + Tlift._val_g[0] - Tdes._val_g[0])*
                          (Tempold._val_g[0] + Tlift._val_g[0] - Tdes._val_g[0]); 

  integral += el_flagdom*wgt_g*Jac_g*deltau_squarenorm_g;
   
    }//gauss loop
     
    }//element loop
    
////////////////////////////////////////        
       std::cout << "integral on processor 0: " << integral << std::endl;

   double J=0.;
#ifdef HAVE_MPI
      MPI_Allreduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  //THIS IS THE RIGHT ONE!!
#else
      J = integral;
#endif
   
    std::cout << "@@@@@@@@@@@@@@@@ functional value: " << J << std::endl;
//=====functional evaluation =======   

    
    ///////// let us also print the functional value in a unique file,
    /////////so that we explore the variation wrt alpha
    
    std::string app_path = "./";
    std::string intgr_fname = app_path + DEFAULT_OUTPUTDIR + "/" + "alpha";
  
	std::ofstream intgr_fstream;
    
    if (paral::get_rank() ==0 ) { 
      intgr_fstream.open(intgr_fname.c_str(),std::ios_base::app); 
      intgr_fstream << output_time << " " << eqn->GetMLProb().GetInputParser().get("alphaT") << " " << eqn->GetMLProb().GetInputParser().get("injsuc")<< " "  << J << " " << std::endl ; 
      intgr_fstream.close();  //you have to close to disassociate the file from the stream
    }
 
   
  return J;  
    
}


/////////////////////////

double ComputeNormControl (const uint Level, const MultiLevelMeshTwo* mesh, const SystemTwo* eqn, const uint reg_ord ) {

  //reg_ord = 0: L2
  //reg_ord = 1: H1

  // processor index
  const uint myproc = mesh->_iproc;
  const uint space_dim =       mesh->get_dim();
  const uint mesh_ord  = (int) mesh->GetRuntimeMap().get("mesh_ord");  
  const uint meshql    = (int) mesh->GetRuntimeMap().get("meshql");    //======== ELEMENT MAPPING =======

  const uint mesh_vb = VV;
  
    CurrentElem       currelem(Level,VV,eqn,*mesh,eqn->GetMLProb().GetElemType());
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,eqn->GetMLProb().GetQrule(currelem.GetDim()));
  
//======Functions in the integrand ============

      //========== 
    CurrentQuantity Tlift(currgp);
    Tlift._qtyptr   =  eqn->GetMLProb().GetQtyMap().GetQuantity("Qty_TempLift"); 
    Tlift.VectWithQtyFillBasic();
    Tlift.Allocate();

//========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = space_dim;
    xyz._FEord    = meshql;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

//========== Quadratic domain, auxiliary  
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = space_dim;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = NVE[ mesh->_geomelem_flag[currelem.GetDim()-1] ][BIQUADR_FE];
  xyz_refbox.Allocate();
  
    
   double integral = 0.;

//loop over the geom el types
      const uint el_ngauss = eqn->GetMLProb().GetQrule(currelem.GetDim()).GetGaussPointsNumber();

//parallel sum
    const uint nel_e = mesh->_off_el[mesh_vb][mesh->_NoLevels*myproc+Level+1];
    const uint nel_b = mesh->_off_el[mesh_vb][mesh->_NoLevels*myproc+Level];
  
    for (int iel=0; iel < (nel_e - nel_b); iel++) {

      currelem.SetDofobjConnCoords(myproc,iel);
      currelem.SetMidpoint();

      currelem.ConvertElemCoordsToMappingOrd(xyz);
      currelem.TransformElemNodesToRef(eqn->GetMLProb().GetMeshTwo().GetDomain(),&xyz_refbox._val_dofs[0]);    
     
     Tlift.GetElemDofs();


  for (uint qp = 0; qp < el_ngauss; qp++) {

     for (uint fe = 0; fe < QL; fe++)     {
       currgp.SetPhiElDofsFEVB_g (fe,qp);
       currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  
    }  
     
      const double  Jac_g = currgp.JacVectVV_g(xyz);  //not xyz_refbox!      
      const double  wgt_g = eqn->GetMLProb().GetQrule(currelem.GetDim()).GetGaussWeight(qp);

  Tlift.val_g();
  Tlift.grad_g();

  double deltau_squarenorm_g = 0.;
   deltau_squarenorm_g += (1-reg_ord) * (Tlift._val_g[0])*(Tlift._val_g[0]); 

   for (uint idim = 0; idim < space_dim; idim++)  {   deltau_squarenorm_g  += reg_ord * (Tlift._grad_g[0][idim])*(Tlift._grad_g[0][idim])  ;   }

  integral += /*el_flagdom**/wgt_g*Jac_g*deltau_squarenorm_g;   //Do it over ALL THE DOMAIN!
   
    }//gauss loop
     
    }//element loop
    
       std::cout << "integral on processor 0: " << integral << std::endl;

   double J = 0.;
#ifdef HAVE_MPI
      MPI_Allreduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  //THIS IS THE RIGHT ONE!!
#else
      J = integral;
#endif
   
    std::cout << "@@@@@@@@@@@@@@@@ control norm: " << reg_ord << " " << J << std::endl;
    
    
  return J;  
  
}
  
  // ========================================================
 //the Box must already be initialized here
 //it must receive an el_xm already in the REFBox frame

///AAA questa funzione NON lavora se tu fai solo DUE SUDDIVISIONI PER LATO e nolevels=1 !!!

 int ElFlagControl(const std::vector<double> el_xm, const MultiLevelMeshTwo* mesh)  {

  Box* box= static_cast<Box*>(mesh->GetDomain());
   
     int el_flagdom=0;

     if (mesh->get_dim() == 2) {
//============== OUTFLOW     
//        if (   el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])
// 	   && el_xm[0] < 0.75*(box->_le[0] - box->_lb[0])
// 	   && el_xm[1] > 0.75*(box->_le[1] - box->_lb[1]) ) {
//                  el_flagdom=1;
//              }
//============== RIGHT
//      if (   el_xm[0] > 0.85*(box->_le[0] - box->_lb[0])
// 	   && el_xm[1] > 0.25*(box->_le[1] - box->_lb[1])
// 	   && el_xm[1] < 0.75*(box->_le[1] - box->_lb[1]) ) {
//                  el_flagdom=1;
//              }
//============== INFLOW
// facciamo 0.2-0.4; 0.4-0.6; 0.6-0.8     
     if (     el_xm[1] > 0.4*(box->_le[1] - box->_lb[1])
           && el_xm[1] < 0.6*(box->_le[1] - box->_lb[1])
	   && el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])
	   && el_xm[0] < 0.75*(box->_le[0] - box->_lb[0]) ) {
                 el_flagdom=1;
             }
// //============== CENTER
//      if (  /*   el_xm[1] > 0.04*(box->_le[1] - box->_lb[1])*/
//           /* &&*/ el_xm[1] < 0.06*(box->_le[1] - box->_lb[1])
// 	   && el_xm[0] > 0.4*(box->_le[0] - box->_lb[0])
// 	   && el_xm[0] < 0.6*(box->_le[0] - box->_lb[0]) ) {
//                  el_flagdom=1;
//              }
     }
 
 else if (mesh->get_dim() == 3) {
 
      if ( el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])  
	&& el_xm[0] < 0.75*(box->_le[0] - box->_lb[0]) 
	&& el_xm[1] > 0.75*(box->_le[1] - box->_lb[1])
	&& el_xm[2] > 0.25*(box->_le[2] - box->_lb[2]) 
	&& el_xm[2] < 0.75*(box->_le[2] - box->_lb[2]) ) {
	el_flagdom=1;
        }
        
 }

return el_flagdom; 

}
 


} //end namespace femus


  
  
#include "OptLoop.hpp"


#include <cmath>

//library headers
#include "NumericVector.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "SystemTwo.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "TimeLoop.hpp"
#include "Files.hpp"
#include "CurrentElem.hpp"
#include "CurrentGaussPointBase.hpp"
#include "CurrentQuantity.hpp"
#include "Quantity.hpp"
#include "ElemType.hpp"
#include "Box.hpp"

#include "paral.hpp"


//application headers
#include "Opt_conf.hpp"
#include "EqnNS.hpp"
#include "EqnMHDCONT.hpp"
#include "EqnMHD.hpp"
#include "EqnNSAD.hpp"
#include "EqnMHDAD.hpp"

namespace femus {

  
  
 //INITIALIZE FSTREAM FOR INTEGRAL
   //  create an ofstream integral.txt
   //the question is: who creates it? if one processor creates it,
//     then can the others create it again?
// or should we create it with processor 0?
//no because the others wouldnt have it.
// every process must have the fstream, 
//but only one process must CREATE the FILE associated to this filestream
//also,every process must associate this fstream with that file.
//I guess that what happens is: the process wants to associate the fstream
//with the file: If the file exists,then it creates it, otherwise it uses it.
//and that should be ok.
//before, the creation of the stream was INSIDE the time loop
//So the program did that: when you create a stream by associating it to a file,
//if the file already exists then you clean the file ad restart from scratch.
//the best thing is:CREATE the stream ONLY ONCE.
//then at every loop step, OPEN and CLOSE the stream, to be perfectly safe.

//what seems to be strange is that every process should have this stream but
//there is only one print in the file... it seems like there is some
//sort of automatic mpi in the streams... but there is no such automation 
// for std::cout, for instance



///This optimization loop is just another way of picking the various equations in a certain order rather than another.
///The equations are stationary, but we have to mimic the algorithmic sequence as a time sequence,
///so that we can see the algorithmic steps in Paraview

//BEFORE, the time loop was actually a COUPLED NONLINEAR LOOP where 
// every equation was advanced by one nonlinear step altogether.
//That brought to a lot of oscillations, that were transferred to 
//all the equations,both the nonlinear ones and the linear ones.
//NOW, the time step/value printed in file will just be an OPTIMIZATION STEP
// so the outer loop will be OPTIMIZATION
//then, inside, we will have NONLINEAR LOOP for the DIRECT equations

//   _t_idx_in and _t_idx_final can be considered as OPTIMIZATION STEPS rather than TIME STEPS
//there is no TIME VALUE involved in the optimization steps
//anyway, it is better to update it as well, so that when you read time in paraview
//you have a different value every time
//Paraview doesnt like time.xmf with files having the same time value, so change it as well.
//once the TIME was a nonlinear loop
//now it is an optimization loop
//question:how do I see the NONLINEAR LOOPS alone?
//you call transient_loop() instead of optimization_loop()
//x_old and _x_oold will be used 
//either for OPTIMIZATION iterations or for NONLINEAR iterations,
//depending on what you need

 OptLoop::OptLoop(Files& files_in): TimeLoop(files_in) { }

 OptLoop::~OptLoop()  {
  
     for (uint Level = 0; Level < _x_oldopt.size(); Level++) {
          delete _x_oldopt[Level];
      }
    
      _x_oldopt.clear();

  }
 
 
 
 
void OptLoop::optimization_loop(MultiLevelProblemTwo& e_map_in)  {
  
  
 #ifdef NS_EQUATIONS
     EqnNS* mgNS = static_cast<EqnNS*>(e_map_in.get_eqs("Eqn_NS"));
  #endif
  #ifdef MHD_EQUATIONS
     EqnMHD* mgMHD = static_cast<EqnMHD*>(e_map_in.get_eqs("Eqn_MHD"));
  #endif
 #ifdef NSAD_EQUATIONS
     EqnNSAD* mgNSAD = static_cast<EqnNSAD*>(e_map_in.get_eqs("Eqn_NSAD"));
 #endif
  #ifdef MHDAD_EQUATIONS
      EqnMHDAD* mgMHDAD = static_cast<EqnMHDAD*>(e_map_in.get_eqs("Eqn_MHDAD"));
  #endif
  #ifdef MHDCONT_EQUATIONS
      EqnMHDCONT* mgMHDCONT = static_cast<EqnMHDCONT*>(e_map_in.get_eqs("Eqn_MHDCONT"));
  #endif


std::string intgr_fname = e_map_in._files._output_path + "/" + "integral.txt";

 std::ofstream intgr_fstream;

//INITIALIZE OPT LOOP
//pseudo time parameters for optimization
    const uint      NoLevels = e_map_in._mesh._NoLevels;
    double pseudo_opttimeval = _time_in;
    int           print_step = _timemap.get("printstep");

     double omega = 1.;
     double Jold = 10.; //TODO AAA
     double J = 0.;

//parameter for nonlinear loops and coupling loops
    const uint N_knonl_NS = 3.;
//     const uint N_knonl_MHD = 3.;
const uint N_faketime = 10.;

double eps_NS_MHD = 1.e-3;
double eps_nl_NS = 1.e-4;
double eps_MHD = 1.e-4;

double eps_ADJCONT = 1.e-3;
double eps_NSAD = 1.e-4;
double eps_MHDAD = 1.e-4;
double eps_MHDCONT = 1.e-4;

double epsJ = 1.e-12;

const uint MaxIterNS = 10; 
const uint MaxIterMHD = 10; 
const uint MaxIterNSMHD = 10; 
const uint MaxIterADJCONT = 10; 
const uint MaxIterNSAD = 20; 
const uint MaxIterMHDAD = 20;
const uint MaxIterMHDCONT = 20;

double nonlin_deltax_NS = 0.;
double nonlin_deltax_MHD = 0.;
double lin_deltax_NSAD = 0.;
double lin_deltax_MHDAD = 0.;
double lin_deltax_MHDCONT = 0.;

//initialize specific data for specific equations
#ifdef MHDCONT_EQUATIONS
     init_equation_data(mgMHDCONT);
#endif
  
//nitialize Becont
  _x_oldopt[NoLevels - 1]->zero();    //initialize Boldopt=0;


  //OPTIMIZATION LOOP
for (uint opt_step = _t_idx_in + 1; opt_step <= _t_idx_final; opt_step++) {
  
    pseudo_opttimeval += 1.; //just to increase the time value

  std::cout << "\n  @@@@@@@@@@@@@@@@ Solving optimization step " << opt_step << std::endl;

#ifdef MHDCONT_EQUATIONS 
   std::cout << "\n  ** Compute the Becont to be given to the STATE equations " << opt_step << std::endl;
// 	B_old = B_oold + omega*(B_old - B_oold) = omega*B_old + (1-omega)*B_oold;

//There is a problem... when you arrive here after it said: reduce the intensity, then in
//xold there was a good value
//xoold there was still the value multiplied
//   omega = 1.; //"omega=1" and "no if" is equal to the old loop
    mgMHDCONT->_x_oold[NoLevels - 1]->close();
    std::cout << "Linfty norm of Becont _x_oold " << mgMHDCONT->_x_oold [NoLevels - 1]->linfty_norm() << std::endl;

    mgMHDCONT->_x_old[NoLevels - 1]->close();
    std::cout << "Linfty norm of Becont _x_old " << mgMHDCONT->_x_old [NoLevels - 1]->linfty_norm() << std::endl;

    _x_oldopt[NoLevels - 1]->close();
    std::cout << "Linfty norm of Becont _x_oldopt " << _x_oldopt [NoLevels - 1]->linfty_norm() << std::endl;

//////////////////

      mgMHDCONT->_x_tmp[NoLevels-1]->zero();
    *(mgMHDCONT->_x_tmp[NoLevels-1]) = *(mgMHDCONT->_x_old[NoLevels-1]);
      mgMHDCONT->Bc_ScaleDofVec(mgMHDCONT->_x_tmp[NoLevels - 1], omega );
      mgMHDCONT->_x_tmp[NoLevels - 1]->close();
    std::cout << "Omega " << omega << std::endl;
    std::cout << "Linfty norm of Becont _x_old*omega "
              << mgMHDCONT->_x_tmp [NoLevels - 1]->linfty_norm() << std::endl;

      mgMHDCONT->_x_oold[NoLevels - 1]->close();
    std::cout << "Linfty norm of Becont _x_oold "
              << mgMHDCONT->_x_oold [NoLevels - 1]->linfty_norm() << std::endl;

     _x_oldopt[NoLevels - 1]->close();
    std::cout << "Linfty norm of Becont _x_oldopt "
              << _x_oldopt [NoLevels - 1]->linfty_norm() << std::endl;

      mgMHDCONT->Bc_AddScaleDofVec(_x_oldopt[NoLevels - 1],mgMHDCONT->_x_tmp [NoLevels - 1],1.- omega);
      mgMHDCONT->_x_tmp[NoLevels - 1]->close();
    std::cout << "Linfty norm of Becont x_old*omega + (1-omega)*xoold " 
              << mgMHDCONT->_x_tmp [NoLevels - 1]->linfty_norm() << std::endl;

    *(mgMHDCONT->_x_old[NoLevels-1]) = *(mgMHDCONT->_x_tmp[NoLevels-1]);

    mgMHDCONT->_x_old[NoLevels - 1]->close();
    std::cout << "Linfty norm of Becont _x_old updated " 
              << mgMHDCONT->_x_old [NoLevels - 1]->linfty_norm() << std::endl;

 
//   mgMHDCONT->x[NoLevels - 1]->close();
//   std::cout << "Linfty norm of Becont x " << mgMHDCONT->x[NoLevels - 1]->linfty_norm() << std::endl;
//NOTICE that at this point x_old for MHDCONT has been filled by ic_read the first time
//for the other times it is updated by every solve

  //i compute B_old because we communicate between eqns with the old values    
//here pay attention, xold will be used by mhdcont again
 #endif


{
  
  uint coupl_NS_MHD = 0;

do {
 
    coupl_NS_MHD++;

    std::cout << "\n @@@@@@@@@@@@@@@@ Solving NS+ MHD in uncoupled algorithm, iteration " << coupl_NS_MHD << std::endl;

#ifdef NS_EQUATIONS 
{ 
  uint knonl_NS =0;
  do {
    
    knonl_NS++;
//      std::cout << "\n >>>>> Solving nonlinear step " << knonl_NS << " for" << mgNS->_eqname << std::endl;
    nonlin_deltax_NS = MGTimeStep(knonl_NS,mgNS);

  } while ( nonlin_deltax_NS > eps_nl_NS && knonl_NS < MaxIterNS );
}

#endif 

#ifdef MHD_EQUATIONS

{
   uint k_MHD =0;
do {
  
  k_MHD++;
//ONE nonlinear step = ONE LINEAR SOLVER  
//    std::cout << "\n >>>>>>>> Solving MHD system (linear in B), " << k_MHD << std::endl;
    nonlin_deltax_MHD = MGTimeStep(k_MHD,mgMHD);

}while (nonlin_deltax_MHD > eps_MHD &&  k_MHD < MaxIterMHD );
  
}
#endif

  std::cout << "\n @@@@@@@@@@@@@@@@ Coupling iteration NS+MHD " << coupl_NS_MHD  << std::endl;
  std::cout << " @@@@@@@@@@@@@@@@  Overall NS+MHD system error " << nonlin_deltax_NS + nonlin_deltax_MHD << std::endl;

 
}//CONVERGENCE CRITERION FOR THE fake coupling LOOP
// if || u - u_old || and || B - u_old || < eps, quit the loop
  while ( /*(nonlin_deltax_NS + nonlin_deltax_MHD) > eps_NS_MHD &&*/ coupl_NS_MHD < MaxIterNSMHD );
 //if you remove the condition on the sum of the norm then it may get the Hartmann already in the first step
 
} //it seems like you need these brackets for the do-while loop...


//************** compute J ********************

double integral = 0.;
#ifdef NS_EQUATIONS      
   integral = ComputeIntegral(e_map_in._mesh._NoLevels - 1,&e_map_in._mesh,mgNS);
 #endif
  
   std::cout << "integral on processor 0: " << integral << std::endl;

  double J=0.;
#ifdef HAVE_MPI
      MPI_Allreduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  //THIS IS THE RIGHT ONE!!
#else
   J = integral;
#endif
   
    std::cout << "@@@@@@@@@@@@@@@@ functional value: " << J << std::endl;

if (paral::get_rank() ==0 ){ 
      intgr_fstream.open(intgr_fname.c_str(),ios_base::app); 
      intgr_fstream << opt_step << " " << pseudo_opttimeval << " " << J << " " << omega /*<< std::endl*/; 
      intgr_fstream.close();  //you have to close to disassociate the file from the stream
}
    
//********* see if J<J_old ********************
 //do with/without scaling (with/without ||1 see what happens)
if ( fabs(J - Jold) > epsJ /*|| 1*/  ) {
  std::cout << "@@@@@@@@@@@@@@@@ Jold= " << Jold << " , J= " << J << std::endl;

  if (J < Jold  /*|| 1*/ ) {  

    if ( paral::get_rank() == 0 ){ 
            intgr_fstream.open(intgr_fname.c_str(),ios_base::app); 
            intgr_fstream << "@@@@@@@@@@@@@@@@ J < Jold" << std::endl;
            intgr_fstream.close();  //you have to close to disassociate the file from the stream
        }
//******* update Jold  //you must update it only here, because here it is the good point to restart from
    Jold = J;
#ifdef MHDCONT_EQUATIONS      
        *(_x_oldopt[NoLevels-1]) = *(mgMHDCONT->_x_old[NoLevels-1]); 
#endif      

	//this will be the new _x_oldopt?
    //if this update is such that J> Jold, then I should not refill _x_oldopt here, but only inside J<Jold
    //this update must be only done inside the J<Jold case
    //in the other cases _x_oldopt is not reupdated and remains the same, just like J_old!
        
    
//The problem is that the criterion J<Jold does not accept any possible oscillations,
//which may occur even though the system has been split into the DIRECT and ADJ/CONT parts.
//So it would often fall in the other part of the loop,
//WHERE THE MAGNITUDE OF THE CONTROL IS REDUCED A LOT!
//So it would not move any further
    
// if abs(delta) < abs (deltaold) try also that

  std::cout << "@@@@@@@@@@@@@@@@ It means that you are closer to the solution" << std::endl;
  //so solve the adjoint and control equations
  //these are LINEAR , so nonlinear loops are not required
  //nevertheless, coupling loops may be needed

//let us do a big LINEAR with PSEUDO-TIME STEPPING algorithm
//the only thing we have to carry on is the COUPLING.
//The idea is to have FEW OPTIMIZATION LOOPS (the outermost one)


{
   uint k_ADJCONT=0;

 do {
     
   k_ADJCONT++;
 
  std::cout << "\n @@@@@@@@@@@@@@@@ Solve ADJOINT and CONTROL altogether in pseudo-time stepping, " << k_ADJCONT << std::endl;

#ifdef NSAD_EQUATIONS
 {
   uint k_NSAD=0;
  do {
   k_NSAD++;
   lin_deltax_NSAD = MGTimeStep(k_NSAD,mgNSAD);
  }
  while (lin_deltax_NSAD > eps_NSAD && k_NSAD < MaxIterNSAD );
 }
#endif     
  #ifdef MHDAD_EQUATIONS
{  
  uint k_MHDAD=0;
 do{
   k_MHDAD++;
lin_deltax_MHDAD = MGTimeStep(k_MHDAD,mgMHDAD);
 }
while(lin_deltax_MHDAD >  eps_MHDAD && k_MHDAD < MaxIterMHDAD );
}
#endif
   #ifdef MHDCONT_EQUATIONS
{
    uint k_MHDCONT=0;
do{
  k_MHDCONT++;
  //not only when k_MHDCONT==1, but also when k_ADJCONT==1
lin_deltax_MHDCONT = MGTimeStep(k_MHDCONT,mgMHDCONT);
 }
while(lin_deltax_MHDCONT >  eps_MHDCONT && k_MHDCONT < MaxIterMHDCONT );
}
#endif
  
  } while( /* (lin_deltax_NSAD + lin_deltax_MHDAD + lin_deltax_MHDCONT) > eps_ADJCONT &&*/ k_ADJCONT< MaxIterADJCONT );
//even if the three systems are linear, their communication is not over right after A SINGLE OVERALL ADJCONT loop
//so remove the condition on the sum of the norms of the residual
}
  
  //update Boold = the previous Bold, in MGTimeStep
  
  //here,update omega
  omega = 1.5*omega;
//  omega=1.5;
    //once, I was not using this overrelaxation. Actually,I was only using the new Becont value.
    //Here, instead, it seems like I want to BOOST... maybe try without boosting first.
//in fact! you see that if you dont increase omega then it's more stable
//are we sure that Omega can be greater than one?



  //then this optimization step is over
     }

   else {
       std::cout << "@@@@@@@@@@@@@@@@ You went too far in following the gradient: reduce the intensity " << std::endl;
    if ( paral::get_rank() == 0 ){ 
            intgr_fstream.open(intgr_fname.c_str(),ios_base::app); 
            intgr_fstream << " J > Jold" << std::endl;
            intgr_fstream.close();  //you have to close to disassociate the file from the stream

    }
       
       
       
//// 
////   //you went too far in following the gradient
////   //update omega: reduce the intensity
   omega = 0.5*omega;
//   omega = 0.5;
   
////   
////   //this optimization step is over

}
   
  }
   else {//abs(J - Jold)

 std::cout << "@@@@@@@@@@@@@@@@ You possibly reached a local minimum " << std::endl;
//in this case you dont have to continue computing other optimization steps
//this is the convergence criterion of your optimization loop

 return; //you should also print... well, the ADJs and CONTR you printed in the previous step is exactly what it takes,
         // but you should reprint the STATE part

  }


 
//   update Jold, print, and restart with a new one
   
   std::cout << "@@@@@@@@@@@@@@@@ The optimization step is over " << std::endl;


    
//at the end of the optimization step, see the result for DIRECT and ADJOINT and CONTROL equations
//now we are printing OUTSIDE the nonlinear loop, so we do not print the nonlinear steps
const uint delta_opt_step = opt_step - _t_idx_in;
     if (delta_opt_step%print_step == 0) e_map_in.PrintSol(opt_step,pseudo_opttimeval);   //print sol.N.h5 and sol.N.xmf
  
      
    }

  
  return;
}



void OptLoop::init_equation_data(const SystemTwo* eqn) {
  
//======equation-specific vectors     =====================
      _x_oldopt.resize(eqn->_NoLevels);
      
        for(uint Level = 0; Level < _x_oldopt.size(); Level++)  {
   uint n_glob = eqn->_dofmap._Dim[Level]; //is it already filled? Now yes!!!!!!!!!
  _x_oldopt[Level] = NumericVector::build().release(); _x_oldopt[Level]->init(n_glob,false, SERIAL);
       }
 
  
 return; 
}




//===============================

//this function should belong to the NS equation, or to MultiLevelProblemTwo
//actually the best place is an OptimalControl framework

//the problem is that this function uses all the structures 
//for dofs and gauss which have been implemented in Eqn,
// therefore i cannot move it so easily;
// in this sense it cannot belong to the OptPhys class.
//Now the point is: what would have happened if the EqnNS
// was a LIBRARY class? We could not add the ComputeIntegral
// function to its member functions, because it is 
// an optimal-control related stuff.
//So, here's another reason why the EqnNS must be application related;
//but on the other hand here it is another reason
//why we have to think of making a NS equation MORE ABSTRACT,
//or at least some of its parts, in such a way that we can 
// SHARE AT LEAST SOME PARTS OF IT TO OTHER "NS FLAVOURS".
//We have to think of operators and boundary conditions,
// in a more general framework

//remember that the mesh is NON-DIMENSIONAL inside the code,
// so you should multiply all coordinates by Lref
 
//The value of this was 3 times in 3D and 2 times in 2D
//it means that a loop is done DIMENSION times instead of 1 time


double ComputeIntegral (const uint Level, const MultiLevelMeshTwo* mesh, const SystemTwo* eqn) {

   const uint mesh_vb = VV;
  
    CurrentElem       currelem(VV,eqn,*mesh,eqn->_eqnmap._elem_type);
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,eqn->_eqnmap, mesh->get_dim());
  
  // processor index
  const uint myproc = mesh->_iproc;
  // geometry -----
  const uint  space_dim =       mesh->get_dim();
  const uint   mesh_ord = (int) mesh->GetRuntimeMap().get("mesh_ord");  
  const uint     meshql = (int) mesh->GetRuntimeMap().get("meshql");    //======== ELEMENT MAPPING =======
  
//======Functions in the integrand ============
  
//========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = DIMENSION;
    xyz._FEord    = meshql;
    xyz._ndof     = eqn->_eqnmap._elem_type[currelem.GetDim()-1][xyz._FEord]->GetNDofs();
    xyz.Allocate();

//========== Quadratic domain, auxiliary  
  CurrentQuantity xyz_refbox(currgp);
  xyz_refbox._dim      = DIMENSION;
  xyz_refbox._FEord    = mesh_ord; //this must be QUADRATIC!!!
  xyz_refbox._ndof     = mesh->GetGeomEl(currelem.GetDim()-1,xyz_refbox._FEord)._elnds;
  xyz_refbox.Allocate();
  
     //========== 
    CurrentQuantity Vel(currgp);
    Vel._qtyptr      = eqn->_eqnmap._qtymap.get_qty("Qty_Velocity");
    Vel.VectWithQtyFillBasic();
    Vel.Allocate();
    
    //========== 
    CurrentQuantity VelDes(currgp);
    VelDes._qtyptr      = eqn->_eqnmap._qtymap.get_qty("Qty_DesVelocity");
    VelDes.VectWithQtyFillBasic();
    VelDes.Allocate();
   
   double integral=0.;
    
      const uint el_ngauss = eqn->_eqnmap._qrule[currelem.GetDim()-1].GetGaussPointsNumber();

//parallel sum
    const uint nel_e = mesh->_off_el[mesh_vb][mesh->_NoLevels*myproc+Level+1];
    const uint nel_b = mesh->_off_el[mesh_vb][mesh->_NoLevels*myproc+Level];
  
    for (uint iel=0; iel < (nel_e - nel_b); iel++) {

    currelem.set_el_nod_conn_lev_subd(Level,mesh->_iproc,iel);
    currelem.SetMidpoint();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);
    mesh->TransformElemNodesToRef(currelem.GetDim(),currelem.GetNodeCoords(),&xyz_refbox._val_dofs[0]);

//======= 
    xyz_refbox.SetElemAverage();
    int el_flagdom = ElFlagControl(xyz_refbox._el_average,mesh);
//=======        

    if ( Vel._eqnptr != NULL )       Vel.GetElemDofs(Level);
    else                             Vel._qtyptr->FunctionDof(Vel,0./*time*/,&xyz_refbox._val_dofs[0]);    //give the Hartmann flow, if not solving NS
    if ( VelDes._eqnptr != NULL ) VelDes.GetElemDofs(Level);
    else                          VelDes._qtyptr->FunctionDof(VelDes,0./*time*/,&xyz_refbox._val_dofs[0]);    

//AAA time is picked as a function pointer of the time C library i think...
    // it doesnt say it was not declared
    //here's why you should remove all unused headers always!
    

    for (uint qp = 0; qp < el_ngauss; qp++) {

for (uint fe = 0; fe < QL; fe++)     {  currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);  }  
     
   const double  Jac_g = currgp.JacVectVV_g(xyz);  //not xyz_refbox!      
   const double  wgt_g = eqn->_eqnmap._qrule[currelem.GetDim()-1].GetGaussWeight(qp);

     for (uint fe = 0; fe < QL; fe++)     {     currgp.SetPhiElDofsFEVB_g (fe,qp);  }

 Vel.val_g();
 VelDes.val_g();


  double deltau_squarenorm_g = 0.;
for (uint j=0; j<space_dim; j++) { deltau_squarenorm_g += (Vel._val_g[j] - VelDes._val_g[j])*(Vel._val_g[j] - VelDes._val_g[j]); }

  //NO for (uint j=0; j<space_dim; j++) { the integral is a scalar!
 
  integral += el_flagdom*wgt_g*Jac_g*deltau_squarenorm_g;

  //}
   
   
    }//gauss loop
     
    }//element loop
    
  return integral;  
  
}



  int ElFlagControl(const std::vector<double> el_xm, const MultiLevelMeshTwo* mesh)  {

  Box* box= static_cast<Box*>(mesh->GetDomain());
   
   
     int el_flagdom=0;

///optimal control
  #if DIMENSION==2
   //flag on the controlled region 2D
       if (   el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])
	   && el_xm[0] < 0.75*(box->_le[0] - box->_lb[0])
	   && el_xm[1] > 0.75*(box->_le[1] - box->_lb[1]) ) {
                 el_flagdom=1;
             }
  #else
   //flag on the controlled region 3D
      if ( el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])  
	&& el_xm[0] < 0.75*(box->_le[0] - box->_lb[0]) 
	&& el_xm[1] > 0.75*(box->_le[1] - box->_lb[1])
	&& el_xm[2] > 0.25*(box->_le[2] - box->_lb[2]) 
	&& el_xm[2] < 0.75*(box->_le[2] - box->_lb[2]) ) {
	el_flagdom=1;
        }
 #endif

return el_flagdom; 
}




} //end namespace femus



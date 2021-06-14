#include "OptLoop.hpp"


#include <cmath>

//library headers
#include "NumericVector.hpp"
#include "MultiLevelProblem.hpp"
#include "SystemTwo.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "TimeLoop.hpp"
#include "Files.hpp"
#include "CurrentElem.hpp"
#include "CurrentQuantity.hpp"
#include "Quantity.hpp"
#include "ElemType.hpp"
#include "Box.hpp"

#include "paral.hpp"
#include "XDMFWriter.hpp"

//application headers
#include "OptLoop.hpp"


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
//x_old and _xoold will be used
//either for OPTIMIZATION iterations or for NONLINEAR iterations,
//depending on what you need

 OptLoop::OptLoop(Files& files_in, const FemusInputParser<double> & map_in): TimeLoop(files_in,map_in) { }

 OptLoop::~OptLoop()  {

          delete _x_oldopt;

  }




void OptLoop::optimization_loop(MultiLevelProblem & e_map_in)  {


 #ifdef NS_EQUATIONS
   SystemTwo & eqnNS = e_map_in.get_system<SystemTwo>("Eqn_NS");
  #endif
  #ifdef MHD_EQUATIONS
   SystemTwo & eqnMHD = e_map_in.get_system<SystemTwo>("Eqn_MHD");
  #endif
 #ifdef NSAD_EQUATIONS
   SystemTwo & eqnNSAD = e_map_in.get_system<SystemTwo>("Eqn_NSAD");
 #endif
  #ifdef MHDAD_EQUATIONS
   SystemTwo & eqnMHDAD = e_map_in.get_system<SystemTwo>("Eqn_MHDAD");
  #endif
  #ifdef MHDCONT_EQUATIONS
   SystemTwo & eqnMHDCONT = e_map_in.get_system<SystemTwo>("Eqn_MHDCONT");
  #endif


std::string intgr_fname = _files.GetOutputPath() + "/" + "integral.txt";

 std::ofstream intgr_fstream;

//INITIALIZE OPT LOOP
//pseudo time parameters for optimization
    const uint      NoLevels = e_map_in.GetMeshTwo()._NoLevels;
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
     init_equation_data(&eqnMHDCONT);
#endif

//nitialize Becont
  _x_oldopt->zero();    //initialize Boldopt=0;


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

// // //     eqnMHDCONT._x_oold->close();
// // //     std::cout << "Linfty norm of Becont _xoold " << eqnMHDCONT._x_oold->linfty_norm() << std::endl;
// // //
// // //     eqnMHDCONT._x_old[NoLevels - 1]->close();
// // //     std::cout << "Linfty norm of Becont _x_old " << eqnMHDCONT._x_old [NoLevels - 1]->linfty_norm() << std::endl;
// // //
// // //     _x_oldopt->close();
// // //     std::cout << "Linfty norm of Becont _x_oldopt " << _x_oldopt->linfty_norm() << std::endl;

//////////////////
        NumericVector* _x_tmp2 = NumericVector::build().release();
        _x_tmp2->init(eqnMHDCONT._dofmap._Dim[eqnMHDCONT.GetGridn()-1],false, SERIAL);

      _x_tmp2->zero();
    *(_x_tmp2) = *(eqnMHDCONT._LinSolver[NoLevels-1]->_EPSC);
      eqnMHDCONT._bcond.Bc_ScaleDofVec(_x_tmp2, omega );
      _x_tmp2->close();
    std::cout << "Omega " << omega << std::endl;
    std::cout << "Linfty norm of Becont _x_old*omega "
              << _x_tmp2->linfty_norm() << std::endl;

// // //       eqnMHDCONT._x_oold->close();
// // //     std::cout << "Linfty norm of Becont _xoold "
// // //               << eqnMHDCONT._x_oold->linfty_norm() << std::endl;

     _x_oldopt->close();
    std::cout << "Linfty norm of Becont _x_oldopt "
              << _x_oldopt->linfty_norm() << std::endl;

      eqnMHDCONT._bcond.Bc_AddScaleDofVec(_x_oldopt,_x_tmp2,1.- omega);
      _x_tmp2->close();
    std::cout << "Linfty norm of Becont x_old*omega + (1-omega)*xoold "
              << _x_tmp2->linfty_norm() << std::endl;

    *(eqnMHDCONT._LinSolver[NoLevels-1]->_EPSC) = *(_x_tmp2);

    eqnMHDCONT._LinSolver[NoLevels-1]->_EPSC->close();
    std::cout << "Linfty norm of Becont _x_old updated "
              << eqnMHDCONT._LinSolver[NoLevels-1]->_EPSC->linfty_norm() << std::endl;


//   eqnMHDCONT.x[NoLevels - 1]->close();
//   std::cout << "Linfty norm of Becont x " << eqnMHDCONT.x[NoLevels - 1]->linfty_norm() << std::endl;
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
//      std::cout << "\n >>>>> Solving nonlinear step " << knonl_NS << " for" << eqnNS->_eqname << std::endl;
    nonlin_deltax_NS = MGTimeStep(knonl_NS,&eqnNS);

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
    nonlin_deltax_MHD = MGTimeStep(k_MHD,&eqnMHD);

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
   integral = ComputeIntegral(e_map_in.GetMeshTwo()._NoLevels - 1,&e_map_in.GetMeshTwo(),&eqnNS);
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
      intgr_fstream.open(intgr_fname.c_str(),std::ios_base::app);
      intgr_fstream << opt_step << " " << pseudo_opttimeval << " " << J << " " << omega /*<< std::endl*/;
      intgr_fstream.close();  //you have to close to disassociate the file from the stream
}

//********* see if J<J_old ********************
 //do with/without scaling (with/without ||1 see what happens)
if ( fabs(J - Jold) > epsJ /*|| 1*/  ) {
  std::cout << "@@@@@@@@@@@@@@@@ Jold= " << Jold << " , J= " << J << std::endl;

  if (J < Jold  /*|| 1*/ ) {

    if ( paral::get_rank() == 0 ){
            intgr_fstream.open(intgr_fname.c_str(),std::ios_base::app);
            intgr_fstream << "@@@@@@@@@@@@@@@@ J < Jold" << std::endl;
            intgr_fstream.close();  //you have to close to disassociate the file from the stream
        }
//******* update Jold  //you must update it only here, because here it is the good point to restart from
    Jold = J;
#ifdef MHDCONT_EQUATIONS
        *(_x_oldopt) = *(eqnMHDCONT._LinSolver[NoLevels-1]->_EPSC);
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
   lin_deltax_NSAD = MGTimeStep(k_NSAD,&eqnNSAD);
  }
  while (lin_deltax_NSAD > eps_NSAD && k_NSAD < MaxIterNSAD );
 }
#endif
  #ifdef MHDAD_EQUATIONS
{
  uint k_MHDAD=0;
 do{
   k_MHDAD++;
lin_deltax_MHDAD = MGTimeStep(k_MHDAD,&eqnMHDAD);
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
lin_deltax_MHDCONT = MGTimeStep(k_MHDCONT,&eqnMHDCONT);
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
            intgr_fstream.open(intgr_fname.c_str(),std::ios_base::app);
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
     if (delta_opt_step%print_step == 0) XDMFWriter::PrintSolLinear(_files.GetOutputPath(),opt_step,pseudo_opttimeval,e_map_in);   //print sol.N.h5 and sol.N.xmf


     delete _x_tmp2;  //delete vector

    }


  return;
}



void OptLoop::init_equation_data(const SystemTwo* eqn) {

   uint n_glob = eqn->_dofmap._Dim[eqn->GetGridn()-1];
  _x_oldopt = NumericVector::build().release();
  _x_oldopt->init(n_glob,false, SERIAL);

 return;
}




//===============================

//this function should belong to the NS equation, or to MultiLevelProblemTwo
//actually the best place is an OptimalControl framework

//remember that the mesh is NON-DIMENSIONAL inside the code,
// so you should multiply all coordinates by Lref

//The value of this was 3 times in 3D and 2 times in 2D
//it means that a loop is done DIMENSION times instead of 1 time

double ComputeIntegral (const uint Level, const MultiLevelMeshTwo* mesh, const SystemTwo* eqn) {

   const uint mesh_vb = VV;

  Mesh		*mymsh		=  eqn->GetMLProb()._ml_msh->GetLevel(Level);
  const unsigned myproc  = mymsh->processor_id();
  // geometry -----
  const uint  space_dim =       mesh->get_dim();


   double integral = 0.;

//parallel sum
    const uint nel_e = mesh->_off_el[mesh_vb][mesh->_NoLevels*myproc+Level+1];
    const uint nel_b = mesh->_off_el[mesh_vb][mesh->_NoLevels*myproc+Level];

    for (uint iel=0; iel < (nel_e - nel_b); iel++) {

    CurrentElem<double>       currelem(iel,myproc,Level,VV,eqn,*mesh,eqn->GetMLProb().GetElemType(),mymsh);
    //   CurrentGaussPointBase & currgp = //   CurrentGaussPointBase::build(currelem,eqn->GetMLProb().GetQuadratureRule(currelem.GetDim()));

//========= DOMAIN MAPPING
    CurrentQuantity xyz(currelem);
    xyz._dim      = DIMENSION;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

     //==========
    CurrentQuantity VelX(currelem);
    VelX._qtyptr      = eqn->GetMLProb().GetQtyMap().GetQuantity("Qty_Velocity0");
    VelX.VectWithQtyFillBasic();
    VelX.Allocate();

    CurrentQuantity VelY(currelem);
    VelY._qtyptr      = eqn->GetMLProb().GetQtyMap().GetQuantity("Qty_Velocity1");
    VelY.VectWithQtyFillBasic();
    VelY.Allocate();

    CurrentQuantity VelZ(currelem);
    VelZ._qtyptr      = eqn->GetMLProb().GetQtyMap().GetQuantity("Qty_Velocity2");
    VelZ.VectWithQtyFillBasic();
    VelZ.Allocate();

    std::vector<CurrentQuantity*> Vel_vec;
    Vel_vec.push_back(&VelX);
    Vel_vec.push_back(&VelY);
    Vel_vec.push_back(&VelZ);

    //==========
  CurrentQuantity VelDesX(currelem);
    VelDesX._dim        = 1;
    VelDesX._FEord      = VelX._FEord;
    VelDesX._ndof       = VelX._ndof;
    VelDesX.Allocate();

  CurrentQuantity VelDesY(currelem);
    VelDesY._dim        = 1;
    VelDesY._FEord      = VelX._FEord;
    VelDesY._ndof       = VelX._ndof;
    VelDesY.Allocate();

  CurrentQuantity VelDesZ(currelem);
    VelDesZ._dim        = 1;
    VelDesZ._FEord      = VelX._FEord;
    VelDesZ._ndof       = VelX._ndof;
    VelDesZ.Allocate();

  std::vector<CurrentQuantity*> VelDes_vec;
    VelDes_vec.push_back(&VelDesX);
    VelDes_vec.push_back(&VelDesY);
    VelDes_vec.push_back(&VelDesZ);


    currelem.SetDofobjConnCoords();


    currelem.ConvertElemCoordsToMappingOrd(xyz);
    xyz.SetElemAverage();
    int el_flagdom = ElFlagControl(xyz._el_average,eqn->GetMLProb()._ml_msh);


 for (uint idim=0; idim < space_dim; idim++)    {
       Vel_vec[idim]->GetElemDofs();
       VelDesired(&eqn->GetMLProb(),*VelDes_vec[idim],currelem,idim);
    }

      const uint el_ngauss = eqn->GetMLProb().GetQuadratureRule(currelem.GetDim()).GetGaussPointsNumber();

    for (uint qp = 0; qp < el_ngauss; qp++) {

for (uint fe = 0; fe < QL; fe++) {
//   currgp.SetPhiElDofsFEVB_g (fe,qp);
//   currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);
}

   const double  Jac_g = 1.; //currgp.JacVectVV_g(xyz);  //not xyz_refbox!
   const double  wgt_g = eqn->GetMLProb().GetQuadratureRule(currelem.GetDim()).GetGaussWeight(qp);


 for (uint idim=0; idim < space_dim; idim++) {
   Vel_vec[idim]->val_g();
   VelDes_vec[idim]->val_g();
 }

  double deltau_squarenorm_g = 0.;
for (uint j=0; j<space_dim; j++) { deltau_squarenorm_g += (Vel_vec[j]->_val_g[0] - VelDes_vec[j]->_val_g[0])*(Vel_vec[j]->_val_g[0] - VelDes_vec[j]->_val_g[0]); }


  integral += el_flagdom*wgt_g*Jac_g*deltau_squarenorm_g;


    }//gauss loop

    }//element loop

  return integral;

}



  int ElFlagControl(const std::vector<double> el_xm, const MultiLevelMesh* mesh)  {

  Box* box= static_cast<Box*>(mesh->GetDomain());


     int el_flagdom=0;

  if (mesh->GetDimension() == 2) {
   //flag on the controlled region 2D
       if (   el_xm[0] > 0.25*(box->_le[0] - box->_lb[0])
	   && el_xm[0] < 0.75*(box->_le[0] - box->_lb[0])
	   && el_xm[1] > 0.75*(box->_le[1] - box->_lb[1]) ) {
                 el_flagdom=1;
             }
  }
      else if (mesh->GetDimension() == 3) {
   //flag on the controlled region 3D
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



 void VelDesired(const MultiLevelProblem * ml_prob, CurrentQuantity& myvect, const CurrentElem<double> & currelem, const uint idim)  {

  const double Lref = ml_prob->GetInputParser().get("Lref");
  const double Uref = ml_prob->GetInputParser().get("Uref");
  double ILref = 1./Lref;

  const double rhof   = ml_prob->GetInputParser().get("rho0");
  const double muvel  = ml_prob->GetInputParser().get("mu0");
  const double MUMHD  = ml_prob->GetInputParser().get("MUMHD");
  const double SIGMHD = ml_prob->GetInputParser().get("SIGMHD");
  const double Bref   = ml_prob->GetInputParser().get("Bref");

  const double DpDz   = 1./*0.5*/;  //AAA: change it according to the pressure distribution

  double DpDzad = DpDz*Lref/(rhof*Uref*Uref);

  double Re  = ml_prob->GetInputParser().get("Re");
  double Rem = ml_prob->GetInputParser().get("Rem");
  double Hm  = ml_prob->GetInputParser().get("Hm");
  double S   = ml_prob->GetInputParser().get("S");


  Box* box= static_cast<Box*>(ml_prob->_ml_msh->GetDomain());


  double Lhalf = 0.5*(box->_le[0] - box->_lb[0]);
  double Lmid  = 0.5*(box->_le[0] + box->_lb[0]);



  //constant for the real reference length in the Hartmann number
  const double LHm =2.;   //this is because the reference length for Hm is HALF THE WIDTH of the domain, which is Lref=1 now



       for (uint d=0; d < myvect._ndof; d++)   {


          if (idim == 0) myvect._val_dofs[d] =  0.;
     else if (idim == 1){
         double xtr = currelem.GetNodeCoords()[d] - Lmid;
         const double magnitude = ml_prob->GetInputParser().get("udes")*DpDzad*Hm/LHm*(cosh(Hm/LHm) - cosh(Hm/LHm*xtr*Lref/Lhalf)) / (SIGMHD*Bref*Bref*sinh(Hm/LHm)*Uref);
                         myvect._val_dofs[d] =  magnitude;
     }
     else if (idim == 2) myvect._val_dofs[d] =  0.;


      }

  return;

 }






//---------------------------------------------------------------------------------------------------------------------

double SetInitialCondition(const MultiLevelProblem * ml_prob, const std::vector <double> &xp, const char name[]) {

//   std::vector<double> xp(ml_prob->_ml_msh->GetDimension());
//   xp[0] = x;
//   xp[1] = y;
//   xp[2] = z;
//
//   if ( ml_prob->_ml_msh->GetDimension() == 3 )    xp[2] = z;

  // defaults ***********
  double value = 0.;
  // defaults ***********

  const double bdry_toll = DEFAULT_BDRY_TOLL;

  Box* box = static_cast<Box*>(ml_prob->_ml_msh->GetDomain());

  std::vector<double> lb(ml_prob->_ml_msh->GetDimension());
  std::vector<double> le(ml_prob->_ml_msh->GetDimension());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (ml_prob->_ml_msh->GetDimension() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }

    std::vector<double> x_rotshift(ml_prob->_ml_msh->GetDimension());
  ml_prob->_ml_msh->GetDomain()->TransformPointToRef(&xp[0],&x_rotshift[0]);






//=================
   if(!strcmp(name,"Qty_Velocity0")) {

      value = 0.;

  }

  else if(!strcmp(name,"Qty_Velocity1")) {

  const double Uref = ml_prob->GetInputParser().get("Uref");

     value = 0.;

   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  {  value = (x_rotshift[0] - lb[0])*(le[0] - x_rotshift[0])*(x_rotshift[2] - lb[2])*(le[2]-x_rotshift[2])/Uref;  }


  }

  else if(!strcmp(name,"Qty_Velocity2")) {

      value = 0.;

  }


//=================
  else if(!strcmp(name,"Qty_VelocityAdj0")) {

      value = 0.;

  }

  else if(!strcmp(name,"Qty_VelocityAdj1")) {

       value = 0.;

  }

  else if(!strcmp(name,"Qty_VelocityAdj2")) {

         value = 0.;

  }

 //=================
  else if(!strcmp(name,"Qty_MagnFieldHomAdj0")) {

      value = 0.;

  }

  else if(!strcmp(name,"Qty_MagnFieldHomAdj1")) {

      value = 0.;

  }

  else if(!strcmp(name,"Qty_MagnFieldHomAdj2")) {

      value = 0.;

  }

 //=================
  else if(!strcmp(name,"Qty_MagnFieldHom0")) {

      value = 0.;

  }

  else if(!strcmp(name,"Qty_MagnFieldHom1")) {

      value = 0.;

  }

  else if(!strcmp(name,"Qty_MagnFieldHom2")) {

      value = 0.;

  }


 //=================
  else if(!strcmp(name,"Qty_MagnFieldExt0")) {

  const double Bref = ml_prob->GetInputParser().get("Bref");

  value = Bref/Bref;

  }

  else if(!strcmp(name,"Qty_MagnFieldExt1")) {

      value = 0.;

  }

  else if(!strcmp(name,"Qty_MagnFieldExt2")) {

      value = 0.;

  }

   //=================
   //=================
   //=================

  else if(!strcmp(name,"Qty_Pressure")) {

      value = 0.;

  }


  else if(!strcmp(name,"Qty_PressureAdj")) {

      value = 0.;

  }

  else if(!strcmp(name,"Qty_MagnFieldHomLagMult")) {

       value = 0.;

  }


  else if(!strcmp(name,"Qty_MagnFieldExtLagMult")) {

      value = 0.;

  }

  else if(!strcmp(name,"Qty_MagnFieldHomLagMultAdj")) {

       value = 0.;

  }

    //=================
   //=================
   //=================
 //=================
  else if(!strcmp(name,"Qty_DesVelocity0")) {

      value = 0.;

  }

  else if(!strcmp(name,"Qty_DesVelocity1")) {

  const double Lref = ml_prob->GetInputParser().get("Lref");
  const double Uref = ml_prob->GetInputParser().get("Uref");
  double ILref = 1./Lref;

  const double rhof   = ml_prob->GetInputParser().get("rho0");
  const double muvel  = ml_prob->GetInputParser().get("mu0");
  const double MUMHD  = ml_prob->GetInputParser().get("MUMHD");
  const double SIGMHD = ml_prob->GetInputParser().get("SIGMHD");
  const double Bref   = ml_prob->GetInputParser().get("Bref");

  const double DpDz   = 1./*0.5*/;  //AAA: change it according to the pressure distribution

  double DpDzad = DpDz*Lref/(rhof*Uref*Uref);

  double Re  = ml_prob->GetInputParser().get("Re");
  double Rem = ml_prob->GetInputParser().get("Rem");
  double Hm  = ml_prob->GetInputParser().get("Hm");
  double S   = ml_prob->GetInputParser().get("S");


  Box* box= static_cast<Box*>(ml_prob->_ml_msh->GetDomain());


  double Lhalf = 0.5*(box->_le[0] - box->_lb[0]);
  double Lmid  = 0.5*(box->_le[0] + box->_lb[0]);

  double xtr = xp[0] - Lmid;


  //constant for the real reference length in the Hartmann number
  const double LHm =2.;   //this is because the reference length for Hm is HALF THE WIDTH of the domain, which is Lref=1 now

  const double magnitude = ml_prob->GetInputParser().get("udes")*DpDzad*Hm/LHm*(cosh(Hm/LHm) - cosh(Hm/LHm*xtr*Lref/Lhalf)) / (SIGMHD*Bref*Bref*sinh(Hm/LHm)*Uref);

  value = magnitude;

  }

  else if(!strcmp(name,"Qty_DesVelocity2")) {

      value = 0.;

  }



  return value;
}

// =====================================================================

bool SetBoundaryCondition(const MultiLevelProblem * ml_prob, const std::vector <double> &xp, const char name[], double &value, const int facename, const double time) {

//   std::vector<double> xp(ml_prob->_ml_msh->GetDimension());
//   xp[0] = x;
//   xp[1] = y;
//
//   if ( ml_prob->_ml_msh->GetDimension() == 3 )    xp[2] = z;

  // defaults ***********
  bool test=1; //dirichlet  // 0 neumann
  value=0.;
  // defaults ***********

  const double bdry_toll = DEFAULT_BDRY_TOLL;

  Box* box = static_cast<Box*>(ml_prob->_ml_msh->GetDomain());

  std::vector<double> lb(ml_prob->_ml_msh->GetDimension());
  std::vector<double> le(ml_prob->_ml_msh->GetDimension());
  lb[0] = box->_lb[0]; //already nondimensionalized
  le[0] = box->_le[0];
  lb[1] = box->_lb[1];
  le[1] = box->_le[1];
  if (ml_prob->_ml_msh->GetDimension() == 3) {
  lb[2] = box->_lb[2];
  le[2] = box->_le[2];
  }

    std::vector<double> x_rotshift(ml_prob->_ml_msh->GetDimension());
  ml_prob->_ml_msh->GetDomain()->TransformPointToRef(&xp[0],&x_rotshift[0]);


  //=================
   if(!strcmp(name,"Qty_Velocity0")) {

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )    test = 1;    //u x n
 if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  test = 1;    //u x n
 if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                test = 1;    //u dot n
 if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  test = 1;    //u dot n

      value = 0.;

  }

  else if(!strcmp(name,"Qty_Velocity1")) {

  const double Uref = ml_prob->GetInputParser().get("Uref");
       value = 0.;

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                     test = 1;    //u x n
  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )     test = 1;   //u x n
  if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )  {                                  test = 0;  //u dot n   //leave this free for VELOCITY INLET
                                                                                                     value = (x_rotshift[0] - lb[0])*(le[0] - x_rotshift[0])*(x_rotshift[2] - lb[2])*(le[2]-x_rotshift[2])/Uref;
  }
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )   test = 0;  //u dot n   //leave this free for outlet
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                 test = 1;   //u x n
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )   test = 1;   //u x n


  }

  else if(!strcmp(name,"Qty_Velocity2")) {

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )    test = 1;    //u x n
 if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                    test = 1;      //u x n
 if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  test = 1;     //u x n
 if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                test = 1;   //u dot n
 if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  test = 1; //u dot n

       value = 0.;

  }


//=================
  else if(!strcmp(name,"Qty_VelocityAdj0")) {

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )    test = 1;    //u x n
 if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                    test = 1;      //u x n
 if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  test = 1;     //u x n
 if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                test = 1;   //u dot n
 if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  test = 1; //u dot n

      value = 0.;

  }

  else if(!strcmp(name,"Qty_VelocityAdj1")) {

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                    test = 1;
  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll)     test = 1;
  if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)                                     test = 0; //INSTEAD THIS MUST CORRESPOND TO THE DIRECT
                                                                                                               //the SPACE of ADJOINT functions is the same as the SPACE for the DIRECT test functions
                                                                                                               //if you fix this then you dont control...
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)    test = 0;   //INSTEAD THIS MUST CORRESPOND TO THE DIRECT
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                test = 1;  //u x n
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)    test = 1;  //u x n

       value = 0.;

  }

  else if(!strcmp(name,"Qty_VelocityAdj2")) {

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )    test = 1;    //u x n
 if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                    test = 1;      //u x n
 if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  test = 1;     //u x n
 if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                test = 1;   //u dot n
 if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  test = 1; //u dot n

         value = 0.;

  }

 //=================
  else if(!strcmp(name,"Qty_MagnFieldHomAdj0")) {

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )    test = 1;    //u x n
 if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                    test = 1;      //u x n
 if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  test = 1;     //u x n
 if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                test = 1;   //u dot n
 if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  test = 1; //u dot n

      value = 0.;

  }

  else if(!strcmp(name,"Qty_MagnFieldHomAdj1")) {

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )    test = 1;    //u x n
 if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                    test = 0;      //u x n
 if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  test = 0;     //u x n
 if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                test = 1;   //u dot n
 if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  test = 1; //u dot n

      value = 0.;

  }

  else if(!strcmp(name,"Qty_MagnFieldHomAdj2")) {

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )    test = 1;    //u x n
 if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  test = 1;    //u x n
 if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                test = 0;    //u dot n
 if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  test = 0;    //u dot n

      value = 0.;

  }

 //=================
  else if(!strcmp(name,"Qty_MagnFieldHom0")) {

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )    test = 1;    //u x n
 if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                    test = 1;      //u x n
 if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  test = 1;     //u x n
 if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                test = 1;   //u dot n
 if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  test = 1; //u dot n

 value = 0.;

  }

  else if(!strcmp(name,"Qty_MagnFieldHom1")) {

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )    test = 1;    //u x n
 if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                    test = 0;      //u x n
 if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  test = 0;     //u x n
 if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                test = 1;   //u dot n
 if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  test = 1; //u dot n

      value = 0.;

  }

  else if(!strcmp(name,"Qty_MagnFieldHom2")) {

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )    test = 1;    //u x n
 if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                    test = 0;    //why isn't this fixed as well?
 if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  test = 0;    //why isn't this fixed as well?
 if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                test = 0;   //u dot n
 if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  test = 0; //u dot n
      value = 0.;

  }


 //=================
  else if(!strcmp(name,"Qty_MagnFieldExt0")) {

  const double Bref = ml_prob->GetInputParser().get("Bref");

  value = Bref/Bref;

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                  test = 1;
  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll)   test = 1;
  if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)                                   test = 1;
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  test = 0; //UseControl
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                              test = 1;
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  test = 1;

  }

  else if(!strcmp(name,"Qty_MagnFieldExt1")) {

  value = 0.;

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                  test = 1;
  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll)   test = 1;
  if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll)                                   test = 1;
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  test = 0; //UseControl
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                              test = 1;
  if ((le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll)  test = 1;

  }

  else if(!strcmp(name,"Qty_MagnFieldExt2")) {

 if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )    test = 1;    //u x n
 if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                    test = 1;    //u x n
 if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )  test = 1;    //u x n
 if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                test = 1;   //u dot n
 if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )  test = 1; //u dot n

 value = 0.;

  }

   //=================
   //=================
   //=================

  else if(!strcmp(name,"Qty_Pressure")) {


  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                     test = 0;  //don't do the integral
  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )     test = 0;
  if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                     test = 0;
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )   test = 1;  //do the integral
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                 test = 0;
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )   test = 0;

      value = 0.;


  }


  else if(!strcmp(name,"Qty_PressureAdj")) {

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                     test = 0;
  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )     test = 0;
  if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                     test = 0;
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )   test = 0;
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                 test = 0;
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )   test = 0;

      value = 0.;

  }

  else if(!strcmp(name,"Qty_MagnFieldHomLagMult")) {

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                     test = 0;
  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )     test = 0;
  if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                     test = 0;
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )   test = 0;
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                 test = 0;
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )   test = 0;

       value = 0.;


  }


  else if(!strcmp(name,"Qty_MagnFieldExtLagMult")) {

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                     test = 0;
  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )     test = 0;
  if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                     test = 0;
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )   test = 0;
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                 test = 0;
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )   test = 0;

      value = 0.;


  }

  else if(!strcmp(name,"Qty_MagnFieldHomLagMultAdj")) {

  if ( x_rotshift[0] > -bdry_toll &&  x_rotshift[0] < bdry_toll )                                     test = 0;
  if ( (le[0]-lb[0])  - x_rotshift[0] > -bdry_toll && (le[0]-lb[0]) - x_rotshift[0] < bdry_toll )     test = 0;
  if ( x_rotshift[1] > -bdry_toll &&  x_rotshift[1] < bdry_toll )                                     test = 0;
  if ( (le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll )   test = 0;
  if ( (x_rotshift[2]) > -bdry_toll && ( x_rotshift[2]) < bdry_toll )                                 test = 0;
  if ( (le[2]-lb[2]) -(x_rotshift[2]) > -bdry_toll &&  (le[2]-lb[2]) -(x_rotshift[2]) < bdry_toll )   test = 0;

       value = 0.;

  }


    return test;

}










} //end namespace femus



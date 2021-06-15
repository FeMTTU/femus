#include "OptLoop.hpp"

#include "MultiLevelProblem.hpp"
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
    void OptLoop::optimization_loop(MultiLevelProblem & ml_prob)  {

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
      OneTimestepEqnLoop(delta_t_step,ml_prob);

#if DEFAULT_PRINT_TIME==1
      std::clock_t    end_time=std::clock();
#endif

      // print solution
      if (delta_t_step%print_step == 0) XDMFWriter::PrintSolLinear(_files.GetOutputPath(),curr_step,curr_time,ml_prob);   //print sol.N.h5 and sol.N.xmf


#if DEFAULT_PRINT_TIME==1
      std::clock_t    end_time2=std::clock();
      std::cout << " Time solver ----->= "   << double(end_time- start_time)/ CLOCKS_PER_SEC
                << " Time printing ----->= " << double(end_time2- end_time) / CLOCKS_PER_SEC <<
                std::endl;
#endif

//=====functional evaluations=======
   SystemTwo & eqnT = ml_prob.get_system<SystemTwo>("Eqn_T");


     double J = 0.;
J = ComputeIntegral    ( ml_prob.GetMeshTwo()._NoLevels - 1,&ml_prob.GetMeshTwo(),&eqnT,_files.GetOutputTime());
J = ComputeNormControl ( ml_prob.GetMeshTwo()._NoLevels - 1,&ml_prob.GetMeshTwo(),&eqnT,0 );
J = ComputeNormControl ( ml_prob.GetMeshTwo()._NoLevels - 1,&ml_prob.GetMeshTwo(),&eqnT,1 );
//=====functional evaluations =======


    }   // end time loop

    return;
  }



// This function computes the integral only for the current processor

double ComputeIntegral (const uint Level, const MultiLevelMeshTwo* mesh, const SystemTwo* eqn, const std::string output_time) {

  Mesh		*mymsh		=  eqn->GetMLProb()._ml_msh->GetLevel(Level);
  //====== processor index
  const uint myproc = mesh->_iproc;
  const uint space_dim =      mesh->get_dim();

  const uint mesh_vb = VV;


   double integral = 0.;


//parallel sum
    const uint nel_e = mesh->_off_el[mesh_vb][mesh->_NoLevels*myproc+Level+1];
    const uint nel_b = mesh->_off_el[mesh_vb][mesh->_NoLevels*myproc+Level];

    for (uint iel=0; iel < (nel_e - nel_b); iel++) {

    CurrentElem<double>       currelem(iel,myproc,Level,VV,eqn,*mesh,eqn->GetMLProb().GetElemType(),mymsh);
    //   CurrentGaussPointBase & currgp = //   CurrentGaussPointBase::build(currelem,eqn->GetMLProb().GetQuadratureRule(currelem.GetDim()));

  //==========
    CurrentQuantity Tempold(currelem);
    Tempold._SolName = "Qty_Temperature";
    Tempold._qtyptr   =  eqn->GetMLProb().GetQtyMap().GetQuantity("Qty_Temperature");
    Tempold.VectWithQtyFillBasic();
    Tempold.Allocate();

  //==========
    CurrentQuantity Tlift(currelem);
    Tlift._SolName = "Qty_TempLift";
    Tlift._qtyptr   =  eqn->GetMLProb().GetQtyMap().GetQuantity("Qty_TempLift");
    Tlift.VectWithQtyFillBasic();
    Tlift.Allocate();

 //===========
    CurrentQuantity Tdes(currelem);
    Tdes._SolName = "Qty_TempDes";
    Tdes._dim      = Tempold._dim;
    Tdes._FEord    = Tempold._FEord;
    Tdes._ndof     = Tempold._ndof;
    Tdes.Allocate();

//========= DOMAIN MAPPING
    CurrentQuantity xyz(currelem);
    xyz._dim      = space_dim;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

      currelem.SetDofobjConnCoords();

     currelem.ConvertElemCoordsToMappingOrd(xyz);
// ===============
     xyz.SetElemAverage();
     int el_flagdom = ElFlagControl(xyz._el_average,eqn->GetMLProb()._ml_msh);
//====================

    Tempold.GetElemDofs();
    Tlift.GetElemDofs();

    TempDesired(Tdes,currelem);

    const uint el_ngauss = eqn->GetMLProb().GetQuadratureRule(currelem.GetDim()).GetGaussPointsNumber();

    for (uint qp = 0; qp < el_ngauss; qp++) {

     for (uint fe = 0; fe < QL; fe++)     {
//        currgp.SetPhiElDofsFEVB_g (fe,qp);
//        currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);
    }

   const double  Jac_g = 1.;// currgp.JacVectVV_g(xyz);
   const double  wgt_g = eqn->GetMLProb().GetQuadratureRule(currelem.GetDim()).GetGaussWeight(qp);

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

  const uint mesh_vb = VV;

  Mesh		*mymsh		=  eqn->GetMLProb()._ml_msh->GetLevel(Level);



   double integral = 0.;


//parallel sum
    const uint nel_e = mesh->_off_el[mesh_vb][mesh->_NoLevels*myproc+Level+1];
    const uint nel_b = mesh->_off_el[mesh_vb][mesh->_NoLevels*myproc+Level];

    for (int iel=0; iel < (nel_e - nel_b); iel++) {

    CurrentElem<double>       currelem(iel,myproc,Level,VV,eqn,*mesh,eqn->GetMLProb().GetElemType(),mymsh);
    //   CurrentGaussPointBase & currgp = //   CurrentGaussPointBase::build(currelem,eqn->GetMLProb().GetQuadratureRule(currelem.GetDim()));

//======Functions in the integrand ============

      //==========
    CurrentQuantity Tlift(currelem);
    Tlift._qtyptr   =  eqn->GetMLProb().GetQtyMap().GetQuantity("Qty_TempLift");
    Tlift._SolName = "Qty_TempLift";
    Tlift.VectWithQtyFillBasic();
    Tlift.Allocate();

//========= DOMAIN MAPPING
    CurrentQuantity xyz(currelem);
    xyz._dim      = space_dim;
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();


      currelem.SetDofobjConnCoords();

      currelem.ConvertElemCoordsToMappingOrd(xyz);

     Tlift.GetElemDofs();

     const uint el_ngauss = eqn->GetMLProb().GetQuadratureRule(currelem.GetDim()).GetGaussPointsNumber();

  for (uint qp = 0; qp < el_ngauss; qp++) {

     for (uint fe = 0; fe < QL; fe++)     {
//        currgp.SetPhiElDofsFEVB_g (fe,qp);
//        currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);
    }

      const double  Jac_g = 1.;//currgp.JacVectVV_g(xyz);
      const double  wgt_g = eqn->GetMLProb().GetQuadratureRule(currelem.GetDim()).GetGaussWeight(qp);

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

 int ElFlagControl(const std::vector<double> el_xm, const MultiLevelMesh* mesh)  {

  Box* box= static_cast<Box*>(mesh->GetDomain());

     int el_flagdom=0;

     if (mesh->GetDimension() == 2) {
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

 else if (mesh->GetDimension() == 3) {

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


 void TempDesired(CurrentQuantity& myvect, const CurrentElem<double> & currelem)  {

   for (uint ivar=0; ivar < myvect._dim; ivar++)
       for (uint d=0; d < myvect._ndof; d++)      myvect._val_dofs[d+ivar*myvect._ndof] =  0.9;


  return;

 }





//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryCondition(const MultiLevelProblem * ml_prob, const std::vector <double> &xp, const char name[], double &value, const int facename, const double time) {
/*
  std::vector<double> xp(ml_prob->_ml_msh->GetDimension());
  xp[0] = x;
  xp[1] = y;

  if ( ml_prob->_ml_msh->GetDimension() == 3 )    xp[2] = z;*/

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


  if(!strcmp(name,"Qty_Temperature")) {

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll )       test=1;
  if ( (le[0]-lb[0])  - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll)        test=1;
  if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)                               test=1;
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)           test=1;

  value = 0.;

  }


  else if(!strcmp(name,"Qty_TempLift")){

     value=1.;

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {  //left of the RefBox
       test=1;
      value=1.;
  }

 if ( (le[0]-lb[0])  - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll)  { //right of the RefBox
       test=1;
      value=1.;
   }


   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox

       if  ( (x_rotshift[0]) < 0.25*(le[0] - lb[0]) || ( x_rotshift[0]) > 0.75*(le[0] - lb[0]) )  {
       test=1;
      value=1.;

      }
       else {
       test=0; /// @todo what is the meaning of value for the Neumann case?
      }

     }

  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox
       test=1;
      value=0.;
    }


  }

  else if(!strcmp(name,"Qty_TempAdj")) {


  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll )                                     test = 1;
  if ( (le[0]-lb[0])  - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll)      test = 1;
  if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)                                      test = 1;
  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)         test = 1;

    value=0.;

  }


  else if(!strcmp(name,"Qty_Velocity0")){

  if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
       test=1;
      value=0.;

       if ( (x_rotshift[1]) > 0.4*(le[1] - lb[1]) && ( x_rotshift[1]) < 0.6*(le[1]-lb[1]) )  {
       value = ml_prob->GetInputParser().get("injsuc");
      }


  }

 if ( (le[0]-lb[0]) - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll  ){ //right of the RefBox
       test=1;
      value=0.;
  }

   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
       test=1;
      value=0.;
  }

  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox
       test=1;
      value=0.;
  } //top RefBox

  }

  else if(!strcmp(name,"Qty_Velocity1")){



      if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
       test=1;
      value=0.;
  }

 if ( (le[0]-lb[0]) - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll  ){ //right of the RefBox
       test=1;
      value=0.;
  }

 if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
        test=1;

      //below, inlet
 if ( (x_rotshift[0]) < 0.25*(le[0] - lb[0]) || ( x_rotshift[0]) > 0.75*(le[0] - lb[0]) ) {
    value = 0.;
  }
  else {
    value = 1.;
    }

  }

  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the RefBox

//outflow only on part of the outlet
  if ( (x_rotshift[0]) > 0.70*(le[0]-lb[0]) ) {
        test=0;
 }  //end part outflow
    else {
       test=1;
      value=0.;
    }

  } //top RefBox




  }

  else if(!strcmp(name,"Qty_Pressure")) {


     value =  1./ ml_prob->GetInputParser().get("pref")*( (le[1] - lb[1]) - xp[1] );




      if ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {//left of the RefBox
        test=0;
  }

 if ( (le[0]-lb[0]) - (x_rotshift[0]) > -bdry_toll && (le[0]-lb[0])  -(x_rotshift[0]) < bdry_toll  ){ //right of the RefBox
        test=0;
  }

   if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox
        test=0;
  }

  if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {  //top of the  of the RefBox

//outflow only on part of the outlet
  if ( (x_rotshift[0]) > 0.70*(le[0]-lb[0]) ) {
        test=1;   /// @todo this does not mean that there is a DIRICHLET condition on PRESSURE... it means DO THE INTEGRAL
 }  //end part outflow
    else {
        test=0;
    }

  } //top RefBox


  }

  return test;
}


//---------------------------------------------------------------------------------------------------------------------

double SetInitialCondition(const MultiLevelProblem * ml_prob, const std::vector <double> &xp, const char name[]) {

//   std::vector<double> xp(ml_prob->_ml_msh->GetDimension());
//   xp[0] = x;
//   xp[1] = y;
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


  if(!strcmp(name,"Qty_Temperature")) {

    value = 0.;
  }


  else if(!strcmp(name,"Qty_TempLift")){

    value = 1.;

      if ((le[1]-lb[1]) -(x_rotshift[1]) > -bdry_toll &&  (le[1]-lb[1]) -(x_rotshift[1]) < bdry_toll)  {
    value = 0.;
  }

  }



  else if(!strcmp(name,"Qty_TempAdj")){

      value = 0.;

  }


  else if(!strcmp(name,"Qty_Velocity0")){

      value = 0.;

     if  ( (x_rotshift[0]) > -bdry_toll && ( x_rotshift[0]) < bdry_toll ) {

 if ( (x_rotshift[1]) > 0.4*(le[1] - lb[1]) && ( x_rotshift[1]) < 0.6*(le[1]-lb[1]) )  {  //left of the refbox
       value = ml_prob->GetInputParser().get("injsuc");
      }
   }

  }

  else if(!strcmp(name,"Qty_Velocity1")){

     //==================================
    if (( x_rotshift[1]) > -bdry_toll && ( x_rotshift[1]) < bdry_toll)  { //bottom  of the RefBox

//below, inlet
  if ( (x_rotshift[0]) < 0.25*(le[0] - lb[0]) || ( x_rotshift[0]) > 0.75*(le[0] - lb[0]) ) {
    value = 0.;
  }
  else {
    value = 1.;
    }

  }


  }

  else if(!strcmp(name,"Qty_Pressure")) {


  value =  1./ ml_prob->GetInputParser().get("pref")*( (le[1] - lb[1]) - xp[1] );


  }

  else if(!strcmp(name,"Qty_TempDes")) {


  value =  0.9;


  }

  return value;
}



} //end namespace femus




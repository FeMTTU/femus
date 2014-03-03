//CONTROL PROGRAM
// launch it with sthg like 'mpiexec -n 6 main-dbg --dt 1 --alphaVel 55 --udes 3 --Bref 0.001'

//C++ includes ============
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <cstdlib>
#include <sstream>

// External library include ( LibMesh, PETSc...) ------------------------------
#include "FemusExtLib_conf.hpp"


#ifdef FEMUS_HAVE_PETSC// Petsc   //TODO remove it later, this is here only for the LOG at the end
#include "petsc.h"
#endif

// library includes
#include "paral.hpp"
#include "FemusDefault.hpp"
#include "FemusInit.hpp"
#include "Utils.hpp"
#include "Files.hpp"
#include "Physics.hpp"
#include "GeomEl.hpp"
#include "mesh.hpp"
#include "FETypeEnum.hpp"
#include "FEElemBase.hpp"
#include "QRule.hpp"
#include "EquationsMap.hpp"
#include "TimeLoop.hpp"
#include "Typedefs.hpp"
#include "CmdLine.hpp"
#include "Quantity.hpp"
#include "QTYnumEnum.hpp"
#include "Box.hpp"  //for the DOMAIN


// application includes ==============================
#include "Opt_conf.hpp"
#include "OptQuantities.hpp"
#include "OptPhysics.hpp"
#include "EqnNS.hpp"
#include "EqnNSAD.hpp"
#include "EqnMHD.hpp"
#include "EqnMHDAD.hpp"
#include "EqnMHDCONT.hpp"



//***************** functions for this application
void optimization_loop(EquationsMap& e_map_in);

double funzione(double t , const double* xyz) {return 1.;} 

// =======================================
// Main program
// =======================================

int main(int argc, char** argv) {

  // ====== FemusInit =====  //put this as the first call because mpi is initialized here
#ifdef LM_INIT
  LibMeshInit init(argc, argv);
#else
  FemusInit init(argc,argv);
#endif

   std::cout << "***** REMEMBER TO PASS THE VALUES TO COMMAND LINE ****** " << std::endl;

  CmdLine::parse(argc,argv);
  
// ======= Files ========================
  Files files; 
  files.get_frtmap().read();
  files.CheckIODirectories();

 //>>>>>>>>> REDIRECT COUT
   std::ofstream file; //if a filestream dies, then also its stream-buffer dies ?!? //So I have to declare it outside? Yes. This seems to work.
   std::streambuf* sbuf = std::cout.rdbuf();  //get the current buffer for cout
   files.RedirectCout(sbuf,file);
// >>>>>>>>>>>>> END REDIRECT COUT

   files.CopyGencaseFiles();
   files.PrintRun(DEFAULT_NEW_RUN);
   files.InitCaseData();

 // =========================================
  // ======= END OF THE INITIALIZATION PART ========================
  // =========================================
 
  // ======= Utils ========================
  Utils utils(files);
  utils._urtmap.read();  
// // // //======== NOW WE DON'T USE THE COMMAND LINE FOR THE VALUES, but the CONFIG ==============
// // // //======== SUBSTITUTE WITH COMMAND LINE VALUES==== modify the just read Utils map==============
// // //     utils.get_utils_map().erase("dt");
// // //     //you must erase before inserting the modified value
// // //     //erase doesnt give error if the parameter is not there
// // //  //PATH TO THE ERROR //when,I get_utils_map , it must be a reference!   
// // //     utils.set_par("dt",CmdLine::get("--dt"));
// // // //======== END SUBSTITUTE WITH COMMAND LINE VALUES================
  utils._urtmap.print();
      
  OptPhysics phys(utils);
  phys._physrtmap.read();
// // // //======== NOW WE DON'T USE THE COMMAND LINE FOR THE VALUES, but the CONFIG ==============
// // // //======== SUBSTITUTE WITH COMMAND LINE VALUES===== modify the just read Phys map ===========
// // //      phys.get_phys_map().erase("alphaVel");
// // //      phys.set_par("alphaVel",CmdLine::get("--alphaVel"));
// // //      phys.get_phys_map().erase("udes");
// // //      phys.set_par("udes",CmdLine::get("--udes"));
// // //      phys.get_phys_map().erase("Bref");
// // //      phys.set_par("Bref",CmdLine::get("--Bref"));
// // // 
// // //      std::cout <<"FROM COMMAND LINE dt " << timeloop........("dt") << std::endl;
// // //      std::cout <<"FROM COMMAND LINE AlphaVel " << phys.get_par("alphaVel") << std::endl;
// // //      std::cout <<"FROM COMMAND LINE udes " << phys.get_par("udes") << std::endl;
// // //      std::cout <<"FROM COMMAND LINE Bref " << phys.get_par("Bref") << std::endl;
// // // //======== END SUBSTITUTE WITH COMMAND LINE VALUES================

 phys.set_nondimgroups();
 phys._physrtmap.print();

//===========================================
  const double Lref  = phys._physrtmap.get("Lref");     // reference L
  Box mybox(utils);
      mybox._boxrtmap.read();
      mybox._boxrtmap.print();
      mybox.init(Lref);

// ====== GeomEl ================================
// ======  Mesh ================================
  uint geomel_type = (uint) utils._urtmap.get("geomel_type");
  uint dimension   = (uint) utils._urtmap.get("dimension");
  GeomEl geomel(dimension,geomel_type);
  Mesh mesh(utils,geomel,Lref,&mybox); 
  mesh.PrintForVisualizationAllLEVAllVB();

  phys.set_mesh(&mesh);

// ======  QRule ================================ //so far we have ONLY ONE quadrature rule for all the equations
  QRule   qrule(&geomel);
  
  // =======FEElems =====  //remember to delete the FE at the end
  std::vector<FEElemBase*> FEElements(QL);
 
  for (int fe=0; fe<QL; fe++) {
    FEElements[fe] = FEElemBase::build(&geomel,fe);       /*VB based*/  //The order of the fe is established by the library
//sort of constructor
    FEElements[fe]->SetOrder(fe);
    FEElements[fe]->AssociateQRule(&qrule);
    FEElements[fe]->SetUtils(&utils);
//end sort of constructor
    FEElements[fe]->init();
  }
  
  // ===== QuantityMap =========================================
  QuantityMap  qty_map(utils,phys);

//================================
// ======= Add QUANTITIES ========  
//================================
  MagnFieldHom bhom("Qty_MagnFieldHom",qty_map,dimension,FE_MAGNFIELDHOM);     qty_map.set_qty(&bhom);  
  MagnFieldExt Bext("Qty_MagnFieldExt",qty_map,dimension,FE_MAGNFIELDEXT);     qty_map.set_qty(&Bext);  

//consistency check
 if (bhom._dim !=  Bext._dim)     {std::cout << "main: inconsistency" << std::endl;abort();}
 if (bhom._FEord !=  Bext._FEord) {std::cout << "main: inconsistency" << std::endl;abort();}

 MagnFieldHomLagMult         bhom_lag_mult("Qty_MagnFieldHomLagMult",qty_map,1,FE_MAGNFIELDHOMLAGMULT);     qty_map.set_qty(&bhom_lag_mult);
 MagnFieldExtLagMult         Bext_lag_mult("Qty_MagnFieldExtLagMult",qty_map,1,FE_MAGNFIELDEXTLAGMULT);     qty_map.set_qty(&Bext_lag_mult);
 MagnFieldHomAdj                  bhom_adj("Qty_MagnFieldHomAdj",qty_map,dimension,FE_MAGNFIELDHOM);        qty_map.set_qty(&bhom_adj);
 MagnFieldHomLagMultAdj  bhom_lag_mult_adj("Qty_MagnFieldHomLagMultAdj",qty_map,1,FE_MAGNFIELDHOMLAGMULT);  qty_map.set_qty(&bhom_lag_mult_adj);

  Pressure  pressure("Qty_Pressure",qty_map,1,FE_PRESSURE);            qty_map.set_qty(&pressure);
  Velocity  velocity("Qty_Velocity",qty_map,dimension,FE_VELOCITY);   qty_map.set_qty(&velocity);  

  VelocityAdj  velocity_adj("Qty_VelocityAdj",qty_map,dimension,FE_VELOCITY);         qty_map.set_qty(&velocity_adj);  
  PressureAdj pressure_adj("Qty_PressureAdj",qty_map,1,FE_PRESSURE);                  qty_map.set_qty(&pressure_adj);
  DesVelocity des_velocity("Qty_DesVelocity",qty_map,dimension,FE_DESVELOCITY);       qty_map.set_qty(&des_velocity);
 
//consistency check
 if (velocity._dim !=  des_velocity._dim) {std::cout << "main: inconsistency" << std::endl;abort();}
 if (velocity._FEord !=  des_velocity._FEord) {std::cout << "main: inconsistency" << std::endl;abort();}

// #if TEMP_DEPS==1
  Temperature       temperature("Qty_Temperature",qty_map,1,FE_TEMPERATURE);      qty_map.set_qty(&temperature);  
  Density               density("Qty_Density",qty_map,1,0);                       qty_map.set_qty(&density);   
  Viscosity           viscosity("Qty_Viscosity",qty_map,1,0);                     qty_map.set_qty(&viscosity);
  HeatConductivity    heat_cond("Qty_HeatConductivity",qty_map,1,0);              qty_map.set_qty(&heat_cond);
  SpecificHeatP      spec_heatP("Qty_SpecificHeatP",qty_map,1,0);                 qty_map.set_qty(&spec_heatP);
// #endif  

  
//================================
//==== END Add QUANTITIES ========
//================================

  // ======== TimeLoop ===================================
  TimeLoop time_loop(utils); 
  time_loop._timemap.read();
  time_loop._timemap.print();
  
  // ====== EquationsMap =================================
  EquationsMap equations_map(utils,phys,qty_map,mesh,FEElements,qrule,time_loop);
  
//===============================================
//================== Add EQUATIONS  AND ======================
//========= associate an EQUATION to QUANTITIES ========
//========================================================

#if NS_EQUATIONS==1
std::vector<Quantity*> InternalVect_NS(2);
InternalVect_NS[QTYZERO] = &velocity;      velocity.SetPosInAssocEqn(0);
InternalVect_NS[QTYONE] = &pressure;       pressure.SetPosInAssocEqn(1);

  EqnNS* eqnNS=new EqnNS(InternalVect_NS,equations_map);
  equations_map.set_eqs(eqnNS);
 
 velocity.set_eqn(eqnNS);
 pressure.set_eqn(eqnNS);

#endif
  
#if NSAD_EQUATIONS==1
std::vector<Quantity*> InternalVect_NSAD(2);
InternalVect_NSAD[QTYZERO] = &velocity_adj;     velocity_adj.SetPosInAssocEqn(0);
InternalVect_NSAD[QTYONE]  = &pressure_adj;     pressure_adj.SetPosInAssocEqn(1);

  EqnNSAD* eqnNSAD=new EqnNSAD(InternalVect_NSAD,equations_map);
  equations_map.set_eqs(eqnNSAD);

  velocity_adj.set_eqn(eqnNSAD);
  pressure_adj.set_eqn(eqnNSAD);

#endif
  
#if MHD_EQUATIONS==1
std::vector<Quantity*> InternalVect_MHD(2);
InternalVect_MHD[QTYZERO] = &bhom;             bhom.SetPosInAssocEqn(0);
InternalVect_MHD[QTYONE]  = &bhom_lag_mult;    bhom_lag_mult.SetPosInAssocEqn(1);

  EqnMHD* eqnMHD=new EqnMHD(InternalVect_MHD,equations_map);
  equations_map.set_eqs(eqnMHD);

             bhom.set_eqn(eqnMHD);
    bhom_lag_mult.set_eqn(eqnMHD);
 
#endif

#if MHDAD_EQUATIONS==1
std::vector<Quantity*> InternalVect_MHDAD(2);
InternalVect_MHDAD[QTYZERO] = &bhom_adj;             bhom_adj.SetPosInAssocEqn(0);
InternalVect_MHDAD[QTYONE]  = &bhom_lag_mult_adj;    bhom_lag_mult_adj.SetPosInAssocEqn(1);
	
  EqnMHDAD* eqnMHDAD=new EqnMHDAD(InternalVect_MHDAD,equations_map);
  equations_map.set_eqs(eqnMHDAD);
  
           bhom_adj.set_eqn(eqnMHDAD);
  bhom_lag_mult_adj.set_eqn(eqnMHDAD);
  
#endif

#if MHDCONT_EQUATIONS==1
std::vector<Quantity*> InternalVect_MHDCONT(2);
InternalVect_MHDCONT[QTYZERO] = &Bext;            Bext.SetPosInAssocEqn(0);
InternalVect_MHDCONT[QTYONE]  = &Bext_lag_mult;   Bext_lag_mult.SetPosInAssocEqn(1);

  EqnMHDCONT* eqnMHDCONT = new EqnMHDCONT(InternalVect_MHDCONT,equations_map);
  equations_map.set_eqs(eqnMHDCONT);

                 Bext.set_eqn(eqnMHDCONT);
        Bext_lag_mult.set_eqn(eqnMHDCONT);
  
#endif  
   
//================================
//========= End add EQUATIONS  and ========
//========= associate an EQUATION to QUANTITIES ========
//================================

  equations_map.setDofBcOpIc();     //  /*TODO fileIO  for  Bc, init, and Ic*/
  equations_map.TransientSetup();  // reset the initial state (if restart) and print the Case   /*TODO fileIO */ 

//initialize specific data for specific equations
//all that happened previously was related to the standard data of EqnBase, basically  
#if MHDCONT_EQUATIONS==1
  eqnMHDCONT->init_equation_data();
#endif
  
//   equations_map.TransientLoop();   // perform the time evolution for all the equations  /*TODO fileIO*/

    optimization_loop(equations_map);  /////

// // //   eqnNS->FunctionIntegral (0,funzione);
// // //   eqnNS->FunctionIntegral (1,funzione);

//============= prepare default for next restart ==========  
// at this point, the run has been completed 
// well, we do not know whether for the whole time range time.N-M.xmf
// or the optimization loop stopped before
// we could also print the last step number
/*(iproc==0)*/  files.PrintRun(DEFAULT_LAST_RUN);  /*TODO fileIO*/

  // ============  log ================================
#ifdef FEMUS_HAVE_PETSC
  std::string petsc_femus_log = "petsc_main.log";
  std::ostringstream petsc_log;
  petsc_log <<  files.get_basepath() + "/" + files.get_frtmap().get("OUTPUT_DIR")
            << "/" << files.get_frtmap().get("OUTTIME_DIR") << "/" << petsc_femus_log;
  PetscViewer my_viewer;
  PetscViewerCreate(MPI_COMM_WORLD, &my_viewer);
  PetscViewerSetType(my_viewer, PETSCVIEWERASCII);
  PetscViewerFileSetName(my_viewer, petsc_log.str().c_str());   
  PetscLogView(my_viewer);
#endif

// ============  clean ================================
  equations_map.clean();  //deallocates the map of equations

  for (int fe=0; fe<QL; fe++) delete FEElements[fe];
  
  files.CloseCaseData();
// >>>>>>>>>>>>> END REDIRECT COUT
  std::cout.rdbuf(sbuf);  //it seems like you have to give the stream buffer
                          //back to cout !!!
                         // http://wordaligned.org/articles/cpp-streambufs
// >>>>>>>>>>>>> END REDIRECT COUT
  
  return 0;
}

//TODO change the QuantityMap constructor so that it receives the mesh directly


//manual breakpoint
// #ifdef HAVE_MPI
// MPI_Barrier(MPI_COMM_WORLD);
// #endif
//   std::cout << "***** Up to this point ****** " << std::endl; abort();
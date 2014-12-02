//C++ includes 
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <cstdlib>
#include <sstream>

// External library include ( LibMesh, PETSc...) ------------------------------
#include "FEMTTUConfig.h"

// FEMuS
#include "paral.hpp" 
#include "FemusInit.hpp"
#include "Files.hpp"
#include "Physics.hpp"
#include "GeomEl.hpp"
#include "MeshTwo.hpp"
#include "FETypeEnum.hpp"
#include "FEElemBase.hpp"
#include "QRule.hpp"
#include "EquationsMap.hpp"
#include "TimeLoop.hpp"
#include "Typedefs.hpp"
#include "Quantity.hpp"
#include "QTYnumEnum.hpp"
#include "Box.hpp"  //for the DOMAIN
#include "LinearSolverM.hpp"


// application 
#include "TempQuantities.hpp"
#include "TempPhysics.hpp"
#include "EqnT.hpp"

 
// =======================================
// TEMPERATURE + NS optimal control problem
// ======================================= 

 int main(int argc, char** argv) {

  // ====== FemusInit =====  //put this as the first call because mpi is initialized here
  FemusInit init(argc,argv);
  
 // ======= Files ========================
  Files files("./"); 
        files.ConfigureRestart();
        files.CheckIODirectories();
        files.RedirectCout();
        files.CopyInputFiles();
   // at this point everything is in the folder of the current run!!!!

  // ======= MyPhysics (implemented as child of Physics) ========================
  RunTimeMap<double> physics_map("Physics",files._output_path);
  TempPhysics phys(physics_map);
  const double Lref  =  phys._physrtmap.get("Lref");     // reference L

  // ======= Mesh =====
  RunTimeMap<double> mesh_map("Mesh",files._output_path);
  Mesh mesh(files,mesh_map,Lref);
  
  // ======= MyDomainShape  (optional, implemented as child of Domain) ====================
  RunTimeMap<double> box_map("Box",files._output_path);
  Box mybox(mesh.get_dim(),box_map);
      mybox.init(mesh.get_Lref());
  
  mesh.SetDomain(&mybox);    
  
  mesh.ReadMeshFile(); 
  mesh.PrintForVisualizationAllLEVAllVB();
  
  phys.set_mesh(&mesh);
  
  
// ======  QRule ================================ //so far we have ONLY ONE quadrature rule for all the equations
  QRule   qrule(mesh._GeomEl);

  // =======Abstract FEElems =====  //remember to delete the FE at the end
  std::vector<FEElemBase*> FEElements(QL);  //TODO what if we dont want to call the default constructor?!? AAA here no constructor is called!!! If you have a pointer the constructor is not called!
                                                     
  for (int fe=0; fe<QL; fe++) {
    FEElements[fe] = FEElemBase::build(mesh._GeomEl,fe);  //The order of the fe is established by the library
//sort of constructor
    FEElements[fe]->SetOrder(fe);
    FEElements[fe]->AssociateQRule(&qrule);
//end sort of constructor
    FEElements[fe]->evaluate_shape_at_qp();
  }

  // ======== TimeLoop ===================================
  TimeLoop time_loop(files); 
           time_loop._timemap.read();
           time_loop._timemap.print();

  // ===== QuantityMap =========================================
  QuantityMap  qty_map(phys);

//===============================================
//================== Add QUANTITIES ======================
//========================================================
  
  Temperature temperature("Qty_Temperature",qty_map,1,0/*biquadratic*/);     qty_map.set_qty(&temperature);
  Temperature temperature2("Qty_Temperature2",qty_map,1,1/*linear*/);        qty_map.set_qty(&temperature2);
  Temperature temperature3("Qty_Temperature3",qty_map,1,2/*constant*/);      qty_map.set_qty(&temperature3);
  // ===== end QuantityMap =========================================

  // ====== EquationsMap =================================
  EquationsMap equations_map(files,phys,qty_map,mesh,FEElements,qrule,time_loop);  //here everything is passed as BASE STUFF, like it should!
                                                                                   //the equations need: physical parameters, physical quantities, Domain, FE, QRule, Time discretization  
//===============================================
//================== Add EQUATIONS AND ======================
//========= associate an EQUATION to QUANTITIES ========
//========================================================
// not all the Quantities need an association with equation
//once you associate one quantity in the internal map of an equation, then it is immediately to be associated to that equation,
//   so this operation of set_eqn could be done right away in the moment when you put the quantity in the equation
 
// // // std::vector<Quantity*> InternalVect_Temp(3); 
std::vector<Quantity*> InternalVect_Temp(1); 

InternalVect_Temp[0] = &temperature;               temperature.SetPosInAssocEqn(0);
// InternalVect_Temp[1] = &temperature2;              temperature2.SetPosInAssocEqn(1);
// InternalVect_Temp[2] = &temperature3;              temperature3.SetPosInAssocEqn(2);

  EqnT* eqnT = new EqnT(InternalVect_Temp,equations_map);
  equations_map.set_eqs(eqnT);  
  
    for (uint l=0; l< mesh._NoLevels; l++)  eqnT->_solver[l]->set_solver_type(GMRES);
    eqnT->_Dir_pen_fl = 0;  //no penalty BC

        temperature.set_eqn(eqnT);
//         temperature2.set_eqn(eqnT);
//         temperature3.set_eqn(eqnT);

//================================ 
//========= End add EQUATIONS  and ========
//========= associate an EQUATION to QUANTITIES ========
//================================

//Ok now that the mesh file is there i want to COMPUTE the MG OPERATORS... but I want to compute them ONCE and FOR ALL,
//not for every equation... but the functions belong to the single equation... I need to make them EXTERNAL
// then I'll have A from the equation, PRL and REST from a MG object.
//So somehow i'll have to put these objects at a higher level... but so far let us see if we can COMPUTE and PRINT from HERE and not from the gencase
	 
  equations_map.setDofBcOpIc();    //once you have the list of the equations, you loop over them to initialize everything
  equations_map.TransientSetup();  // reset the initial state (if restart) and print the Case

  equations_map.TransientLoop();
//   phys.transient_loopPlusJ(equations_map);

// at this point, the run has been completed 
  files.PrintRunForRestart(DEFAULT_LAST_RUN);/*(iproc==0)*/  //============= prepare default for next restart ==========  
  files.log_petsc();
  
// ============  clean ================================
  // here we clean all that we allocated as new in the main
  equations_map.clean();  //deallocates the map of equations
  for (int fe=0; fe<QL; fe++)  {  delete FEElements[fe]; }

  mesh.clear();
  
  
  return 0;
}

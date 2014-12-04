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
#include "GaussPoints.hpp"
#include "EquationsMap.hpp"
#include "TimeLoop.hpp"
#include "Typedefs.hpp"
#include "Quantity.hpp"
#include "QTYnumEnum.hpp"
#include "Box.hpp"  //for the DOMAIN


// application 
#include "Temp_conf.hpp"
#include "TempQuantities.hpp"
#include "TempPhysics.hpp"
#include "EqnNS.hpp"
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
              phys.set_nondimgroups(); //implement it
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
  std::vector<Gauss>   qrule;
  qrule.reserve(VB);
  for (int vb=0;vb < VB; vb++) { 
          Gauss qrule_temp(mesh._GeomEl[vb]._geomel_id.c_str(),"fifth");  //TODO THEN WE CAN TRY WITH THE OTHER ORDERS AS WELL 
         qrule.push_back(qrule_temp);
  }
  
  // =======Abstract FEElems =====  //remember to delete the FE at the end
  std::vector< std::vector<FEElemBase*> >  FEElements_vec(VB); 
  std::vector<FEElemBase*> FEElements(QL);
  const std::string  FEFamily[QL] = {"biquadratic","linear","constant"}; 
  
  for (int fe=0; fe<QL; fe++) {
    FEElements[fe] = FEElemBase::build(mesh._GeomEl[VV],fe);  //The order of the fe is established by the library
    FEElements[fe]->AssociateQRule(qrule);
    FEElements[fe]->evaluate_shape_at_qp(fe);
  }

  for (int vb=0;vb < VB; vb++) { 
  FEElements_vec[vb] = FEElements;
  }
  
  // ======== TimeLoop ===================================
  TimeLoop time_loop(files); 
           time_loop._timemap.read();
           time_loop._timemap.print();

  // ===== QuantityMap =========================================
  QuantityMap  qty_map(phys);

  Temperature temperature("Qty_Temperature",qty_map,1,FE_TEMPERATURE);     qty_map.set_qty(&temperature);
  TempLift       templift("Qty_TempLift",qty_map,1,FE_TEMPERATURE);        qty_map.set_qty(&templift);  
  TempAdj         tempadj("Qty_TempAdj",qty_map,1,FE_TEMPERATURE);         qty_map.set_qty(&tempadj);  
  TempDes         tempdes("Qty_TempDes",qty_map,1,FE_TEMPERATURE);         qty_map.set_qty(&tempdes);  
  Pressure       pressure("Qty_Pressure",qty_map,1,FE_PRESSURE);           qty_map.set_qty(&pressure);
  Velocity       velocity("Qty_Velocity",qty_map,mesh.get_dim(),FE_VELOCITY);   qty_map.set_qty(&velocity);  

#if FOURTH_ROW==1
  Pressure pressure_2("Qty_Pressure_2",qty_map,1,T4_ORD);            qty_map.set_qty(&pressure_2);
#endif 
  
  // ===== end QuantityMap =========================================

  // ====== EquationsMap =================================
  EquationsMap equations_map(files,phys,qty_map,mesh,FEElements_vec,qrule,time_loop);  //here everything is passed as BASE STUFF, like it should!
                                                                                   //the equations need: physical parameters, physical quantities, Domain, FE, QRule, Time discretization  
  
//===============================================
//================== Add EQUATIONS AND ======================
//========= associate an EQUATION to QUANTITIES ========
//========================================================
// not all the Quantities need an association with equation
//once you associate one quantity in the internal map of an equation, then it is immediately to be associated to that equation,
//   so this operation of set_eqn could be done right away in the moment when you put the quantity in the equation
 
#if NS_EQUATIONS==1
  //to retrieve a quantity i can take it from the qtymap of the problem  //but here, in the main, i can take that quantity directly...
// std::vector<Quantity*> InternalVect_NS;
// InternalVect_NS.push_back(&velocity);      
// InternalVect_NS.push_back(&pressure);        
std::vector<Quantity*> InternalVect_NS(2); 
InternalVect_NS[0] = &velocity;  velocity.SetPosInAssocEqn(0);
InternalVect_NS[1] = &pressure;  pressure.SetPosInAssocEqn(1);

  EqnNS* eqnNS = new EqnNS(InternalVect_NS,equations_map);  equations_map.set_eqs(eqnNS);
  
           velocity.set_eqn(eqnNS);
           pressure.set_eqn(eqnNS);
#endif
  
#if T_EQUATIONS==1
std::vector<Quantity*> InternalVect_Temp( 3 + FOURTH_ROW );  //of course this must be exactly equal to the following number of get_qty
                                                             // can I do this dynamic? 
                                                             // well, the following order is essential, because it is the same order 
                                                             // as the BLOCKS in the MATRIX, but at least with add you can avoid setting also the SIZE explicitly.
                                                             //The order in which you put the push_back instructions is essential and it gives you 
                                                             //the order in the std::vector!!!
InternalVect_Temp[0] = &temperature;       temperature.SetPosInAssocEqn(0);
InternalVect_Temp[1] = &templift;             templift.SetPosInAssocEqn(1);
InternalVect_Temp[2] = &tempadj;               tempadj.SetPosInAssocEqn(2);

#if FOURTH_ROW==1
InternalVect_Temp[3] = &pressure_2;         pressure_2.SetPosInAssocEqn(3);
#endif

  EqnT* eqnT = new EqnT(InternalVect_Temp,equations_map);
  equations_map.set_eqs(eqnT);  

        temperature.set_eqn(eqnT);
           templift.set_eqn(eqnT);
            tempadj.set_eqn(eqnT);
   #if FOURTH_ROW==1
         pressure_2.set_eqn(eqnT);
   #endif

#endif   

//================================ 
//========= End add EQUATIONS  and ========
//========= associate an EQUATION to QUANTITIES ========
//================================

//Ok now that the mesh file is there i want to COMPUTE the MG OPERATORS... but I want to compute them ONCE and FOR ALL,
//not for every equation... but the functions belong to the single equation... I need to make them EXTERNAL
// then I'll have A from the equation, PRL and REST from a MG object.
//So somehow i'll have to put these objects at a higher level... but so far let us see if we can COMPUTE and PRINT from HERE and not from the gencase
	 
//   eqnNS->ComputeMatrix();  //CLEARLY THIS FUNCTION DOES NOT WORK AT THIS POINT, because not all the data in the mesh class are filled here! 
                           //In fact, part of them is only filled by gencase
	 
  equations_map.setDofBcOpIc();    //once you have the list of the equations, you loop over them to initialize everything
  equations_map.TransientSetup();  // reset the initial state (if restart) and print the Case

  phys.transient_loopPlusJ(equations_map);

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

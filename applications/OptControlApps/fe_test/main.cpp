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
#include "GeomEl.hpp"
#include "MeshTwo.hpp"
#include "GenCase.hpp"
#include "FETypeEnum.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "TimeLoop.hpp"
#include "Typedefs.hpp"
#include "Quantity.hpp"
#include "QTYnumEnum.hpp"
#include "Box.hpp"  //for the DOMAIN
#include "LinearEquationSolver.hpp"


// application 
#include "TempQuantities.hpp"
#include "EqnT.hpp"

#ifdef HAVE_LIBMESH
#include "libmesh/libmesh.h"
#endif 

// =======================================
// Test for finite element families
// ======================================= 

 int main(int argc, char** argv) {

#ifdef HAVE_LIBMESH
   libMesh::LibMeshInit init(argc,argv);
#else   
   FemusInit init(argc,argv);
#endif
  
 // ======= Files ========================
  Files files; 
        files.ConfigureRestart();
        files.CheckIODirectories();
        files.CopyInputFiles();   // at this point everything is in the folder of the current run!!!!
        files.RedirectCout();

  // ======= MyPhysics (implemented as child of Physics) ========================
  FemusInputParser<double> physics_map("Physics",files._output_path);
  const double Lref  =  physics_map.get("Lref");     // reference L

  // ======= Mesh =====
  FemusInputParser<double> mesh_map("Mesh",files._output_path);

  GenCase mesh(files,mesh_map,"");
          mesh.SetLref(1.);  
	  
  // ======= MyDomainShape  (optional, implemented as child of Domain) ====================
  FemusInputParser<double> box_map("Box",files._output_path);
  Box mybox(mesh.get_dim(),box_map);
      mybox.InitAndNondimensionalize(mesh.get_Lref());

          mesh.SetDomain(&mybox);    
	  
          mesh.GenerateCase();

          mesh.SetLref(Lref);
      mybox.InitAndNondimensionalize(mesh.get_Lref());
	  
          mesh.ReadMeshFileAndNondimensionalize(); 
	  mesh.PrintMultimeshXdmf();
          mesh.PrintForVisualizationAllLEVAllVB();

  // ======  QRule ================================
  std::vector<Gauss>   qrule;
  qrule.reserve(mesh.get_dim());
  for (int idim=0;idim < mesh.get_dim(); idim++) {
         Gauss qrule_temp(mesh.GetGeomEl(idim,mesh._mesh_order)._geomel_id.c_str(),"fifth");
         qrule.push_back(qrule_temp);
  }
  
  
  // =======Abstract FEElems =====  //remember to delete the FE at the end
  const std::string  FEFamily[QL] = {"biquadratic","linear","constant"}; 
  std::vector< std::vector<elem_type*> > FEElemType_vec(mesh.get_dim());
  for (int idim=0;idim < mesh.get_dim(); idim++)    FEElemType_vec[idim].resize(QL);

  for (int idim=0;idim < mesh.get_dim(); idim++) {
    for (int fe=0; fe<QL; fe++) {
       FEElemType_vec[idim][fe] = elem_type::build(mesh.GetGeomEl(idim,mesh._mesh_order)._geomel_id.c_str(),fe, qrule[idim].GetGaussOrderString().c_str());
       FEElemType_vec[idim][fe]->EvaluateShapeAtQP(mesh.GetGeomEl(idim,mesh._mesh_order)._geomel_id.c_str(),fe);
      }
    }
                                                     
  // ======== TimeLoop ===================================
  TimeLoop time_loop(files); 
           time_loop._timemap.read();
           time_loop._timemap.print();

  // ===== QuantityMap =========================================
  QuantityMap  qty_map(mesh,&physics_map);

//===============================================
//================== Add QUANTITIES ======================
//========================================================
  
  Temperature temperature("Qty_Temperature",qty_map,1,0/*biquadratic*/);     qty_map.set_qty(&temperature);
  Temperature temperature2("Qty_Temperature2",qty_map,1,1/*linear*/);        qty_map.set_qty(&temperature2);
  Temperature temperature3("Qty_Temperature3",qty_map,1,2/*constant*/);      qty_map.set_qty(&temperature3);
  // ===== end QuantityMap =========================================

  // ====== MultiLevelProblemTwo =================================
  MultiLevelProblemTwo equations_map(files,physics_map,qty_map,mesh,FEElemType_vec,qrule);  //here everything is passed as BASE STUFF, like it should!
                                                                                   //the equations need: physical parameters, physical quantities, Domain, FE, QRule, Time discretization  
//===============================================
//================== Add EQUATIONS AND ======================
//========= associate an EQUATION to QUANTITIES ========
//========================================================
// not all the Quantities need an association with equation
//once you associate one quantity in the internal map of an equation, then it is immediately to be associated to that equation,
//   so this operation of set_eqn could be done right away in the moment when you put the quantity in the equation
 
std::vector<Quantity*> InternalVect_Temp(3); 
// std::vector<Quantity*> InternalVect_Temp(1); 

InternalVect_Temp[0] = &temperature;               temperature.SetPosInAssocEqn(0);
InternalVect_Temp[1] = &temperature2;              temperature2.SetPosInAssocEqn(1);
InternalVect_Temp[2] = &temperature3;              temperature3.SetPosInAssocEqn(2);

  EqnT* eqnT = new EqnT(InternalVect_Temp,equations_map);
  equations_map.set_eqs(eqnT);  
  
    for (uint l=0; l< mesh._NoLevels; l++)  eqnT->_solver[l]->set_solver_type(GMRES);
    eqnT->_Dir_pen_fl = 0;  //no penalty BC

        temperature.set_eqn(eqnT);
        temperature2.set_eqn(eqnT);
        temperature3.set_eqn(eqnT);

//================================ 
//========= End add EQUATIONS  and ========
//========= associate an EQUATION to QUANTITIES ========
//================================

//Ok now that the mesh file is there i want to COMPUTE the MG OPERATORS... but I want to compute them ONCE and FOR ALL,
//not for every equation... but the functions belong to the single equation... I need to make them EXTERNAL
// then I'll have A from the equation, PRL and REST from a MG object.
//So somehow i'll have to put these objects at a higher level... but so far let us see if we can COMPUTE and PRINT from HERE and not from the gencase
	 
  equations_map.setDofBcOpIc();    //once you have the list of the equations, you loop over them to initialize everything
  time_loop.TransientSetup(equations_map);  // reset the initial state (if restart) and print the Case

  time_loop.TransientLoop(equations_map);

// at this point, the run has been completed 
  files.PrintRunForRestart(DEFAULT_LAST_RUN);/*(iproc==0)*/  //============= prepare default for next restart ==========  
  files.log_petsc();
  
// ============  clean ================================
  // here we clean all that we allocated as new in the main
  equations_map.clean();  //deallocates the map of equations
  mesh.clear();
  
  
  return 0;
}

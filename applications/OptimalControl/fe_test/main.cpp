//C++ includes 
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <cstdlib>
#include <sstream>

// External library include ( LibMesh, PETSc...) ------------------------------
#include "FemusConfig.hpp"

// FEMuS
#include "paral.hpp" 
#include "FemusInit.hpp"
#include "Files.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "GenCase.hpp"
#include "FETypeEnum.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "TimeLoop.hpp"
#include "Typedefs.hpp"
#include "Quantity.hpp"
#include "Box.hpp"  //for the DOMAIN
#include "LinearEquationSolver.hpp"
#include "XDMFWriter.hpp"

// application 
#include "TempQuantities.hpp"

#ifdef HAVE_LIBMESH
#include "libmesh/libmesh.h"
#endif 

 void GenMatRhsT(MultiLevelProblem &ml_prob);


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
  FemusInputParser<double> physics_map("Physics",files.GetOutputPath());
  const double Lref  =  physics_map.get("Lref");     // reference L

  // ======= Mesh =====
  const unsigned NoLevels = 3;
  const unsigned dim = 2;
  const GeomElType geomel_type = QUAD;

  GenCase mesh(NoLevels,dim,geomel_type,"");
          mesh.SetLref(1.);  
	  
  // ======= MyDomainShape  (optional, implemented as child of Domain) ====================
  FemusInputParser<double> box_map("Box",files.GetOutputPath());

  Box mybox(mesh.get_dim(),box_map);
      mybox.InitAndNondimensionalize(mesh.get_Lref());

          mesh.SetDomain(&mybox);    
	  
          mesh.GenerateCase(files.GetOutputPath());

          mesh.SetLref(Lref);
      mybox.InitAndNondimensionalize(mesh.get_Lref());
	  
          XDMFWriter::ReadMeshAndNondimensionalizeBiquadraticHDF5(files.GetOutputPath(),mesh); 
	  XDMFWriter::PrintMeshXDMF(files.GetOutputPath(),mesh,BIQUADR_FE);
          XDMFWriter::PrintMeshLinear(files.GetOutputPath(),mesh);

	  
  // ===== QuantityMap =========================================
  QuantityMap  qty_map;
  qty_map.SetMeshTwo(&mesh);
  qty_map.SetInputParser(&physics_map);

//===============================================
//================== Add QUANTITIES ======================
//========================================================
  
  Temperature temperature("Qty_Temperature",qty_map,1,0/*biquadratic*/);     qty_map.AddQuantity(&temperature);
//   Temperature temperature2("Qty_Temperature2",qty_map,1,1/*linear*/);        qty_map.AddQuantity(&temperature2);
//   Temperature temperature3("Qty_Temperature3",qty_map,1,2/*constant*/);      qty_map.AddQuantity(&temperature3);
  // ===== end QuantityMap =========================================

  // ====== Start new main =================================
  MultiLevelMesh ml_msh;
  ml_msh.GenerateCoarseBoxMesh(8,8,0,0,1,0,1,0,1,QUAD9,"fifth"); //   ml_msh.GenerateCoarseBoxMesh(numelemx,numelemy,numelemz,xa,xb,ya,yb,za,zb,elemtype,"seventh");
  ml_msh.RefineMesh(NoLevels,NoLevels,NULL);
  ml_msh.PrintInfo();

  ml_msh.SetDomain(&mybox);    
  
  MultiLevelSolution ml_sol(&ml_msh);
  ml_sol.AddSolution("FAKE",LAGRANGE,SECOND,0);

  MultiLevelProblem ml_prob(&ml_sol);  
  ml_prob.SetMeshTwo(&mesh);
  ml_prob.SetQuadratureRuleAllGeomElems("fifth");
//   ml_prob.SetElemTypeAllDims();
  ml_prob.SetInputParser(&physics_map); 
  ml_prob.SetQtyMap(&qty_map);

//===============================================
//================== Add EQUATIONS AND ======================
//========= associate an EQUATION to QUANTITIES ========
//========================================================

  SystemTwo &  eqnT = ml_prob.add_system<SystemTwo>("Eqn_T");
          eqnT.AddSolutionToSystemPDE("FAKE");
          eqnT.AddUnknownToSystemPDE(&temperature); 
//           eqnT.AddUnknownToSystemPDE(&temperature2); 
//           eqnT.AddUnknownToSystemPDE(&temperature3); 
          eqnT.SetAssembleFunction(GenMatRhsT);

//================================ 
//========= End add EQUATIONS  and ========
//========= associate an EQUATION to QUANTITIES ========
//================================

//Ok now that the mesh file is there i want to COMPUTE the MG OPERATORS... but I want to compute them ONCE and FOR ALL,
//not for every equation... but the functions belong to the single equation... I need to make them EXTERNAL
// then I'll have A from the equation, PRL and REST from a MG object.
//So somehow i'll have to put these objects at a higher level... but so far let us see if we can COMPUTE and PRINT from HERE and not from the gencase
	 
   for (MultiLevelProblem::const_system_iterator eqn = ml_prob.begin(); eqn != ml_prob.end(); eqn++) {
     
        SystemTwo* sys = static_cast<SystemTwo*>(eqn->second);
// //=====================
    sys -> init_two();
    sys -> _LinSolver[0]->set_solver_type(GMRES);  //if I keep PREONLY it doesn't run

//=====================
    sys -> init_unknown_vars();
//=====================
    sys -> _dofmap.ComputeMeshToDof();
//=====================
    sys -> initVectors();
//=====================
    sys -> Initialize();
//=====================
    sys -> _bcond.GenerateBdc();
//=====================
    GenCase::ReadMGOps(files.GetOutputPath(),sys);
    
    }
    
  // ======== Loop ===================================
  FemusInputParser<double> loop_map("TimeLoop",files.GetOutputPath());
  TimeLoop time_loop(files,loop_map); 

  time_loop.TransientSetup(ml_prob);  // reset the initial state (if restart) and print the Case

  time_loop.TransientLoop(ml_prob);

// at this point, the run has been completed 
  files.PrintRunForRestart(DEFAULT_LAST_RUN);
  files.log_petsc();
  
// ============  clean ================================
  ml_prob.clear();
  mesh.clear();
  
  
  return 0;
}

//C++ includes 
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <cstdlib>
#include <sstream>

// External library include ( LibMesh, PETSc...)
#include "FEMTTUConfig.h"

// FEMuS
#include "paral.hpp" 
#include "FemusInit.hpp"
#include "Files.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "GenCase.hpp"
#include "FETypeEnum.hpp"
#include "GaussPoints.hpp"
#include "MultiLevelProblem.hpp"
#include "ElemType.hpp"
#include "TimeLoop.hpp"
#include "Typedefs.hpp"
#include "Quantity.hpp"
#include "QTYnumEnum.hpp"
#include "Box.hpp"  //for the DOMAIN
#include "XDMFWriter.hpp"

// application 
#include "TempQuantities.hpp"
#include "OptLoop.hpp"


#ifdef HAVE_LIBMESH
#include "libmesh/libmesh.h"
#endif


void  GenMatRhsT(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assemble_matrix);
void  GenMatRhsNS(MultiLevelProblem &ml_prob, unsigned Level, const unsigned &gridn, const bool &assemble_matrix);

// =======================================
// TEMPERATURE + NS optimal control problem
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
        files.CopyInputFiles();
        files.RedirectCout();

  // ======= Physics Input Parser ========================
  FemusInputParser<double> physics_map("Physics",files.GetOutputPath());

  const double rhof   = physics_map.get("rho0");
  const double Uref   = physics_map.get("Uref");
  const double Lref   = physics_map.get("Lref");
  const double  muf   = physics_map.get("mu0");

  const double  _pref = rhof*Uref*Uref;           physics_map.set("pref",_pref);
  const double   _Re  = (rhof*Uref*Lref)/muf;     physics_map.set("Re",_Re);
  const double   _Fr  = (Uref*Uref)/(9.81*Lref);  physics_map.set("Fr",_Fr);
  const double   _Pr  = muf/rhof;                 physics_map.set("Pr",_Pr);

  // ======= Mesh =====
  FemusInputParser<double> mesh_map("Mesh",files.GetOutputPath());
  GenCase mesh(mesh_map,"inclQ2D2x2.gam");
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
	  
  //gencase is dimensionalized, meshtwo is nondimensionalized
  //since the meshtwo is nondimensionalized, all the BC and IC are gonna be implemented on a nondimensionalized mesh
  //now, a mesh may or may not have an associated domain
  //moreover, a mesh may or may not be read from file
  //the generation is dimensional, the nondimensionalization occurs later
  //Both the Mesh and the optional domain must be nondimensionalized
  //first, we have to say if the mesh has a shape or not
  //that depends on the application, it must be put at the main level
  //then, after you know the shape, you may or may not generate the mesh with that shape 
  //the two things are totally independent, and related to the application, not to the library

  // ======== Loop ===================================
  FemusInputParser<double> loop_map("TimeLoop",files.GetOutputPath());
  OptLoop opt_loop(files, loop_map); 
   
  // ===== QuantityMap : this is like the MultilevelSolution =========================================
  QuantityMap  qty_map;
  qty_map.SetMeshTwo(&mesh);
  qty_map.SetInputParser(&physics_map);

  Temperature temperature("Qty_Temperature",qty_map,1,QQ);          qty_map.AddQuantity(&temperature);
  TempLift       templift("Qty_TempLift",qty_map,1,QQ,opt_loop);    qty_map.AddQuantity(&templift);  
  TempAdj         tempadj("Qty_TempAdj",qty_map,1,QQ);              qty_map.AddQuantity(&tempadj);  
  TempDes         tempdes("Qty_TempDes",qty_map,1,QQ);              qty_map.AddQuantity(&tempdes);  //this is not going to be an Unknown!
  Pressure       pressure("Qty_Pressure",qty_map,1,LL);                qty_map.AddQuantity(&pressure);
  Velocity       velocity("Qty_Velocity",qty_map,mesh.get_dim(),QQ);   qty_map.AddQuantity(&velocity);  

#if FOURTH_ROW==1
  Pressure2 pressure_2("Qty_Pressure_2",qty_map,1,KK);            qty_map.AddQuantity(&pressure_2);
#endif 
  
  // ===== end QuantityMap =========================================
  
  // ====== Start new main =================================
  MultiLevelMesh ml_msh;
  ml_msh.GenerateCoarseBoxMesh(8,8,0,0,1,0,1,0,0,QUAD9,"fifth"); //   ml_msh.GenerateCoarseBoxMesh(numelemx,numelemy,numelemz,xa,xb,ya,yb,za,zb,elemtype,"seventh");
  ml_msh.RefineMesh(mesh_map.get("nolevels"),mesh_map.get("nolevels"),NULL);
  ml_msh.PrintInfo();
  
  MultiLevelSolution ml_sol(&ml_msh);
  ml_sol.AddSolution("FAKE",LAGRANGE,SECOND,0);
  
  MultiLevelProblem ml_prob(&ml_sol);
  ml_prob.SetMeshTwo(&mesh);
  ml_prob.SetQruleAndElemType("fifth");
  ml_prob.SetInputParser(&physics_map);
  ml_prob.SetQtyMap(&qty_map); 

//===============================================
//================== Add EQUATIONS AND ======================
//========= associate an EQUATION to QUANTITIES ========
//========================================================
// not all the Quantities need to be unknowns of an equation

  SystemTwo & eqnNS = ml_prob.add_system<SystemTwo>("Eqn_NS");
          eqnNS.AddSolutionToSystemPDE("FAKE");
          eqnNS.AddUnknownToSystemPDE(&velocity); 
          eqnNS.AddUnknownToSystemPDE(&pressure);
	  eqnNS.SetAssembleFunction(GenMatRhsNS);
  
  SystemTwo & eqnT = ml_prob.add_system<SystemTwo>("Eqn_T");
         eqnT.AddSolutionToSystemPDE("FAKE");
         eqnT.AddUnknownToSystemPDE(&temperature);
         eqnT.AddUnknownToSystemPDE(&templift);
         eqnT.AddUnknownToSystemPDE(&tempadj);
#if FOURTH_ROW==1
         eqnT.AddUnknownToSystemPDE(&pressure_2);   //the order in which you add defines the order in the matrix as well, so it is in tune with the assemble function
#endif
	 eqnT.SetAssembleFunction(GenMatRhsT);
  
//================================ 
//========= End add EQUATIONS  and ========
//========= associate an EQUATION to QUANTITIES ========
//================================

//Ok now that the mesh file is there i want to COMPUTE the MG OPERATORS... but I want to compute them ONCE and FOR ALL,
//not for every equation... but the functions belong to the single equation... I need to make them EXTERNAL
// then I'll have A from the equation, PRL and REST from a MG object.
//So somehow i'll have to put these objects at a higher level... but so far let us see if we can COMPUTE and PRINT from HERE and not from the gencase
	 
 
//once you have the list of the equations, you loop over them to initialize everything

   for (MultiLevelProblem::const_system_iterator eqn = ml_prob.begin(); eqn != ml_prob.end(); eqn++) {
     
   SystemTwo* sys = static_cast<SystemTwo*>(eqn->second);
//=====================
    sys -> init();     //the dof map is built here based on all the solutions associated with that system
    sys -> _LinSolver[0]->set_solver_type(GMRES);  //if I keep PREONLY it doesn't run

//=====================
    sys -> init_sys();
//=====================
    sys -> _dofmap.ComputeMeshToDof();
//=====================
    sys -> initVectors();
//=====================
    sys -> Initialize();         //why do they do this BEFORE the dofmap?
///=====================
    sys -> _bcond.GenerateBdc(); //why do they do this BEFORE the dofmap?
//=====================
    sys -> ReadMGOps(files.GetOutputPath());
    
    } 
	 
	 
  opt_loop.TransientSetup(ml_prob);  // reset the initial state (if restart) and print the Case

  opt_loop.optimization_loop(ml_prob);

// at this point, the run has been completed 
  files.PrintRunForRestart(DEFAULT_LAST_RUN);
  files.log_petsc();
  
// ============  clean ================================
  ml_prob.clear();
  mesh.clear();
  
  
  return 0;
}
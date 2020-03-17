//C++ includes
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <cstdlib>
#include <sstream>

#include "FemusConfig.hpp"
#include "paral.hpp"
#include "FemusInit.hpp"
#include "Files.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "GenCase.hpp"
#include "FETypeEnum.hpp"
#include "GaussPoints.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "ElemType.hpp"
#include "TimeLoop.hpp"
#include "Typedefs.hpp"
#include "Quantity.hpp"
#include "Box.hpp"  //for the DOMAIN
#include "XDMFWriter.hpp"

// application
#include "TempQuantities.hpp"
#include "OptLoop.hpp"


#ifdef HAVE_LIBMESH
#include "libmesh/libmesh.h"
#endif

void  GenMatRhsT(MultiLevelProblem &ml_prob);
void  GenMatRhsNS(MultiLevelProblem &ml_prob);

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
  const unsigned NoLevels = 3;
  const unsigned dim = 2;
  const GeomElType geomel_type = QUAD;
  GenCase mesh(NoLevels,dim,geomel_type,"inclQ2D2x2.gam");
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

  // ===== QuantityMap : this is like the MultilevelSolution =========================================
  QuantityMap  qty_map;
  qty_map.SetMeshTwo(&mesh);
  qty_map.SetInputParser(&physics_map);

  Pressure       pressure("Qty_Pressure",qty_map,1,LL);             qty_map.AddQuantity(&pressure);
  VelocityX     velocityX("Qty_Velocity0",qty_map,1,QQ);            qty_map.AddQuantity(&velocityX);
  VelocityY     velocityY("Qty_Velocity1",qty_map,1,QQ);            qty_map.AddQuantity(&velocityY);
  Temperature temperature("Qty_Temperature",qty_map,1,QQ);          qty_map.AddQuantity(&temperature);
  TempLift       templift("Qty_TempLift",qty_map,1,QQ);             qty_map.AddQuantity(&templift);
  TempAdj         tempadj("Qty_TempAdj",qty_map,1,QQ);              qty_map.AddQuantity(&tempadj);
  // ===== end QuantityMap =========================================

  // ====== Start new main =================================

  MultiLevelMesh ml_msh;
  ml_msh.GenerateCoarseBoxMesh(8,8,0,0,1,0,2,0,0,QUAD9,"fifth"); //   ml_msh.GenerateCoarseBoxMesh(numelemx,numelemy,numelemz,xa,xb,ya,yb,za,zb,elemtype,"fifth");
  ml_msh.RefineMesh(NoLevels,NoLevels,NULL);
  ml_msh.PrintInfo();

  ml_msh.SetWriter(XDMF);
  //ml_msh.GetWriter()->write(files.GetOutputPath(),"biquadratic");

  ml_msh.SetDomain(&mybox);

  MultiLevelSolution ml_sol(&ml_msh);
  ml_sol.AddSolution("Qty_Temperature",LAGRANGE,SECOND,0);
  ml_sol.AddSolution("Qty_TempLift",LAGRANGE,SECOND,0);
  ml_sol.AddSolution("Qty_TempAdj",LAGRANGE,SECOND,0);
  ml_sol.AddSolutionVector(ml_msh.GetDimension(),"Qty_Velocity",LAGRANGE,SECOND,0);
  ml_sol.AddSolution("Qty_Pressure",LAGRANGE,FIRST,0);
  ml_sol.AddSolution("Qty_TempDes",LAGRANGE,SECOND,0,false); //this is not going to be an Unknown! //moreover, this is not going to need any BC (i think they are excluded with "false") // I would like to have a Solution that is NOT EVEN related to the mesh at all... just like a function "on-the-fly"

  // ******* Set problem *******
  MultiLevelProblem ml_prob(&ml_sol);
  ml_prob.SetMeshTwo(&mesh);
  ml_prob.SetQuadratureRuleAllGeomElems("fifth");
//   ml_prob.SetElemTypeAllDims();
  ml_prob.SetInputParser(&physics_map);
  ml_prob.SetQtyMap(&qty_map);

  // ******* Initial condition *******
  ml_sol.Initialize("Qty_Temperature",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_TempLift",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_TempAdj",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_Velocity0",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_Velocity1",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_Pressure",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_TempDes",SetInitialCondition, &ml_prob);

  /// @todo you have to call this before you can print
  /// @todo I can also call it after instantiation MLProblem
  /// @todo I cannot call it with "All" and with a FUNCTION, because I need the string for "All" as a variable
  /// @todo Have to say that you have to call this initialize BEFORE the generation of the boundary conditions;
  /// if you called this after, it would superimpose the BOUNDARY VALUES
  /// @todo you have to initialize also those variables which are NOT unknowns!

  // ******* Set boundary function function *******
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Generate boundary conditions *******
  ml_sol.GenerateBdc("Qty_Temperature","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_TempLift","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_TempAdj","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_Velocity0","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_Velocity1","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_Pressure","Steady",&ml_prob);


  // ******* Debug *******
  ml_sol.SetWriter(VTK);
  std::vector<std::string> print_vars(1); print_vars[0] = "All"; // we should find a way to make this easier
  ml_sol.GetWriter()->Write(files.GetOutputPath(),"biquadratic",print_vars);


//===============================================
//================== Add EQUATIONS AND ======================
//========= associate an EQUATION to QUANTITIES ========
//========================================================
// not all the Quantities need to be unknowns of an equation

  SystemTwo & eqnNS = ml_prob.add_system<SystemTwo>("Eqn_NS");

          eqnNS.AddSolutionToSystemPDEVector(ml_msh.GetDimension(),"Qty_Velocity");
          eqnNS.AddSolutionToSystemPDE("Qty_Pressure");

          eqnNS.AddUnknownToSystemPDE(&velocityX);
          eqnNS.AddUnknownToSystemPDE(&velocityY);
          eqnNS.AddUnknownToSystemPDE(&pressure);

	  eqnNS.SetAssembleFunction(GenMatRhsNS);


  SystemTwo & eqnT = ml_prob.add_system<SystemTwo>("Eqn_T");

         eqnT.AddSolutionToSystemPDE("Qty_Temperature");
         eqnT.AddSolutionToSystemPDE("Qty_TempLift");
         eqnT.AddSolutionToSystemPDE("Qty_TempAdj");

         eqnT.AddUnknownToSystemPDE(&temperature);
         eqnT.AddUnknownToSystemPDE(&templift);
         eqnT.AddUnknownToSystemPDE(&tempadj); //the order in which you add defines the order in the matrix as well, so it is in tune with the assemble function

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

  // ******* set MG-Solver *******
  sys->SetMgType(F_CYCLE);
  sys->SetAbsoluteLinearConvergenceTolerance(1.e-10);
  sys->SetNonLinearConvergenceTolerance(1.e-10);//1.e-5
  sys->SetNumberPreSmoothingStep(1);
  sys->SetNumberPostSmoothingStep(1);
  sys->SetMaxNumberOfLinearIterations(8);     //2
  sys->SetMaxNumberOfNonLinearIterations(15); //10

  // ******* Set Preconditioner *******
  sys->SetLinearEquationSolverType(FEMuS_DEFAULT);//FEMuS_ASM,VANKA_SMOOTHER

  // ******* init *******
  sys->init();

  // ******* Set Smoother *******
  sys->SetSolverFineGrids(GMRES);
  sys->SetPreconditionerFineGrids(ILU_PRECOND);
  sys->SetTolerances(1.e-12,1.e-20,1.e+50,20);   /// what the heck do these parameters mean?

    // ******* Add variables to be solved *******  /// do we always need this?
  sys->ClearVariablesToBeSolved();
  sys->AddVariableToBeSolved("All");


  // ******* For Gmres Preconditioner only *******
  sys->SetDirichletBCsHandling(ELIMINATION);

  // ******* For Gmres Preconditioner only *******
//   sys->solve();

//=====================
    sys -> init_two();     //the dof map is built here based on all the solutions associated with that system
    sys -> _LinSolver[0]->set_solver_type(GMRES);  //if I keep PREONLY it doesn't run

//=====================
    sys -> init_unknown_vars();
//=====================
    sys -> _dofmap.ComputeMeshToDof();
//=====================
    sys -> initVectors();
//=====================
    sys -> Initialize();
///=====================
    sys -> _bcond.GenerateBdc();
//=====================
    GenCase::ReadMGOps(files.GetOutputPath(),sys);

    }

  // ======== Loop ===================================
  FemusInputParser<double> loop_map("TimeLoop",files.GetOutputPath());
  OptLoop opt_loop(files, loop_map);

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

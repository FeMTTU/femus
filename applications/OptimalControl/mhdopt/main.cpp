//C++ includes
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <cstdlib>
#include <sstream>

// FEMuS
#include "FemusConfig.hpp"
#include "paral.hpp"
#include "FemusDefault.hpp"
#include "FemusInit.hpp"
#include "Files.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"
#include "GenCase.hpp"
#include "FETypeEnum.hpp"
#include "ElemType.hpp"
#include "TimeLoop.hpp"
#include "Typedefs.hpp"
#include "CmdLine.hpp"
#include "Quantity.hpp"
#include "Box.hpp"  //for the DOMAIN
#include "XDMFWriter.hpp"


// application includes
#include "OptLoop.hpp"
#include "OptQuantities.hpp"

#ifdef HAVE_LIBMESH
#include "libmesh/libmesh.h"
#endif

using namespace femus;

  void GenMatRhsNS(MultiLevelProblem &ml_prob);
  void GenMatRhsNSAD(MultiLevelProblem &ml_prob);
  void GenMatRhsMHD(MultiLevelProblem &ml_prob);
  void GenMatRhsMHDAD(MultiLevelProblem &ml_prob);
  void GenMatRhsMHDCONT(MultiLevelProblem &ml_prob);

// =======================================
// MHD optimal control problem
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

  // ======= Physics ========================
  FemusInputParser<double> physics_map("Physics",files.GetOutputPath());

  const double rhof   = physics_map.get("rho0");
  const double Uref   = physics_map.get("Uref");
  const double Lref   = physics_map.get("Lref");
  const double  muf   = physics_map.get("mu0");
  const double MUMHD  = physics_map.get("MUMHD");
  const double SIGMHD = physics_map.get("SIGMHD");
  const double   Bref = physics_map.get("Bref");
  const double sigma  = physics_map.get("sigma");

  const double   _pref = rhof*Uref*Uref;             physics_map.set("pref",_pref);
  const double   _Re  = (rhof*Uref*Lref)/muf;        physics_map.set("Re",_Re);
  const double   _Fr  = (Uref*Uref)/(9.81*Lref);     physics_map.set("Fr",_Fr);
  const double   _Pr=muf/rhof;                       physics_map.set("Pr",_Pr);

  const double   _Rem = MUMHD*SIGMHD*Uref*Lref;      physics_map.set("Rem",_Rem);
  const double   _Hm  = Bref*Lref*sqrt(SIGMHD/muf);  physics_map.set("Hm",_Hm);
  const double   _S   = _Hm*_Hm/(_Re*_Rem);          physics_map.set("S",_S);

  const double   _We  = (Uref*Uref*Lref*rhof)/sigma; physics_map.set("We",_We);

  // ======= Mesh =====
  const unsigned NoLevels = 1;
  const unsigned dim = 3;
  const GeomElType geomel_type = HEX;
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

//================================
// ======= Add QUANTITIES ========
//================================
  VelocityX  velocityX("Qty_Velocity0",qty_map,1,QQ);   qty_map.AddQuantity(&velocityX);
  VelocityY  velocityY("Qty_Velocity1",qty_map,1,QQ);   qty_map.AddQuantity(&velocityY);
  VelocityZ  velocityZ("Qty_Velocity2",qty_map,1,QQ);   qty_map.AddQuantity(&velocityZ);

  VelocityAdjX  velocity_adjX("Qty_VelocityAdj0",qty_map,1,QQ);   qty_map.AddQuantity(&velocity_adjX);
  VelocityAdjY  velocity_adjY("Qty_VelocityAdj1",qty_map,1,QQ);   qty_map.AddQuantity(&velocity_adjY);
  VelocityAdjZ  velocity_adjZ("Qty_VelocityAdj2",qty_map,1,QQ);   qty_map.AddQuantity(&velocity_adjZ);

  MagnFieldHomAdjX    bhom_adjX("Qty_MagnFieldHomAdj0",qty_map,1,QQ);        qty_map.AddQuantity(&bhom_adjX);
  MagnFieldHomAdjY    bhom_adjY("Qty_MagnFieldHomAdj1",qty_map,1,QQ);        qty_map.AddQuantity(&bhom_adjY);
  MagnFieldHomAdjZ    bhom_adjZ("Qty_MagnFieldHomAdj2",qty_map,1,QQ);        qty_map.AddQuantity(&bhom_adjZ);

  MagnFieldHomX bhomX("Qty_MagnFieldHom0",qty_map,1,QQ);     qty_map.AddQuantity(&bhomX);
  MagnFieldHomY bhomY("Qty_MagnFieldHom1",qty_map,1,QQ);     qty_map.AddQuantity(&bhomY);
  MagnFieldHomZ bhomZ("Qty_MagnFieldHom2",qty_map,1,QQ);     qty_map.AddQuantity(&bhomZ);

  MagnFieldExtX BextX("Qty_MagnFieldExt0",qty_map,1,QQ);     qty_map.AddQuantity(&BextX);
  MagnFieldExtY BextY("Qty_MagnFieldExt1",qty_map,1,QQ);     qty_map.AddQuantity(&BextY);
  MagnFieldExtZ BextZ("Qty_MagnFieldExt2",qty_map,1,QQ);     qty_map.AddQuantity(&BextZ);

  MagnFieldHomLagMult         bhom_lag_mult("Qty_MagnFieldHomLagMult",qty_map,1,LL);     qty_map.AddQuantity(&bhom_lag_mult);
  MagnFieldExtLagMult         Bext_lag_mult("Qty_MagnFieldExtLagMult",qty_map,1,LL);     qty_map.AddQuantity(&Bext_lag_mult);
  MagnFieldHomLagMultAdj  bhom_lag_mult_adj("Qty_MagnFieldHomLagMultAdj",qty_map,1,LL);  qty_map.AddQuantity(&bhom_lag_mult_adj);
  Pressure                         pressure("Qty_Pressure",qty_map,1,LL);                qty_map.AddQuantity(&pressure);
  PressureAdj                  pressure_adj("Qty_PressureAdj",qty_map,1,LL);             qty_map.AddQuantity(&pressure_adj);

//consistency check
//  if (bhom._dim !=  Bext._dim)     {std::cout << "main: inconsistency" << std::endl;abort();}
//  if (bhom._FEord !=  Bext._FEord) {std::cout << "main: inconsistency" << std::endl;abort();}
//  if (velocity._dim !=  des_velocity._dim) {std::cout << "main: inconsistency" << std::endl; abort();}
//  if (velocity._FEord !=  des_velocity._FEord) {std::cout << "main: inconsistency" << std::endl; abort();}

//================================
//==== END Add QUANTITIES ========
//================================

  // ====== Start new main =================================
  MultiLevelMesh ml_msh;
  ml_msh.GenerateCoarseBoxMesh(8,8,8,0,1,0,2,0,1,HEX27,"fifth"); //   ml_msh.GenerateCoarseBoxMesh(numelemx,numelemy,numelemz,xa,xb,ya,yb,za,zb,elemtype,"seventh");
  ml_msh.RefineMesh(NoLevels,NoLevels,NULL);
  ml_msh.PrintInfo();

  ml_msh.SetDomain(&mybox);

  MultiLevelSolution ml_sol(&ml_msh);

  ml_sol.AddSolutionVector(ml_msh.GetDimension(),"Qty_MagnFieldHom",LAGRANGE,SECOND,0);
  ml_sol.AddSolutionVector(ml_msh.GetDimension(),"Qty_MagnFieldExt",LAGRANGE,SECOND,0);
  ml_sol.AddSolutionVector(ml_msh.GetDimension(),"Qty_MagnFieldHomAdj",LAGRANGE,SECOND,0);
  ml_sol.AddSolutionVector(ml_msh.GetDimension(),"Qty_Velocity",LAGRANGE,SECOND,0);
  ml_sol.AddSolutionVector(ml_msh.GetDimension(),"Qty_VelocityAdj",LAGRANGE,SECOND,0);

  ml_sol.AddSolutionVector(ml_msh.GetDimension(),"Qty_DesVelocity",LAGRANGE,SECOND,0,false);

  ml_sol.AddSolution("Qty_Pressure",LAGRANGE,FIRST,0);
  ml_sol.AddSolution("Qty_PressureAdj",LAGRANGE,FIRST,0);
  ml_sol.AddSolution("Qty_MagnFieldHomLagMult",LAGRANGE,FIRST,0);
  ml_sol.AddSolution("Qty_MagnFieldExtLagMult",LAGRANGE,FIRST,0);
  ml_sol.AddSolution("Qty_MagnFieldHomLagMultAdj",LAGRANGE,FIRST,0);

  MultiLevelProblem ml_prob(&ml_sol);
  ml_prob.SetMeshTwo(&mesh);
  ml_prob.SetQuadratureRuleAllGeomElems("fifth");
//   ml_prob.SetElemTypeAllDims();
  ml_prob.SetInputParser(&physics_map);
  ml_prob.SetQtyMap(&qty_map);


   // ******* Initial condition *******
  ml_sol.Initialize("Qty_Velocity0",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_Velocity1",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_Velocity2",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_VelocityAdj0",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_VelocityAdj1",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_VelocityAdj2",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_MagnFieldHom0",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_MagnFieldHom1",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_MagnFieldHom2",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_MagnFieldExt0",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_MagnFieldExt1",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_MagnFieldExt2",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_MagnFieldHomAdj0",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_MagnFieldHomAdj1",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_MagnFieldHomAdj2",SetInitialCondition, &ml_prob);

  ml_sol.Initialize("Qty_DesVelocity0",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_DesVelocity1",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_DesVelocity2",SetInitialCondition, &ml_prob);

  ml_sol.Initialize("Qty_Pressure",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_PressureAdj",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_MagnFieldHomLagMult",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_MagnFieldExtLagMult",SetInitialCondition, &ml_prob);
  ml_sol.Initialize("Qty_MagnFieldHomLagMultAdj",SetInitialCondition, &ml_prob);

  // ******* Set boundary function function *******
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

  // ******* Generate boundary conditions *******
  ml_sol.GenerateBdc("Qty_Velocity0","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_Velocity1","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_Velocity2","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_VelocityAdj0","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_VelocityAdj1","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_VelocityAdj2","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_MagnFieldHom0","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_MagnFieldHom1","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_MagnFieldHom2","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_MagnFieldExt0","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_MagnFieldExt1","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_MagnFieldExt2","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_MagnFieldHomAdj0","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_MagnFieldHomAdj1","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_MagnFieldHomAdj2","Steady",&ml_prob);

  ml_sol.GenerateBdc("Qty_Pressure","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_PressureAdj","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_MagnFieldHomLagMult","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_MagnFieldExtLagMult","Steady",&ml_prob);
  ml_sol.GenerateBdc("Qty_MagnFieldHomLagMultAdj","Steady",&ml_prob);

  // ******* Debug *******
  ml_sol.SetWriter(VTK);
  std::vector<std::string> print_vars(1); print_vars[0] = "All"; // we should find a way to make this easier
  ml_sol.GetWriter()->Write(files.GetOutputPath(),"biquadratic",print_vars);

//===============================================
//================== Add EQUATIONS  AND ======================
//========= associate an EQUATION to QUANTITIES ========
//========================================================

#if NS_EQUATIONS==1
  SystemTwo & eqnNS = ml_prob.add_system<SystemTwo>("Eqn_NS");

          eqnNS.AddSolutionToSystemPDEVector(ml_msh.GetDimension(),"Qty_Velocity");
	  eqnNS.AddSolutionToSystemPDE("Qty_Pressure");

          eqnNS.AddUnknownToSystemPDE(&velocityX);
          eqnNS.AddUnknownToSystemPDE(&velocityY);
          eqnNS.AddUnknownToSystemPDE(&velocityZ);
          eqnNS.AddUnknownToSystemPDE(&pressure);

          eqnNS.SetAssembleFunction(GenMatRhsNS);
#endif

#if NSAD_EQUATIONS==1
  SystemTwo & eqnNSAD = ml_prob.add_system<SystemTwo>("Eqn_NSAD");

            eqnNSAD.AddSolutionToSystemPDEVector(ml_msh.GetDimension(),"Qty_VelocityAdj");
            eqnNSAD.AddSolutionToSystemPDE("Qty_PressureAdj");

            eqnNSAD.AddUnknownToSystemPDE(&velocity_adjX);
            eqnNSAD.AddUnknownToSystemPDE(&velocity_adjY);
            eqnNSAD.AddUnknownToSystemPDE(&velocity_adjZ);
            eqnNSAD.AddUnknownToSystemPDE(&pressure_adj);

            eqnNSAD.SetAssembleFunction(GenMatRhsNSAD);
#endif

#if MHD_EQUATIONS==1
  SystemTwo & eqnMHD = ml_prob.add_system<SystemTwo>("Eqn_MHD");

           eqnMHD.AddSolutionToSystemPDEVector(ml_msh.GetDimension(),"Qty_MagnFieldHom");
           eqnMHD.AddSolutionToSystemPDE("Qty_MagnFieldHomLagMult");

           eqnMHD.AddUnknownToSystemPDE(&bhomX);
           eqnMHD.AddUnknownToSystemPDE(&bhomY);
           eqnMHD.AddUnknownToSystemPDE(&bhomZ);
           eqnMHD.AddUnknownToSystemPDE(&bhom_lag_mult);

           eqnMHD.SetAssembleFunction(GenMatRhsMHD);
#endif

#if MHDAD_EQUATIONS==1
  SystemTwo & eqnMHDAD = ml_prob.add_system<SystemTwo>("Eqn_MHDAD");

             eqnMHDAD.AddSolutionToSystemPDEVector(ml_msh.GetDimension(),"Qty_MagnFieldHomAdj");
             eqnMHDAD.AddSolutionToSystemPDE("Qty_MagnFieldHomLagMultAdj");

             eqnMHDAD.AddUnknownToSystemPDE(&bhom_adjX);
             eqnMHDAD.AddUnknownToSystemPDE(&bhom_adjY);
             eqnMHDAD.AddUnknownToSystemPDE(&bhom_adjZ);
             eqnMHDAD.AddUnknownToSystemPDE(&bhom_lag_mult_adj);

             eqnMHDAD.SetAssembleFunction(GenMatRhsMHDAD);
#endif

#if MHDCONT_EQUATIONS==1
  SystemTwo & eqnMHDCONT = ml_prob.add_system<SystemTwo>("Eqn_MHDCONT");

               eqnMHDCONT.AddSolutionToSystemPDEVector(ml_msh.GetDimension(),"Qty_MagnFieldExt");
               eqnMHDCONT.AddSolutionToSystemPDE("Qty_MagnFieldExtLagMult");

               eqnMHDCONT.AddUnknownToSystemPDE(&BextX);
               eqnMHDCONT.AddUnknownToSystemPDE(&BextY);
               eqnMHDCONT.AddUnknownToSystemPDE(&BextZ);
	       eqnMHDCONT.AddUnknownToSystemPDE(&Bext_lag_mult);

               eqnMHDCONT.SetAssembleFunction(GenMatRhsMHDCONT);
#endif

//================================
//========= End add EQUATIONS  and ========
//========= associate an EQUATION to QUANTITIES ========
//================================

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

  // ======== OptLoop ===================================
  FemusInputParser<double> loop_map("TimeLoop",files.GetOutputPath());
  OptLoop opt_loop(files,loop_map);

  opt_loop.TransientSetup(ml_prob);  // reset the initial state (if restart) and print the Case   /*TODO fileIO */

  opt_loop.optimization_loop(ml_prob);

//============= prepare default for next restart ==========
// at this point, the run has been completed
// well, we do not know whether for the whole time range time.N-M.xmf
// or the optimization loop stopped before
// we could also print the last step number
  files.PrintRunForRestart(DEFAULT_LAST_RUN); /*(iproc==0)*/
  files.log_petsc();

// ============  clean ================================
  ml_prob.clear();
  mesh.clear();

  return 0;

}

//C++ includes 
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <cstdlib>
#include <sstream>

// External library include 
#include "FEMTTUConfig.h"

// FEMuS
#include "paral.hpp"
#include "FemusDefault.hpp"
#include "FemusInit.hpp"
#include "Files.hpp"
#include "GeomEl.hpp"
#include "MultiLevelMeshTwo.hpp"
#include "GenCase.hpp"
#include "FETypeEnum.hpp"
#include "MultiLevelProblemTwo.hpp"
#include "ElemType.hpp"
#include "TimeLoop.hpp"
#include "Typedefs.hpp"
#include "CmdLine.hpp"
#include "Quantity.hpp"
#include "QTYnumEnum.hpp"
#include "Box.hpp"  //for the DOMAIN


// application includes
#include "Opt_conf.hpp"
#include "OptLoop.hpp"
#include "OptQuantities.hpp"
#include "EqnNS.hpp"
#include "EqnNSAD.hpp"
#include "EqnMHD.hpp"
#include "EqnMHDAD.hpp"
#include "EqnMHDCONT.hpp"

#ifdef HAVE_LIBMESH
#include "libmesh/libmesh.h"
#endif

using namespace femus;


// double funzione(double t , const double* xyz) {return 1.;} 

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
  FemusInputParser<double> physics_map("Physics",files._output_path);
  
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
  FemusInputParser<double> mesh_map("Mesh",files._output_path);
    GenCase mesh(files,mesh_map,"straightQ3D2x2x2ZERO.gam");
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
  
  // =======FEElems =====  //remember to delete the FE at the end
  const std::string  FEFamily[QL] = {"biquadratic","linear","constant"}; 
  std::vector< std::vector<elem_type*> > FEElemType_vec(mesh.get_dim());
  for (int idim=0;idim < mesh.get_dim(); idim++)   FEElemType_vec[idim].resize(QL);
  
  for (int idim=0;idim < mesh.get_dim(); idim++) { 
    for (int fe=0; fe<QL; fe++) {
       FEElemType_vec[idim][fe] = elem_type::build(mesh.GetGeomEl(idim,mesh._mesh_order)._geomel_id.c_str(),fe,
						            qrule[idim].GetGaussOrderString().c_str());
       FEElemType_vec[idim][fe]->EvaluateShapeAtQP(mesh.GetGeomEl(idim,mesh._mesh_order)._geomel_id.c_str(),fe);
     }
    }  
  
  // ===== QuantityMap =========================================
  QuantityMap  qty_map(mesh,&physics_map);

//================================
// ======= Add QUANTITIES ========  
//================================
  MagnFieldHom bhom("Qty_MagnFieldHom",qty_map,mesh.get_dim(),FE_MAGNFIELDHOM);     qty_map.set_qty(&bhom);  
  MagnFieldExt Bext("Qty_MagnFieldExt",qty_map,mesh.get_dim(),FE_MAGNFIELDEXT);     qty_map.set_qty(&Bext);  

//consistency check
 if (bhom._dim !=  Bext._dim)     {std::cout << "main: inconsistency" << std::endl;abort();}
 if (bhom._FEord !=  Bext._FEord) {std::cout << "main: inconsistency" << std::endl;abort();}

 MagnFieldHomLagMult         bhom_lag_mult("Qty_MagnFieldHomLagMult",qty_map,1,FE_MAGNFIELDHOMLAGMULT);     qty_map.set_qty(&bhom_lag_mult);
 MagnFieldExtLagMult         Bext_lag_mult("Qty_MagnFieldExtLagMult",qty_map,1,FE_MAGNFIELDEXTLAGMULT);     qty_map.set_qty(&Bext_lag_mult);
 MagnFieldHomAdj                  bhom_adj("Qty_MagnFieldHomAdj",qty_map,mesh.get_dim(),FE_MAGNFIELDHOM);        qty_map.set_qty(&bhom_adj);
 MagnFieldHomLagMultAdj  bhom_lag_mult_adj("Qty_MagnFieldHomLagMultAdj",qty_map,1,FE_MAGNFIELDHOMLAGMULT);  qty_map.set_qty(&bhom_lag_mult_adj);

  Pressure  pressure("Qty_Pressure",qty_map,1,FE_PRESSURE);            qty_map.set_qty(&pressure);
  Velocity  velocity("Qty_Velocity",qty_map,mesh.get_dim(),FE_VELOCITY);   qty_map.set_qty(&velocity);  

  VelocityAdj  velocity_adj("Qty_VelocityAdj",qty_map,mesh.get_dim(),FE_VELOCITY);         qty_map.set_qty(&velocity_adj);  
  PressureAdj pressure_adj("Qty_PressureAdj",qty_map,1,FE_PRESSURE);                  qty_map.set_qty(&pressure_adj);
  DesVelocity des_velocity("Qty_DesVelocity",qty_map,mesh.get_dim(),FE_DESVELOCITY);       qty_map.set_qty(&des_velocity);
 
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

  // ====== MultiLevelProblemTwo =================================
  MultiLevelProblemTwo equations_map(files,physics_map,qty_map,mesh,FEElemType_vec,qrule);
  
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

  // ======== OptLoop ===================================
  OptLoop opt_loop(files); 
           opt_loop._timemap.read();
           opt_loop._timemap.print();

  opt_loop.TransientSetup(equations_map);  // reset the initial state (if restart) and print the Case   /*TODO fileIO */ 

  opt_loop.optimization_loop(equations_map);
    
// // //   eqnNS->FunctionIntegral (0,funzione);
// // //   eqnNS->FunctionIntegral (1,funzione);

//============= prepare default for next restart ==========  
// at this point, the run has been completed 
// well, we do not know whether for the whole time range time.N-M.xmf
// or the optimization loop stopped before
// we could also print the last step number
  files.PrintRunForRestart(DEFAULT_LAST_RUN); /*(iproc==0)*/ 
  files.log_petsc();

// ============  clean ================================
  equations_map.clean();
  mesh.clear();
  
  return 0;
  
}

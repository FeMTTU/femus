
#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "VTKWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "Marker.hpp"
#include "Line.hpp"

#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"

#include "NumericVector.hpp"
#include "adept.h"

#include "../include/mpmFem.hpp"


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
//double NeumannFactor = 0.;

using namespace femus;

// OLD BEST RESULT WITH E = 4.2 * 1.e6, 5 levels, dt= 0.01, NR = 300, R0 = 1.5, factor = 1.3
// MOST BEST RESULT WITH E = 4.2 * 1.e6, 4 levels, dt= 0.01, NR = 300, R0 = 1.4, factor = 1.14,  beta = 0.3, Gamma = 0.5

double DT = 0.01;

double SetVariableTimeStep(const double time)
{
  double dt =  DT;
  return dt;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int facename, const double time)
{
    
  bool test = 1; //dirichlet
  value = 0.;

  if(!strcmp(name, "DY")) {
    if(2 == facename || 3 == facename || 4 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "DX")) {
    if(2 == facename || 3 == facename || 4 == facename) {
      test = 0;
      value = 0;
    }
  }
  
  if(!strcmp(name, "VY")) {
    if(2 == facename || 3 == facename || 4 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "VX")) {
    if(2 == facename || 3 == facename || 4 == facename ) {
      test = 0;
      value = 0;
    }
  }

  return test;

}

int main(int argc, char** args)
{
  
  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  std::string soft = "soft";
  std::string  stiff = "stiff";
  
  std::string  material[2] ={soft, stiff};
  
  std::map < std::pair < std::string, unsigned > , double > YC; 
  
//   YC[std::make_pair (soft, 2u)] = 0.2;  
//   YC[std::make_pair (soft, 3u)] = 0.1;
//   YC[std::make_pair (soft, 4u)] = 0.075;
//   YC[std::make_pair (soft, 5u)] = 0.05;
// 
//   YC[std::make_pair (stiff, 2u)] = 0.3;
//   YC[std::make_pair (stiff, 3u)] = 0.15;
//   YC[std::make_pair (stiff, 4u)] = 0.09;
//   YC[std::make_pair (stiff, 5u)] = 0.05;

  YC[std::make_pair (soft, 2u)] = 0.2;  
  YC[std::make_pair (soft, 3u)] = 0.2 / 2.;
  YC[std::make_pair (soft, 4u)] = 0.2 / 4.;
  YC[std::make_pair (soft, 5u)] = 0.2 / 8.;

  YC[std::make_pair (stiff, 2u)] = 0.3;  
  YC[std::make_pair (stiff, 3u)] = 0.3 / 2.;
  YC[std::make_pair (stiff, 4u)] = 0.3 / 4.;
  YC[std::make_pair (stiff, 5u)] = 0.3 / 8.;    
  
//   std::map < std::pair < std::string, unsigned > , double > NF; 
//   
//   NF[std::make_pair (soft, 2u)] = 0.;
//   NF[std::make_pair (soft, 3u)] = 0.;
//   NF[std::make_pair (soft, 4u)] = 0.;
//   NF[std::make_pair (soft, 5u)] = 0.3;
// 
//   NF[std::make_pair (stiff, 2u)] = 0.3;
//   NF[std::make_pair (stiff, 3u)] = 0.3;
//   NF[std::make_pair (stiff, 4u)] = 0.3;
//   NF[std::make_pair (stiff, 5u)] = 0.6;  
  
  std::map < std::pair < std::string, unsigned > , double > SF1; 
  
  std::map < std::pair < std::string, unsigned > , double > SF2; 
  
  SF1[std::make_pair (soft, 2u)] = 1.e-3;
  SF1[std::make_pair (soft, 3u)] = 1.e-3;
  SF1[std::make_pair (soft, 4u)] = 1.e-3;
  SF1[std::make_pair (soft, 5u)] = 1.e-3;
    
  SF2[std::make_pair (soft, 2u)] = 1.e-7;
  SF2[std::make_pair (soft, 3u)] = 1.e-7;
  SF2[std::make_pair (soft, 4u)] = 1.e-7;
  SF2[std::make_pair (soft, 5u)] = 1.e-7;
  
  SF1[std::make_pair (stiff, 2u)] = 1.e-5;
  SF1[std::make_pair (stiff, 3u)] = 1.e-5;
  SF1[std::make_pair (stiff, 4u)] = 1.e-5;
  SF1[std::make_pair (stiff, 5u)] = 1.e-5;  

  SF2[std::make_pair (stiff, 2u)] = 1.e-9;
  SF2[std::make_pair (stiff, 3u)] = 1.e-9;
  SF2[std::make_pair (stiff, 4u)] = 1.e-9;
  SF2[std::make_pair (stiff, 5u)] = 1.e-9;  
  
  std::map < std::pair < std::string, unsigned > , double > YM; 
  
  YM[std::make_pair (soft, 2u)] = 4.2 * 1.e6;
  YM[std::make_pair (soft, 3u)] = 4.2 * 1.e6;
  YM[std::make_pair (soft, 4u)] = 4.2 * 1.e6;
  YM[std::make_pair (soft, 5u)] = 4.2 * 1.e6;

  YM[std::make_pair (stiff, 2u)] = 4.2 * 1.e8;
  YM[std::make_pair (stiff, 3u)] = 4.2 * 1.e8;
  YM[std::make_pair (stiff, 4u)] = 4.2 * 1.e8;
  YM[std::make_pair (stiff, 5u)] = 4.2 * 1.e8;  
  
  std::pair <std::string, unsigned > simulation;
  
  unsigned timestep0 = 0;
  unsigned n_timesteps = 250 + timestep0;
    
  std::vector < std::map < std::pair < std::string, unsigned > , double > > CM(n_timesteps +1 - timestep0); 
    
  unsigned mat0 = 0;
  unsigned matn = 2;
  
  unsigned nl0 = 2;
  unsigned nln = 6;
  
  for(unsigned mat = mat0; mat< matn; mat++){
    for(unsigned nl = nl0; nl < nln; nl++) {
      
      simulation = std::make_pair (material[mat], nl);
             
      unsigned numberOfUniformLevels = simulation.second; 
      unsigned numberOfSelectiveLevels = 0;

      double Lref = 1.;
      double Uref = 1.;
      double rhos = 1000;
      double nu = 0.4;
      double E = YM[simulation];
  
      beta = 0.3; //was 0.25 
      Gamma = 0.5;

      Parameter par(Lref, Uref);

      // Generate Solid Object
      Solid solid;
      solid = Solid(par, E, nu, rhos, "Neo-Hookean");
            
      MultiLevelMesh mlMsh;
      mlMsh.ReadCoarseMesh("../input/inclined_plane_2D.neu", "fifth", 1.);
      mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);

      mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);
      numberOfUniformLevels = 1;

      unsigned dim = mlMsh.GetDimension();
      
      std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<dim<<std::endl;

      MultiLevelSolution mlSol(&mlMsh);
      // add variables to mlSol
 
      mlSol.AddSolution("DX", LAGRANGE, SECOND, 2);
      if(dim > 1) mlSol.AddSolution("DY", LAGRANGE, SECOND, 2);
      if(dim > 2) mlSol.AddSolution("DZ", LAGRANGE, SECOND, 2);

//       mlSol.AddSolution("VX", LAGRANGE, SECOND, 2, false);
//       if(dim > 1) mlSol.AddSolution("VY", LAGRANGE, SECOND, 2, false);
//       if(dim > 2) mlSol.AddSolution("VZ", LAGRANGE, SECOND, 2, false);
            
      mlSol.AddSolution("M", LAGRANGE, SECOND, 2);
      mlSol.AddSolution("Mat", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);
      
      //mlSol.AddSolution("NF", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

      mlSol.Initialize("All");

      mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

      // ******* Set boundary conditions *******
      mlSol.GenerateBdc("DX", "Steady");
      if(dim > 1) mlSol.GenerateBdc("DY", "Steady");
      if(dim > 2) mlSol.GenerateBdc("DZ", "Steady");
      
//       mlSol.GenerateBdc("VX", "Steady");
//       if(dim > 1) mlSol.GenerateBdc("VY", "Steady");
//       if(dim > 2) mlSol.GenerateBdc("VZ", "Steady");
      
      mlSol.GenerateBdc("M", "Steady");

      MultiLevelProblem ml_prob(&mlSol);

      ml_prob.parameters.set<Solid> ("SolidMPM") = solid;
      ml_prob.parameters.set<Solid> ("SolidFEM") = solid;

      // ******* Add MPM system to the MultiLevel problem *******
      TransientNonlinearImplicitSystem& system = ml_prob.add_system < TransientNonlinearImplicitSystem > ("MPM_FEM");
      system.AddSolutionToSystemPDE("DX");
      if(dim > 1)system.AddSolutionToSystemPDE("DY");
      if(dim > 2) system.AddSolutionToSystemPDE("DZ");
      
//       system.AddSolutionToSystemPDE("VX");
//       if(dim > 1)system.AddSolutionToSystemPDE("VY");
//       if(dim > 2) system.AddSolutionToSystemPDE("VZ");

      // ******* System MPM Assembly *******
      system.SetAssembleFunction(AssembleMPMSys);
      //system.SetAssembleFunction(AssembleFEM);
      // ******* set MG-Solver *******
      system.SetMgType(V_CYCLE);


      system.SetAbsoluteLinearConvergenceTolerance(1.0e-10);
      system.SetMaxNumberOfLinearIterations(1);
      system.SetNonLinearConvergenceTolerance(1.e-9);
      system.SetMaxNumberOfNonLinearIterations(20);

      system.SetNumberPreSmoothingStep(1);
      system.SetNumberPostSmoothingStep(1);

      // ******* Set Preconditioner *******
      system.SetLinearEquationSolverType(FEMuS_DEFAULT);

      system.init();

      // ******* Set Smoother *******
      system.SetSolverFineGrids(GMRES);

      system.SetPreconditionerFineGrids(ILU_PRECOND);

      system.SetTolerances(1.e-10, 1.e-15, 1.e+50, 40, 40);


      //BEGIN init particles   
      unsigned size = 1;
      std::vector < std::vector < double > > x; // marker
      double xc = .0;
      double yc = YC[simulation];
  
      x.resize(size);
      x[0].resize(dim, 0.);
      x[0][0] = xc;
      x[0][1] = yc;
 
      double R = 1.6;
      double R0 = 1.4; 
  
      double PI = acos(-1.);
      unsigned NR = 300;
      unsigned NL = NR / (2 * PI);
      double DL = R0 / NL;

      for(unsigned i = 0; i < NL; i++) {
        double  r = R0 - i * DL;
        unsigned Nr = static_cast <unsigned>(ceil(NR * r / R0));
        double dtheta = 2 * PI / Nr;
        unsigned sizeOld = x.size();
        x.resize(sizeOld + Nr);
        for(unsigned s = sizeOld; s < x.size(); s++) {
          x[s].resize(dim);
        }
        for(unsigned j = 0; j < Nr; j++) {
          x[sizeOld + j][0] = xc + r * cos(j * dtheta);
          x[sizeOld + j][1] = yc + r * sin(j * dtheta);
        }
      }
      double MASS = PI * R0 * R0 * rhos;
      size = x.size();
      std::vector < double > mass(x.size(), MASS / x.size()); // uniform marker volume
  
      if( fabs(R-R0) > 1.0e-10 ) {
    
        double factor = 1.14; 
        unsigned NL = getNumberOfLayers((R-R0)/DL, factor);
        std::cout << NL <<std::endl;
      
        double  r = R0;
        for(unsigned i = 1; i <= NL; i++) {
          DL = DL / factor;
          r += DL;
          NR = static_cast <unsigned>(ceil (NR * factor) );
          double dtheta = 2 * PI / NR;
          unsigned sizeOld = x.size();
          x.resize(sizeOld + NR);
          for(unsigned s = sizeOld; s < x.size(); s++) {
            x[s].resize(dim);
          }
          for(unsigned j = 0; j < NR; j++) {
            x[sizeOld + j][0] = xc + r * cos(j * dtheta);
            x[sizeOld + j][1] = yc + r * sin(j * dtheta);
          }
          mass.resize(x.size(), rhos * r * dtheta * DL);
        }
        size = x.size();
      }
  
//       double totalMass = 0;
//       for(unsigned i = 0; i < mass.size(); i++){
//         totalMass += mass[i];
//       }
//       std::cout << totalMass<<" "<< rhos * PI * R * R << std::endl;
  
      std::vector < MarkerType > markerType(size,VOLUME);
      
      unsigned solType = 2;
      linea = new Line(x, mass, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);
      
      
      //END init particles 
      
      //BEGIN Print solution 0 *******
      
      std::ostringstream outputFolder;
      outputFolder << "outputE" << simulation.first  << "Level"<<simulation.second;
      
      mlSol.SetWriter(VTK);
      
      std::vector<std::string> mov_vars;
      mov_vars.push_back("DX");
      mov_vars.push_back("DY");
      //mov_vars.push_back("DZ");
      mlSol.GetWriter()->SetMovingMesh(mov_vars);
      
      std::vector<std::string> print_vars;
      print_vars.push_back("All");
      
      
      //mlSol.GetWriter()->SetDebugOutput(true);
      mlSol.GetWriter()->Write(outputFolder.str(), "biquadratic", print_vars, 0);
      
      std::vector < std::vector < std::vector < double > > > line(1);
      
      linea->GetLine(line[0]);
      PrintLine(outputFolder.str(), "line", line, 0);
      
      //END Print solution 0 *******
      
      linea->GetParticlesToGridMaterial();
            
      system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
     
      
      std::ostringstream filename;
      filename <<outputFolder.str()<< "/centerOfMass.txt";
      std::ofstream fout;
      fout.open( filename.str().c_str() );
      fout << "iteration, time, analytic, ";
      fout << "E" << simulation.first  << "Level"<<simulation.second << std::endl;
      fout << "0, 0., 0., 0." << std::endl;
      
      for(unsigned time_step = 1; time_step <= n_timesteps; time_step++) {

        double theta = PI / 4;
        gravity[0] =  9.81 * sin(theta);
        gravity[1] = -9.81 * cos(theta);  
        scalingFactor1 = SF1[simulation];
        scalingFactor2 = SF2[simulation];
    
//         //NeumannFactor = NF[simulation];
//     
//         if ( time_step <= timestep0 ){
//           gravity[0] = 0.;
//           NeumannFactor = 0.;
//         }
//         
//         //SetNeumannFactor(&mlSol);
        
        
        system.MGsolve();

        // ******* Print solution *******
        mlSol.GetWriter()->Write(outputFolder.str(), "biquadratic", print_vars, time_step);

        GridToParticlesProjection(ml_prob, *linea);
        linea->GetLine(line[0]); 
        
        PrintLine(outputFolder.str(), "line", line, time_step);

        if(time_step >= timestep0){
          CM[time_step - timestep0][simulation] = line[0][0][0];  
          
          unsigned it = time_step - timestep0;
          double time = it * DT;
          fout << it <<", "<< time <<", "<<1./3.*9.81*time*time*sqrt(2.)/2. << ", ";
          fout << CM[it][simulation] <<std::endl;
       }
    }  
    fout.close();
    delete linea;
    }
  }
  
  std::ostringstream filename;
  filename << "./centerOfMass.txt";
  std::ofstream fout;
  fout.open( filename.str().c_str() );
  fout << "iteration, time, analytic, ";
  for(unsigned mat = mat0; mat < matn; mat++){
    for(unsigned nl = nl0; nl < nln; nl++) {
      simulation = std::make_pair (material[mat], nl);
      fout << "E" << simulation.first  << "Level"<<simulation.second<<", "; 
    }
  }
  fout << std::endl;
  for(unsigned it = 0;it < CM.size();it++){
    double time = it * DT;
    fout << it <<", "<< time <<", "<<1./3.*9.81*time*time*sqrt(2.)/2. << ", ";
    for(unsigned mat = mat0; mat< matn; mat++){
      for(unsigned nl = nl0; nl < nln; nl++) {
        simulation = std::make_pair (material[mat], nl);
        fout << CM[it][simulation] <<", ";
      }
    }
    fout << std::endl;
  }
  
  fout.close();
  
  return 0;

} //end main



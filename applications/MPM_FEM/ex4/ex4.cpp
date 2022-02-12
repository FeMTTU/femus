
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

//double yMin;

using namespace femus;

double DT = 0.001;
unsigned numberOfUniformLevels0;

double GetTimeStep(const double &y){
  
  double dt;
  double DY =  0.184 / pow(2,numberOfUniformLevels0-1);
  double y0 = -1.52 + DY;
  
  std::cout << y0 << std::endl;
  double DT0 = 0.001 / pow(2,numberOfUniformLevels0-1);
  
  double y1 = -0.16;
  double DT1 = 0.005;
  
  if(y < y0){
    dt = DT0;
  }
  else if(y < y1){
    const double y2 = 0.5 * (y0 + y1);
    const double DT2 = 0.01;
    double f0 = (y - y1) / (y0 - y1) * (y - y2) / (y0 - y2);
    double f1 = (y - y0) / (y1 - y0) * (y - y2) / (y1 - y2);
    double f2 = (y - y0) / (y2 - y0) * (y - y1) / (y2 - y1);
   
    dt = f0 * DT0 + f1 * DT1 + f2 * DT2;
  }
  else {
    dt = DT1;
  }
  return dt;
  
}

double SetVariableTimeStep(const double time)
{
  return DT;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  if(1 == facename) test = 0;
  
  if(!strcmp(name,"DY") && 3 == facename){
    test = 0;    
  }

  return test;

}

int main(int argc, char** args)
{

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  
  std::string stiff = "stiff";
  std::string soft = "soft";
  
  //std::string material = stiff;
  std::string material = stiff;
  
  std::map < std::pair < std::string, unsigned > , double > SF1; 
  std::map < std::pair < std::string, unsigned > , double > SF2; 
  std::map < std::pair < std::string, unsigned > , double > NF; 
  
  //stiff soft matrix
  
  SF1[std::make_pair (stiff, 1u)] = 1.e-2;
  SF1[std::make_pair (stiff, 2u)] = 1.e-2;
  SF1[std::make_pair (stiff, 3u)] = 1.e-2;
  SF1[std::make_pair (stiff, 4u)] = 1.e-2;
  
  SF2[std::make_pair (stiff, 1u)] = 1.e-6;
  SF2[std::make_pair (stiff, 2u)] = 1.e-6;
  SF2[std::make_pair (stiff, 3u)] = 1.e-6;
  SF2[std::make_pair (stiff, 4u)] = 1.e-6;
  
  NF[std::make_pair (stiff, 1u)] = 0.14; //0.13 very little little
  NF[std::make_pair (stiff, 2u)] = 0.06;  //0.07 very little up
  NF[std::make_pair (stiff, 3u)] = 0.1; //0.1 very little up 
  NF[std::make_pair (stiff, 4u)] = 0.2; //0.2 very little up  
  
  //soft soft matrix
  
  SF1[std::make_pair (soft, 1u)] = 1.e-5;
  SF1[std::make_pair (soft, 2u)] = 1.e-5;
  SF1[std::make_pair (soft, 3u)] = 1.e-5;  
  SF1[std::make_pair (soft, 4u)] = 1.e-5;  
    
  SF2[std::make_pair (soft, 1u)] = 1.e-9;
  SF2[std::make_pair (soft, 2u)] = 1.e-9;
  SF2[std::make_pair (soft, 3u)] = 1.e-9;  
  SF2[std::make_pair (soft, 4u)] = 1.e-9;  
    
  NF[std::make_pair (soft, 1u)] = 0.;
  NF[std::make_pair (soft, 2u)] = 0.;
  NF[std::make_pair (soft, 3u)] = 0.1;
  NF[std::make_pair (soft, 4u)] = 0.2;
      
  unsigned n_timesteps = 3500;
  
  std::vector < std::map < std::pair < std::string, unsigned > , double > > CM(n_timesteps + 1);
  std::vector < std::map < std::pair < std::string, unsigned > , double > > time(n_timesteps + 1); 
  
  unsigned nli[3]={2,3,4};
  for(unsigned kk = 0; kk < 3; kk++){
      
    unsigned nl = nli[kk];  
      
    std::pair <std::string, unsigned > simulation = std::make_pair (material, nl);
    
    scalingFactor1 = SF1[simulation];
    scalingFactor2 = SF2[simulation];
  
    NeumannFactor = NF[simulation];
  
    MultiLevelMesh mlMsh;
    double scalingFactor = 1.;
    unsigned numberOfUniformLevels = simulation.second; //for refinement in 3D
    numberOfUniformLevels0 = numberOfUniformLevels;
    //unsigned numberOfUniformLevels = 1;
    unsigned numberOfSelectiveLevels = 0;

    double Lref = 1.;
    double Uref = 1.;

    //initialize parameters for rolling ball (MPM)
    double rho_MPM = 1000.;
    double nu_MPM = 0.4;
    double E_MPM = 5.91 * 1.e6;

    //initialize parameters for plate (FEM)
    double rho_FEM = 1000.;
    double nu_FEM = 0.4;
    double E_FEM = 4.2 * 1.e7;

    beta = 0.3; //was 0.25
    Gamma = 0.5;

    Parameter par(Lref, Uref);

    // Generate Solid Object
    Solid solidMPM;
    Solid solidFEM;

    solidMPM = Solid(par, E_MPM, nu_MPM, rho_MPM, "Neo-Hookean");
    solidFEM = Solid(par, E_FEM, nu_FEM, rho_FEM, "Neo-Hookean");

    mlMsh.ReadCoarseMesh("../input/basket.neu", "fifth", scalingFactor);
    mlMsh.RefineMesh(numberOfUniformLevels + numberOfSelectiveLevels, numberOfUniformLevels , NULL);

    mlMsh.EraseCoarseLevels(numberOfUniformLevels - 1);
    numberOfUniformLevels = 1;

    unsigned dim = mlMsh.GetDimension();

    MultiLevelSolution mlSol(&mlMsh);
    // add variables to mlSol
    mlSol.AddSolution("DX", LAGRANGE, SECOND, 2);
    if(dim > 1) mlSol.AddSolution("DY", LAGRANGE, SECOND, 2);
    if(dim > 2) mlSol.AddSolution("DZ", LAGRANGE, SECOND, 2);

    mlSol.AddSolution("VX", LAGRANGE, SECOND, 2);
    if(dim > 1) mlSol.AddSolution("VY", LAGRANGE, SECOND, 2);
    if(dim > 2) mlSol.AddSolution("VZ", LAGRANGE, SECOND, 2);

    mlSol.AddSolution("AX", LAGRANGE, SECOND, 2);
    if(dim > 1) mlSol.AddSolution("AY", LAGRANGE, SECOND, 2);
    if(dim > 2) mlSol.AddSolution("AZ", LAGRANGE, SECOND, 2);

    mlSol.AddSolution("M", LAGRANGE, SECOND, 2);
    mlSol.AddSolution("Mat", DISCONTINUOUS_POLYNOMIAL, ZERO, 0, false);

    mlSol.Initialize("All");

    mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);

    mlSol.SetIfFSI(true);

    // ******* Set boundary conditions *******
    mlSol.GenerateBdc("DX", "Steady");
    if(dim > 1) mlSol.GenerateBdc("DY", "Steady");
    if(dim > 2) mlSol.GenerateBdc("DZ", "Steady");
    mlSol.GenerateBdc("M", "Steady");

    MultiLevelProblem ml_prob(&mlSol);

    ml_prob.parameters.set<Solid> ("SolidMPM") = solidMPM;
    ml_prob.parameters.set<Solid> ("SolidFEM") = solidFEM;

    // ******* Add MPM system to the MultiLevel problem *******
    TransientNonlinearImplicitSystem& system = ml_prob.add_system < TransientNonlinearImplicitSystem > ("MPM_FEM");
    system.AddSolutionToSystemPDE("DX");
    if(dim > 1)system.AddSolutionToSystemPDE("DY");
    if(dim > 2) system.AddSolutionToSystemPDE("DZ");

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

    // ******* Add MPM system to the MultiLevel problem *******
    NonLinearImplicitSystem& system2 = ml_prob.add_system < NonLinearImplicitSystem > ("DISP");
    system2.AddSolutionToSystemPDE("DX");
    if(dim > 1)system2.AddSolutionToSystemPDE("DY");
    if(dim > 2) system2.AddSolutionToSystemPDE("DZ");
      
    // ******* System MPM Assembly *******
    system2.SetAssembleFunction(AssembleSolidDisp);
    //system2.SetAssembleFunction(AssembleFEM);
    // ******* set MG-Solver *******
    system2.SetMgType(V_CYCLE);
      
      
    system2.SetAbsoluteLinearConvergenceTolerance(1.0e-10);
    system2.SetMaxNumberOfLinearIterations(1);
    system2.SetNonLinearConvergenceTolerance(1.e-9);
    system2.SetMaxNumberOfNonLinearIterations(1);
      
    system2.SetNumberPreSmoothingStep(1);
    system2.SetNumberPostSmoothingStep(1);
      
    // ******* Set Preconditioner *******
    system2.SetLinearEquationSolverType(FEMuS_DEFAULT);
      
    system2.init();
      
    // ******* Set Smoother *******
    system2.SetSolverFineGrids(GMRES);
      
    system2.SetPreconditionerFineGrids(ILU_PRECOND);
      
    system2.SetTolerances(1.e-10, 1.e-15, 1.e+50, 2, 2);
    

    //BEGIN init particles
    unsigned size = 1;
    std::vector < std::vector < double > > x; // marker
    double yc = 0.;  // FOR E = 4.2 * 1.e8 --> 0.115 (for 3 refinements) 0.09 (for 4) and 0.05  (for 5, this one maybe to be changed)
    // FOR E = 4.2 * 1.e6 --> 0.1. (for 3 refinements) 0.075 (for 4) and 0.05  (for 5)

    x.resize(size);
    x[0].resize(dim, 0.);
    x[0][1] = yc;

    double R = 0.16;
    double R0 = 0.14;

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
        x[sizeOld + j][0] = r * cos(j * dtheta);
        x[sizeOld + j][1] = yc + r * sin(j * dtheta);
      }
    }
    double MASS = PI * R0 * R0 * rho_MPM;
    size = x.size();
    std::vector < double > mass(x.size(), MASS / x.size()); // uniform marker volume

    if( fabs(R - R0) > 1.0e-10 ) {

      double factor = 1.14;
      unsigned NL = getNumberOfLayers((R - R0) / DL, factor);
      //std::cout << NL << std::endl;

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
          x[sizeOld + j][0] = r * cos(j * dtheta);
          x[sizeOld + j][1] = yc + r * sin(j * dtheta);
        }
        mass.resize(x.size(), rho_MPM * r * dtheta * DL);
      }
      size = x.size();
      std::cout<< "Number of Particles = " << x.size() << " Number of Layers = "<< NL <<std::endl;
    }

    double totalMass = 0;
    for(unsigned i = 0; i < mass.size(); i++) {
      totalMass += mass[i];
    }

    std::cout << totalMass << " " << rho_MPM * PI * R * R << std::endl;

    //return 1;

    std::vector < MarkerType > markerType;
    markerType.resize(size);

    for(unsigned j = 0; j < size; j++) {
      markerType[j] = VOLUME;
    }
   
    
    unsigned solType = 2;
    linea = new Line(x, mass, markerType, mlSol.GetLevel(numberOfUniformLevels - 1), solType);
    
    linea->GetParticlesToGridMaterial();

    //END init particles

    std::ostringstream outputFolder;
    outputFolder << "outputE" << simulation.first  << "Level"<<simulation.second;
    
    // ******* Print solution *******
    mlSol.SetWriter(VTK);

    std::vector<std::string> mov_vars;
    mov_vars.push_back("DX");
    mov_vars.push_back("DY");
    mov_vars.push_back("DZ");
    mlSol.GetWriter()->SetMovingMesh(mov_vars);

    std::vector<std::string> print_vars;
    print_vars.push_back("All");


    mlSol.GetWriter()->SetDebugOutput(true);
    mlSol.GetWriter()->Write(outputFolder.str(), "biquadratic", print_vars, 0);
    
    std::vector < std::vector < std::vector < double > > > line(1);
    linea->GetLine(line[0]);
    PrintLine(outputFolder.str(), "line", line, 0);
    
   
    
    
    std::ostringstream filename;
    filename <<outputFolder.str()<< "/centerOfMass.txt";
    std::ofstream fout;
    fout.open( filename.str().c_str() );
    fout << "iteration, time, ";
    fout << "E" << simulation.first  << "Level"<<simulation.second << std::endl;
    
    CM[0][simulation] = line[0][0][1];  
    time[0][simulation] = 0.;  
    fout << "0, 0., "<< CM[0][simulation] << std::endl;
    
    
    double theta = 0;
    gravity[0] = 9.81 * sin(theta);
    gravity[1] = -9.81 * cos(theta);

    system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
    
    for(unsigned time_step = 1; time_step <= n_timesteps; time_step++) {

      std::vector< double > xMin(3);
      std::vector< double > xMax(3);

      linea->GetExtrema(xMin, xMax);

      DT = GetTimeStep(xMin[1]);
      
      //std::cout << "ymin = " << xMin[1] << " DT = " << DT << std::endl;
    
      system.CopySolutionToOldSolution();

      system.MGsolve();
      
      GridToParticlesProjection(ml_prob, *linea);
       
      if(xMin[1] < -1.336){
        system2.MGsolve();
      }
        
      // ******* Print solution *******
      mlSol.GetWriter()->Write(outputFolder.str(), "biquadratic", print_vars, time_step);

//       // ******* Print solution *******
//       mlSol.GetWriter()->Write(outputFolder.str(), "biquadratic", print_vars, time_step);
// 
//       GridToParticlesProjection(ml_prob, *linea);

      linea->GetLine(line[0]);
      PrintLine(outputFolder.str(),"line", line, time_step);

      CM[time_step][simulation] = line[0][0][1];  
        
      time[time_step][simulation] = system.GetTime();
      
      std::cout << "it = " << time_step <<" time = "<< time[time_step][simulation] <<" DT = " << DT << " ymin = " << xMin[1] << std::endl;
      
      fout << time_step <<", "<< time[time_step][simulation] <<", ";
      fout << CM[time_step][simulation] <<std::endl;
      
      if( time[time_step][simulation] > 3.5) break;
    }
    
    fout.close();
    delete linea;
     
  }
  
  
//   std::ostringstream filename;
//   filename << "./centerOfMass"<< material <<".txt";
//   std::ofstream fout;
//   fout.open( filename.str().c_str() );
//   fout << "iteration, time, ";
//   for(unsigned nl = 1; nl < 4; nl++) {
//     std::pair <std::string, unsigned > simulation = std::make_pair (material, nl);
//      fout << "E" << simulation.first  << "Level"<<simulation.second<<", "; 
//   }
//   
//   fout << std::endl;
//   for(unsigned it = 0;it < CM.size();it++){
//     
//     fout << it <<", "<< time[it][simulation] << ", ";
//     for(unsigned nl = 1; nl < 4; nl++) {
//       std::pair <std::string, unsigned > simulation = std::make_pair (material, nl);
//       fout << CM[it][simulation] <<", ";
//     }
//     fout << std::endl;
//   }
//   
//   fout.close();
  
  
  
  return 0;

} //end main



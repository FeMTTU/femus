#include "MultiLevelProblem.hpp"
#include "TransientSystem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"
#include "FemusInit.hpp"
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "FElemTypeEnum.hpp"
#include "Files.hpp"
#include "../include/FSIassembly.hpp"

using std::cout;
using std::endl;
using namespace femus;

double SetVariableTimeStep(const double time);

bool SetBoundaryCondition(const std::vector < double >& x, const char name[],
			  double &value, const int FaceName, const double = 0.);

bool SetRefinementFlag(const std::vector < double >& x, const int &ElemGroupNumber,const int &level);

//------------------------------------------------------------------------------------------------------------------

int main(int argc,char **args) {
  
  /// Init Petsc-MPI communicator
  FemusInit mpinit(argc,args,MPI_COMM_WORLD);

  Files files; 
        files.CheckIODirectories();
        files.RedirectCout();
	
  unsigned short nm,nr;
  std::cout<<"#MULTIGRID levels? (>=1) \n";
  //std::cin>>nm;
  nm=4;

  std::cout<<"#MAX_REFINEMENT levels? (>=0) \n";
  //std::cin>>nr;
  nr=0;
  int tmp=nm;
  nm+=nr;
  nr=tmp;

  const std::string infile = "./input/fsifirst.neu";
  
  double Lref = 1.;
  double Uref = 1.;
  double rhof = 1000.;
  double muf = 1.;
  double rhos = 1000;
  double ni = 0.4;
  double E = 5600000;
  
  MultiLevelMesh ml_msh(nm,nr,infile.c_str(),"fifth",Lref,SetRefinementFlag);
  
  MultiLevelSolution ml_sol(&ml_msh);
  
  //Start System Variables
  ml_sol.AddSolution("DX",LAGRANGE,SECOND,2);
  ml_sol.AddSolution("DY",LAGRANGE,SECOND,2);
//   ml_sol.AssociatePropertyToSolution("DX","Displacement"); // Add this line
//   ml_sol.AssociatePropertyToSolution("DY","Displacement"); // Add this line 
  ml_sol.AddSolution("U",LAGRANGE,SECOND,2);
  ml_sol.AddSolution("V",LAGRANGE,SECOND,2);
  
  //   ml_sol.PairSolution("U","DX"); // Add this line
  //   ml_sol.PairSolution("V","DY"); // Add this line 
  
  ml_sol.AddSolution("AX",LAGRANGE,SECOND,1,0);
  ml_sol.AddSolution("AY",LAGRANGE,SECOND,1,0);
  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("P",DISCONTINOUS_POLYNOMIAL,FIRST,1);
  ml_sol.AssociatePropertyToSolution("P","Pressure"); // Add this line

  //Initialize (update Init(...) function)
  ml_sol.Initialize("All");

  //Set Boundary (update Dirichlet(...) function)
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  ml_sol.GenerateBdc("DX","Steady");
  ml_sol.GenerateBdc("DY","Steady");
  ml_sol.GenerateBdc("U","Time_dependent");
  ml_sol.GenerateBdc("V","Steady");
  ml_sol.GenerateBdc("AX","Steady");
  ml_sol.GenerateBdc("AY","Steady");
  ml_sol.GenerateBdc("P","Steady");

  
  MultiLevelProblem ml_prob(&ml_sol);
  

  Parameter par(Lref,Uref);
  
  // Generate Solid Object
  Solid solid(par,E,ni,rhos,"Neo-Hookean");
  cout << "Solid properties: " << endl;
  cout << solid << endl;
  
  // Generate Fluid Object
  Fluid fluid(par,muf,rhof,"Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  // Add fluid object
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
  
  // Add Solid Object
  ml_prob.parameters.set<Solid>("Solid") = solid;

  ml_msh.MarkStructureNode();
   
  //create systems
  // add the system FSI to the MultiLevel problem
  TransientMonolithicFSINonlinearImplicitSystem & system = ml_prob.add_system<TransientMonolithicFSINonlinearImplicitSystem> ("Fluid-Structure-Interaction");
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");
  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");
  system.AddSolutionToSystemPDE("P");
  
  // init all the systems
  system.init();
   
  // System Fluid-Structure-Interaction
  system.SetAssembleFunction(AssembleMatrixResFSI);  
  system.SetMaxNumberOfLinearIterations(1);
  system.SetLinearConvergenceTolerance(1.e-8);  
  system.SetMgType(V_CYCLE);
  system.SetMaxNumberOfNonLinearIterations(4);
  system.SetNonLinearConvergenceTolerance(1.e-5);
  system.SetDirichletBCsHandling(PENALTY);
  
  //system.SetDirichletBCsHandling(ELIMINATION);
  
  // time loop parameter
  system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
  const unsigned int n_timesteps = 5;
  const unsigned int write_interval = 1;
  
  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  VTKWriter vtkio(&ml_sol);
  vtkio.SetMovingMesh(mov_vars);
  
  for (unsigned time_step = 0; time_step < n_timesteps; time_step++) {
   
    // Solving Fluid-Structure-Interaction system
    std::cout << std::endl;
    std::cout << " *********** Fluid-Structure-Interaction ************  " << std::endl;
    system.solve();
   
    //The update of the acceleration must be done before the update of the other variables
    system.NewmarkAccUpdate();
    
    //update Solution
    system.UpdateSolution();

    // print solution
    if ( !(time_step%write_interval) ) {
        
      //print solution 
      std::vector<std::string> print_vars;
      print_vars.push_back("DX");
      print_vars.push_back("DY");
      print_vars.push_back("U");
      print_vars.push_back("V");
      print_vars.push_back("P");
      
//       ml_prob.printsol_vtu_inline("biquadratic",print_vars,time_step);
      vtkio.write(files.GetOutputPath(),"biquadratic",print_vars,time_step);
    }
  
  } //end loop timestep
  

  // Destroy all the new systems
  ml_prob.clear();
   
  return 0;
}


bool SetRefinementFlag(const std::vector < double >& x, const int &elemgroupnumber,const int &level) {
  bool refine=0;

  //refinemenet based on elemen group number
  if (elemgroupnumber==5) refine=1;
  if (elemgroupnumber==6) refine=1;
  if (elemgroupnumber==7 && level<5) refine=1;

  return refine;

}

//---------------------------------------------------------------------------------------------------------------------

double SetVariableTimeStep(const double time) {
 if(time < 4.) {
   return 0.01;
 } 
 else {
   return 0.001; 
 }
}

//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryCondition(const std::vector < double >& x,const char name[], double &value, const int facename, const double time) {
  bool test=1; //dirichlet
  value=0.;
  if(!strcmp(name,"U")) {
    if(1==facename){   //inflow
      test=1;
      double um = 2.0;
      if(time < 2.0) {
        value=1.5*um*4.0/0.1681*x[1]*(0.41-x[1])*0.5*(1. - cos(0.5*3.141592653589793*time));
       }
      else {
        value=1.5*um*4.0/0.1681*x[1]*(0.41-x[1]);
      }
    }  
    else if(2==facename ){  //outflow
     test=0;
 //    test=1;
     value=0.;
    }
    else if(3==facename ){  // no-slip fluid wall
      test=1;
      value=0.;	
    }
    else if(4==facename ){  // no-slip solid wall
      test=1;
      value=0.;
    }
  }  
  else if(!strcmp(name,"V")){
    if(1==facename){            //inflow
      test=1;
      value=0.;
    }  
    else if(2==facename ){      //outflow
     test=0;
 //    test=1;
     value=0.;
    }
    else if(3==facename ){      // no-slip fluid wall
      test=1;
      value=0;
    }
    else if(4==facename ){      // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"W")){
    if(1==facename){
      test=1;
      value=0.;
    }  
    else if(2==facename ){  
      test=1;
      value=0.;
    }
    else if(3==facename ){  
      test=1;
      value=0.;
    }
    else if(4==facename ){  
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"P")){
    if(1==facename){
      test=0;
      value=0.;
    }  
    else if(2==facename ){  
      test=0;
      value=0.;
    }
    else if(3==facename ){  
      test=0;
      value=0.;
    }
    else if(4==facename ){  
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"DX")){
    if(1==facename){         //inflow
      test=1;
      value=0.;
    }  
    else if(2==facename ){   //outflow
     test=1;
     value=0.;
    }
    else if(3==facename ){   // no-slip fluid wall
      test=0; //0
      value=0.;	
    }
    else if(4==facename ){   // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"DY")){
    if(1==facename){         //inflow
      test=0; // 0
      value=0.;
    }  
    else if(2==facename ){   //outflow
     test=0; // 0
     value=0.;
    }
    else if(3==facename ){   // no-slip fluid wall
      test=1;
      value=0.;	
    }
    else if(4==facename ){   // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"DZ")){
    if(1==facename){         //inflow
      test=1;
      value=0.;
    }  
    else if(2==facename ){   //outflow
     test=1;
     value=0.;
    }
    else if(3==facename ){   // no-slip fluid wall
      test=1;
      value=0.;	
    }
    else if(4==facename ){   // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if (!strcmp(name,"AX")) {
    test=0;
    value=0;
  }
  else if (!strcmp(name,"AY")) {
    test=0;
    value=0;
  }
  else if (!strcmp(name,"AZ")) {
    test=0;
    value=0;
  }
  return test;
}


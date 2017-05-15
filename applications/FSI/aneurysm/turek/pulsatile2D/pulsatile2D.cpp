#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"
#include "FemusInit.hpp"
#include "SparseMatrix.hpp"
#include "FElemTypeEnum.hpp"
#include "Files.hpp"
#include "MonolithicFSINonLinearImplicitSystem.hpp"
#include "TransientSystem.hpp"
#include "VTKWriter.hpp"
#include "../../include/FSITimeDependentAssemblySupg.hpp"
#include <cmath>
double scale = 1000.;

using namespace std;
using namespace femus;

double SetVariableTimeStep(const double time);

bool SetBoundaryConditionTurek2D(const std::vector < double >& x, const char name[],
                                 double &value, const int facename, const double time);

bool SetBoundaryConditionThrombus2D(const std::vector < double >& x, const char name[],
                                    double &value, const int facename, const double time);

bool SetBoundaryConditionAorticBifurcation(const std::vector < double >& x, const char name[],
    double &value, const int facename, const double time);

bool SetBoundaryConditionThrombus2DPorous(const std::vector < double >& x, const char name[],
    double &value, const int facename, const double time);

bool SetBoundaryConditionVeinValve(const std::vector < double >& x, const char name[],
                                   double &value, const int facename, const double time);

void GetSolutionNorm(MultiLevelSolution& mlSol, const unsigned & group, std::vector <double> &data);
//------------------------------------------------------------------------------------------------------------------


const double leaflet[129][2] = {
  {   -0.00402782    ,  0.0601179     },  {   -0.00391292    ,	0.0601527     },  {   -0.00379967    ,	0.0601912     },  {   -0.00368805    ,	0.0602334     },
  {   -0.00357807    ,	0.0602794     },  {   -0.00346973    ,	0.0603291     },  {   -0.00336303    ,	0.0603826     },  {   -0.00325798    ,	0.0604397     },
  {   -0.00315456    ,	0.0605006     },  {   -0.00305329    ,	0.0605649     },  {   -0.00295415    ,	0.0606319     },  {   -0.00285713    ,	0.0607016     },
  {   -0.00276223    ,	0.060774      },  {   -0.00266945    ,	0.060849      },  {   -0.00257879    ,	0.0609268     },  {   -0.00249026    ,	0.0610072     },
  {   -0.00240385    ,	0.0610904     },  {   -0.00231958    ,	0.0611756     },  {   -0.00223735    ,	0.0612626     },  {   -0.00215715    ,	0.0613513     },
  {   -0.00207899    ,	0.0614417     },  {   -0.00200286    ,	0.0615339     },  {   -0.00192876    ,	0.0616277     },  {   -0.0018567     ,  0.0617232     },
  {   -0.00178666    ,	0.0618205     },  {   -0.00171853    ,	0.061919      },  {   -0.00165215    ,	0.0620186     },  {   -0.00158754    ,	0.0621193     },
  {   -0.00152469    ,	0.062221      },  {   -0.00146361    ,	0.0623238     },  {   -0.00140429    ,	0.0624276     },  {   -0.00134673    ,	0.0625325     },
  {   -0.00129093    ,	0.0626385     },  {   -0.00123672    ,	0.0627453     },  {   -0.00118402    ,	0.0628528     },  {   -0.00113285    ,	0.062961      },
  {   -0.00108319    ,	0.0630698     },  {   -0.00103505    ,  0.0631794     },  {   -0.000988434   ,	0.0632896     },  {   -0.000943336   ,	0.0634005     },
  {   -0.000899757   ,	0.063512      },  {   -0.000857621   ,	0.0636239     },  {   -0.000816835   ,	0.0637363     },  {   -0.000777398   ,	0.0638492     },
  {   -0.000739311   ,	0.0639624     },  {   -0.000702573   ,	0.0640761     },  {   -0.000667185   ,	0.0641903     },  {   -0.000633146   ,	0.0643049     },
  {   -0.000600456   ,	0.0644199     },  {   -0.000568998   ,	0.0645354     },  {   -0.000538812   ,	0.0646513     },  {   -0.000509899   ,	0.0647674     },
  {   -0.000482258   ,	0.0648837     },  {   -0.00045589    ,  0.0650004     },  {   -0.000430794   ,	0.0651174     },  {   -0.000406971   ,	0.0652347     },
  {   -0.00038442    ,  0.0653523     },  {   -0.000363115   ,	0.06547       },  {   -0.000343096   ,	0.0655879     },  {   -0.000324363   ,	0.0657059     },
  {   -0.000306917   ,	0.0658242     },  {   -0.000290758   ,	0.0659426     },  {   -0.000275886   ,	0.0660612     },  {   -0.000262299   ,	0.06618       },
  {   -0.00025	     ,  0.066299      },  {   -0.000239534   ,	0.0664144     },  {   -0.000230166   ,	0.0665298     },  {   -0.000221896   ,	0.0666454     },
  {   -0.000214725   ,	0.0667609     },  {   -0.000208652   ,	0.0668766     },  {   -0.000203678   ,	0.0669923     },  {   -0.000199802   ,	0.067108      },
  {   -0.000197024   ,	0.0672238     },  {   -0.000195205   ,	0.0673395     },  {   -0.000194055   ,	0.0674552     },  {   -0.000193574   ,	0.0675708     },
  {   -0.000193763   ,	0.0676864     },  {   -0.000194621   ,	0.0678021     },  {   -0.000196148   ,	0.0679177     },  {   -0.000198345   ,	0.0680334     },
  {   -0.000201211   ,	0.068149      },  {   -0.000204526   ,	0.0682648     },  {   -0.000207968   ,	0.0683806     },  {   -0.000211535   ,	0.0684963     },
  {   -0.000215229   ,	0.0686121     },  {   -0.000219049   ,	0.0687279     },  {   -0.000222994   ,	0.0688437     },  {   -0.000227066   ,	0.0689594     },
  {   -0.000231264   ,	0.0690752     },  {   -0.000235294   ,	0.0691908     },  {   -0.000238841   ,	0.0693063     },  {   -0.000241906   ,	0.0694219     },
  {   -0.000244489   ,	0.0695375     },  {   -0.00024659    ,  0.0696531     },  {   -0.000248209   ,	0.0697687     },  {   -0.000249345   ,	0.0698844     },
  {   -0.00025       ,  0.07          },  {   -0.000255613   ,	0.0701159     },  {   -0.000261276   ,	0.0702317     },  {   -0.00026699    ,  0.0703476     },
  {   -0.000272754   ,	0.0704635     },  {   -0.000278569   ,	0.0705793     },  {   -0.000284435   ,	0.0706952     },  {   -0.00029035    ,  0.070811      },
  {   -0.000296317   ,	0.0709269     },  {   -0.000302377   ,	0.0710426     },  {   -0.000308591   ,	0.0711583     },  {   -0.000314961   ,	0.0712741     },
  {   -0.000321485   ,	0.0713898     },  {   -0.000328164   ,	0.0715055     },  {   -0.000334998   ,	0.0716212     },  {   -0.000341987   ,	0.0717368     },
  {   -0.000349131   ,	0.0718525     },  {   -0.000356484   ,	0.0719682     },  {   -0.000364109   ,	0.072084      },  {   -0.000372005   ,	0.0721997     },
  {   -0.000380173   ,	0.0723154     },  {   -0.000388611   ,	0.072431      },  {   -0.000397321   ,	0.0725467     },  {   -0.000406302   ,	0.0726623     },
  {   -0.000415554   ,	0.0727779     },  {   -0.000425128   ,	0.0728934     },  {   -0.00043511    ,  0.0730089     },  {   -0.0004455     ,  0.0731244     },
  {   -0.000456299   ,	0.0732398     },  {   -0.000467505   ,	0.0733551     },  {   -0.000479121   ,	0.0734705     },  {   -0.000491144   ,	0.0735857     },
  {   -0.000503576   ,	0.073701      }
};

bool MoveMesh(const std::vector < double >& x, std::vector < double >& dX) {
  bool movedMesh = false;
  double epsilon = 1.0e-5; // maximum leaflet colsing
  double  delta = 0.000071; //leaflet thicknes;
  unsigned N = 128;
  if(x[1] >= leaflet[0][1] && x[1] <= leaflet[N][1]) {
    unsigned i0 = 0;
    unsigned i1 = N;
    while(i1 - i0 > 1) {
      unsigned i2 = i0 + (i1 - i0) / 2;
      if(leaflet[i2][1] > x[1]) {
        i1 = i2;
      }
      else {
        i0 = i2;
      }
    }
    double s = (x[1] - leaflet[i0][1]) / (leaflet[i1][1] - leaflet[i0][1]);
    double xl = leaflet[i0][0] * (1. - s) + leaflet[i1][0] * s;

    if(x[0] - xl >= 0) {  //on the right of the lealflet
      if(x[0] + dX[0] > -epsilon * x[0] / xl) {
        dX[0] = -epsilon * x[0] / xl - x[0] ;
        movedMesh = true;
      }
    }
    else if(x[0] - xl >= -delta) {  //inside the lealflet
      if(x[0] + dX[0] > -epsilon + (x[0] - xl)) {
        dX[0] = -epsilon - (xl - x[0]) - x[0];
        movedMesh = true;
      }
    }
    else { //on the right of the lealflet
      if(x[0] + dX[0] > -epsilon - delta - 3.* epsilon * (xl - delta - x[0]) / 0.001) {
        dX[0] = -epsilon - delta - 3.* epsilon * (xl - delta - x[0]) / 0.001 - x[0];
      }
    }
    std::cout << i0 << " " << i1 << " " << xl << " " << s << std::endl;
    std::cout << x[0] + dX[0] << std::endl;
    
  }
  return movedMesh;
}

int main(int argc, char **args)
{

  // ******* Init Petsc-MPI communicator *******
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  std::vector < double > x(2);
  x[0] = -0.000425128 - 0.000071 - 0.001;
  x[1] =  0.0728934;
  std::vector < double > dX(2);
  dX[0] = 1;
  dX[1] = 0.;
  MoveMesh(x, dX);
 


  unsigned simulation = 0;

  return 1;


  if(argc >= 2) {
    if(!strcmp("0", args[1])) {       /** FSI Turek2D no stent */
      simulation = 0;
    }
    else if(!strcmp("1", args[1])) {       /** FSI Turek porous */
      simulation = 1;
    }
    else if(!strcmp("2", args[1])) {       /** FSI Turek stents 60 micron */
      simulation = 2;
    }
    else if(!strcmp("3", args[1])) {       /** FSI Turek 11 stents 60 micron */
      simulation = 3;
    }
    else if(!strcmp("4", args[1])) {       /** FSI AAA thrombus 2D */
      simulation = 4;
    }
    else if(!strcmp("5", args[1])) {       /** FSI Aortic Bifurcation */
      simulation = 5;
    }
    else if(!strcmp("6", args[1])) {       /** FSI AAA thrombus 2D Porous */
      simulation = 6;
    }
    else if(!strcmp("7", args[1])) {       /** FSI Vein Valve */
      simulation = 7;
    }
  }

  //Files files;
  //files.CheckIODirectories();
  //files.RedirectCout();

  // ******* Extract the problem dimension and simulation identifier based on the inline input *******


  // ******* Extract the mesh.neu file name based on the simulation identifier *******
//   std::string infile = "./input/aneurysm_Sara_5.neu";
  //std::string infile = "./input/Turek_porous_60micron.neu";
  //std::string infile = "./input/Turek_stents_60micron.neu";
  //std::string infile = "./input/Turek_11stents_60micron.neu";
  std::string infile;
  if(simulation == 0) {
    infile = "./input/Turek.neu";
  }
  else if(simulation == 1) {
    infile = "./input/Turek_porous_60micron.neu";
  }
  else if(simulation == 2) {
    infile = "./input/Turek_stents_60micron.neu";
  }
  else if(simulation == 3) {
    infile = "./input/Turek_11stents_60micron.neu";
  }
  else if(simulation == 4) {
    infile = "./input/AAA_thrombus_2D.neu";
  }
  else if(simulation == 5) {
    infile = "./input/aortic_bifurcation.neu";
  }
  else if(simulation == 6) {
    infile = "./input/AAA_thrombus_2D_porous.neu";
  }
  else if(simulation == 7) {
    //infile = "./input/vein_valve_closed.neu";
    //infile = "./input/vein_valve_thiner.neu";
    infile = "./input/vein_valve_modifiedFluid.neu";
  }

  // ******* Set physics parameters *******
  double Lref, Uref, rhof, muf, rhos, ni, E, E1;

  Lref = 1.;
  Uref = 1.;


  if(simulation == 7) {
    rhof = 1060.;
    muf = 2.2 * 1.0e-3;
    rhos = 960;
    ni = 0.5;
    //E = 3.3 * 1.0e6; //vein young modulus
    E = 4.3874951 * 1.0e12;
    E1 = 15 * 1.0e5; //leaflet young modulus
  }
  else {
    rhof = 1035.;
    muf = 3.5 * 1.0e-3; //wrong=3.38*1.0e-4*rhof, note:3.38*1.0e-6*rhof=3.5*1.0e-3
    rhos = 1120;
    ni = 0.5;
    //E = 100 * 1.e6; //simulation=5
    E = 1. * 1.0e6; //simulation=0,1,2,3
    E1 = 50000;
  }

  Parameter par(Lref, Uref);

  // Generate Solid Object
  Solid solid;
  solid = Solid(par, E, ni, rhos, "Mooney-Rivlin");

  Solid solid1;
  solid1 = Solid(par, E1, ni, rhos, "Mooney-Rivlin");

  cout << "Solid properties: " << endl;
  cout << solid << endl;

  // Generate Fluid Object
  Fluid fluid(par, muf, rhof, "Newtonian");
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;

  // ******* Init multilevel mesh from mesh.neu file *******
  unsigned short numberOfUniformRefinedMeshes, numberOfAMRLevels;

  numberOfUniformRefinedMeshes = 3;
  numberOfAMRLevels = 0;

  std::cout << 0 << std::endl;

  MultiLevelMesh ml_msh(numberOfUniformRefinedMeshes + numberOfAMRLevels, numberOfUniformRefinedMeshes,
                        infile.c_str(), "fifth", Lref, NULL);

  //ml_msh.EraseCoarseLevels(numberOfUniformRefinedMeshes - 2);

  ml_msh.PrintInfo();

  // mark Solid nodes

  //ml_msh.MarkStructureNode();

  // ******* Init multilevel solution ******
  MultiLevelSolution ml_sol(&ml_msh);

  // ******* Add solution variables to multilevel solution and pair them *******
  ml_sol.AddSolution("DX", LAGRANGE, SECOND, 2);
  ml_sol.AddSolution("DY", LAGRANGE, SECOND, 2);

  ml_sol.AddSolution("U", LAGRANGE, SECOND, 2);
  ml_sol.AddSolution("V", LAGRANGE, SECOND, 2);

  // Pair each velocity variable with the corresponding displacement variable
  ml_sol.PairSolution("U", "DX");    // Add this line
  ml_sol.PairSolution("V", "DY");    // Add this line

  // Since the Pressure is a Lagrange multiplier it is used as an implicit variable
  ml_sol.AddSolution("P", DISCONTINOUS_POLYNOMIAL, FIRST, 2);
  ml_sol.AssociatePropertyToSolution("P", "Pressure", false);    // Add this line

  ml_sol.AddSolution("lmbd", DISCONTINOUS_POLYNOMIAL, ZERO, 0, false);

  // ******* Initialize solution *******
  ml_sol.Initialize("All");

  if(simulation == 0 || simulation == 1 || simulation == 2 || simulation == 3) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionTurek2D);
  }
  else if(simulation == 4 || simulation == 6) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionThrombus2D);
  }
  else if(simulation == 5) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionAorticBifurcation);
  }
  else if(simulation == 7) {
    ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionVeinValve);
  }

  // ******* Set boundary conditions *******
  ml_sol.GenerateBdc("DX", "Steady");
  ml_sol.GenerateBdc("DY", "Steady");

  if(simulation == 4 || simulation == 5 || simulation == 6) {
    ml_sol.GenerateBdc("U", "Steady");
    ml_sol.GenerateBdc("V", "Time_dependent");
  }
  else if(simulation == 7) {
    ml_sol.GenerateBdc("U", "Steady");
    ml_sol.GenerateBdc("V", "Steady");
  }
  else {
    ml_sol.GenerateBdc("U", "Time_dependent");
    ml_sol.GenerateBdc("V", "Steady");
  }
  if(simulation == 7) {
    ml_sol.GenerateBdc("P", "Steady");
  }
  else {
    ml_sol.GenerateBdc("P", "Time_dependent");
  }

//   for(unsigned level = 0; level < numberOfUniformRefinedMeshes; level++ ){
//     SetLambda(ml_sol, level , SECOND, ELASTICITY);
//   }

  // ******* Define the FSI Multilevel Problem *******

  MultiLevelProblem ml_prob(&ml_sol);
  // Add fluid object
  ml_prob.parameters.set<Fluid> ("Fluid") = fluid;
  // Add Solid Object
  ml_prob.parameters.set<Solid> ("Solid") = solid;
  ml_prob.parameters.set<Solid> ("Solid1") = solid1;

  // ******* Add FSI system to the MultiLevel problem *******
  TransientMonolithicFSINonlinearImplicitSystem & system = ml_prob.add_system<TransientMonolithicFSINonlinearImplicitSystem> ("Fluid-Structure-Interaction");
  system.AddSolutionToSystemPDE("DX");
  system.AddSolutionToSystemPDE("DY");

  system.AddSolutionToSystemPDE("U");
  system.AddSolutionToSystemPDE("V");

  system.AddSolutionToSystemPDE("P");

  // ******* System Fluid-Structure-Interaction Assembly *******
  system.SetAssembleFunction(FSITimeDependentAssemblySupg);

  // ******* set MG-Solver *******
  system.SetMgType(F_CYCLE);

  system.SetNonLinearConvergenceTolerance(1.e-7);
  //system.SetResidualUpdateConvergenceTolerance ( 1.e-15 );
  system.SetMaxNumberOfNonLinearIterations(10);
  //system.SetMaxNumberOfResidualUpdatesForNonlinearIteration ( 4 );

  system.SetMaxNumberOfLinearIterations(6);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-13);

  system.SetNumberPreSmoothingStep(0);
  system.SetNumberPostSmoothingStep(2);

  // ******* Set Preconditioner *******

  system.SetMgSmoother(ASM_SMOOTHER);

  system.init();

  // ******* Set Smoother *******
  system.SetSolverFineGrids(RICHARDSON);
  //system.SetSolverFineGrids(GMRES);

  system.SetPreconditionerFineGrids(ILU_PRECOND);

  system.SetTolerances(1.e-12, 1.e-20, 1.e+50, 20, 10);

  // ******* Add variables to be solved *******
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  // ******* Set the last (1) variables in system (i.e. P) to be a schur variable *******
  system.SetNumberOfSchurVariables(1);

  // ******* Set block size for the ASM smoothers *******
  system.SetElementBlockNumber(2);

  // ******* Print solution *******
  ml_sol.SetWriter(VTK);


  std::vector<std::string> mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  //mov_vars.push_back("DZ");
  ml_sol.GetWriter()->SetMovingMesh(mov_vars);

  std::vector<std::string> print_vars;
  print_vars.push_back("All");

  ml_sol.GetWriter()->SetDebugOutput(true);
  ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, 0);


  // ******* Solve *******
  std::cout << std::endl;
  std::cout << " *********** Fluid-Structure-Interaction ************  " << std::endl;

  // time loop parameter
  system.AttachGetTimeIntervalFunction(SetVariableTimeStep);
  const unsigned int n_timesteps = 1024;


  std::vector < std::vector <double> > data(n_timesteps);

  for(unsigned time_step = 0; time_step < n_timesteps; time_step++) {
    for(unsigned level = 0; level < numberOfUniformRefinedMeshes; level++) {
      SetLambda(ml_sol, level , SECOND, ELASTICITY);
    }
    data[time_step].resize(5);
    if(time_step > 0)
      system.SetMgType(V_CYCLE);
    system.CopySolutionToOldSolution();

    system.MGsolve();
    //data[time_step][0] = time_step / 16.;
    //data[time_step][0] = time_step / 20.;
    data[time_step][0] = time_step / 32.;
    //data[time_step][0] = time_step / ( 64 * 1.4 );
    if(simulation == 0 || simulation == 1 || simulation == 2 || simulation == 3) {
      GetSolutionNorm(ml_sol, 9, data[time_step]);
    }
    else if(simulation == 4) {    //AAA_thrombus, 15=thrombus
      GetSolutionNorm(ml_sol, 7, data[time_step]);
    }
    else if(simulation == 6) {    //AAA_thrombus_porous, 15=thrombus
      GetSolutionNorm(ml_sol, 7, data[time_step]);
    }
    ml_sol.GetWriter()->Write(DEFAULT_OUTPUTDIR, "biquadratic", print_vars, time_step + 1);
  }


  int  iproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  if(iproc == 0) {
    std::ofstream outf;
    if(simulation == 0) {
      outf.open("DataPrint_Turek.txt");
    }
    else if(simulation == 1) {
      outf.open("DataPrint_TurekPorous.txt");
    }
    else if(simulation == 2) {
      outf.open("DataPrint_TurekStents.txt");
    }
    else if(simulation == 3) {
      outf.open("DataPrint_Turek11Stents.txt");
    }
    else if(simulation == 4) {
      outf.open("DataPrint_AAA_thrombus_2D.txt");
    }
    else if(simulation == 6) {
      outf.open("DataPrint_AAA_thrombus_2D_porous.txt");
    }


    if(!outf) {
      std::cout << "Error in opening file DataPrint.txt";
      return 1;
    }
    for(unsigned k = 0; k < n_timesteps; k++) {
      outf << data[k][0] << "\t" << data[k][1] << "\t" << data[k][2] << "\t" << data[k][3] << "\t" << data[k][4] << std::endl;
    }
    outf.close();
  }



  // ******* Clear all systems *******
  ml_prob.clear();
  return 0;
}

double SetVariableTimeStep(const double time)
{
  //double dt = 1. / ( 64 * 1.4 );
  double dt = 1. / 32;
  double shiftedTime = time - floor(time);
  if(time > 1 && shiftedTime >= 0.125 && shiftedTime < 0.25) {
//     dt = 1./256;
    dt = 1. / 32;
  }
  std::cout << " Shifted Time = " << shiftedTime << " dt = " << dt << std::endl;

//   double PI = acos(-1.);
//   double dt = 1./32.*(1.25 + 0.75 * sin(2.*PI*time - 3./4. * PI));
//   std::cout << "time = " << time <<  " dt * 32 = " << dt * 32<<std::endl;
  //double dt = 1. / 20;
  //double dt = 1./16.;
  //double dt = 1.0e6;

//   if( turek_FSI == 2 ){
//     if ( time < 9 ) dt = 0.05;
//     else dt = 0.025;
//   }
//   else if ( turek_FSI == 3 ){
//     //if	    ( time < 5. ) dt = 0.1;
//     //else
//     if ( time < 6. ) dt = 0.01;
//     else             dt = 0.01;
//   }
//   else if ( simulation == 3 ) dt=0.001;
//   else if ( simulation == 4 ) dt=0.1;
//   else if ( simulation == 5 ) dt=0.1;
//   else if ( simulation == 6 ) dt=0.1;
//   else if ( simulation == 7 ) dt=0.001;
//   else{
//     std::cout << "Warning this simulation case has not been considered yet for the time dependent case"<<std::endl;
//     abort();
//   }
  return dt;
}


//---------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionTurek2D(const std::vector < double >& x, const char name[], double &value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;


//   std::ifstream inf;
//   inf.open ( "./input/womersleyProfile_velMax65cms.txt" );
//   if ( !inf ) {
//     std::cout << "velocity file ./input/womersleyProfile_velMax65cms.txt can not be opened\n";
//     exit ( 0 );
//   }
//
//   std::vector<double> vel ( 64 );
//
//   for ( unsigned i = 0; i < 64; i++ ) {
//     inf >> vel[i];
//   }
//   inf.close();
//
//   double period = 1. / 1.4;
//   double dt = period / 64;
//
//   double time1 = time - floor ( time / period ) * period;
//
//   unsigned j = static_cast < unsigned > ( floor ( time1 / dt ) );

  //git pstd::cout<< name << " " << time <<" "<< j <<" "<<  vel[j] << std::endl;

  double PI = acos(-1.);
  double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;
  //double ramp = ( time < period ) ? sin ( PI / 2 * time / period ) : 1.;

  if(!strcmp(name, "U")) {
    if(1 == facename) {
      //value = 0.05 * (x[1] * 1000 - 6) * ( x[1] * 1000 - 8); //inflow
      value = 0.05 * (x[1] * 1000 - 6) * (x[1] * 1000 - 8) * (1. + 0.75 * sin(2. * PI * time)) * ramp; //inflow
      //value = ( x[1] * 1000 - 6 ) * ( x[1] * 1000 - 8 ) * vel[j] * ramp; //inflow
    }
    else if(2 == facename || 5 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if(!strcmp(name, "V")) {
    //if ( 2 == facename ) {
    //if ( 5 == facename ) {
    if(2 == facename || 5 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if(!strcmp(name, "P")) {
    test = 0;
    value = 0.;
  }
  else if(!strcmp(name, "DX")) {
    //if(2 == facename || 4 == facename || 5 == facename || 6 == facename) {
    if(5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "DY")) {
    //if(1 == facename || 3 == facename || 5 == facename || 6 == facename) {
    if(5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;

}


bool SetBoundaryConditionThrombus2D(const std::vector < double >& x, const char name[], double &value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);

  double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;

  if(!strcmp(name, "V")) {
    if(1 == facename) {
      double r2 = (x[0] * 100.) * (x[0] * 100.);
      //value = -0.01/.9 * (.9 - r2); //inflow
      value = -0.05 / .81 * (.81 - r2) * (1. + 0.75 * sin(2.*PI * time)) * ramp;        //inflow
    }
    if(2 == facename || 5 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if(!strcmp(name, "U")) {
    if(2 == facename) {
      test = 0;
      value = (10000 + 2500 * sin(2 * PI * time)) * ramp;;
    }
    else if(5 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "P")) {
    test = 0;
    value = 0.;
  }
  else if(!strcmp(name, "DX")) {
    if(5 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "DY")) {
    if(5 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;
}


bool SetBoundaryConditionAorticBifurcation(const std::vector < double >& x, const char name[], double &value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);

  double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;

  if(!strcmp(name, "V")) {
    if(1 == facename) {
      double r2 = (x[0] * 100.) * (x[0] * 100.);
      //value = -0.01/.9 * (.9 - r2); //inflow
      value = -0.1 / .81 * (.81 - r2) * (1. + 0.25 * sin(2.*PI * time)) * ramp;        //inflow
    }
    if(2 == facename || 3 == facename || 7 == facename) {
      test = 0;
      value = 0.;
    }
  }
  else if(!strcmp(name, "U")) {
    if(2 == facename || 3 == facename || 7 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "P")) {
    test = 0;
    value = 0.;
//     if ( 2 == facename || 3 == facename ) {
//       //value = ( 12500 + 2500 * sin ( 2 * PI * time ) ) * ramp;
//       value = 50 * ramp;
//     }
  }
  else if(!strcmp(name, "DX")) {
    if(7 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "DY")) {
    if(7 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;
}

bool SetBoundaryConditionVeinValve(const std::vector < double >& x, const char name[], double &value, const int facename, const double time)
{
  bool test = 1; //dirichlet
  value = 0.;

  double PI = acos(-1.);
  double ramp = (time < 1) ? sin(PI / 2 * time) : 1.;

  if(!strcmp(name, "V")) {
    if(1 == facename || 2 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "P")) {
    test = 0;
    value = 0.;
    if(1 == facename) {
      //value = -1;
      //value = ( /*2.5*/ + 2.5 * sin ( 2 * PI * time ) ) * ramp;
      //value = ( 5 + 3 * sin ( 2 * PI * time ) ) * ramp; //+ 4.5
      //value = ( 6 + 3 * sin ( 2 * PI * time ) ) * ramp; //+ 4.5
      //value = ( 12 + 9 * sin ( 2 * PI * time ) ) * ramp; //runna
      //value = ( 24 + 21 * sin ( 2 * PI * time ) ) * ramp; //runna
      value = (0 + 2 * sin(2 * PI * time)) * ramp;      //+ 4.5
    }
    else if(2 == facename) {
      //value = 1;
      //value = ( /*2.5*/ - 2.5 * sin ( 2 * PI * time ) ) * ramp;
      //value = ( 4 - 1 * sin ( 2 * PI * time ) ) * ramp; //- 4.5
      //value = ( 5 - 3 * sin ( 2 * PI * time ) ) * ramp; //non runna
      value = (0 - 2 * sin(2 * PI * time)) * ramp;      //- 4.5
    }
  }
  else if(!strcmp(name, "DX")) {
    if(5 == facename) {
      test = 0;
      value = 0;
    }
  }
  else if(!strcmp(name, "DY")) {
    if(5 == facename || 6 == facename) {
      test = 0;
      value = 0;
    }
  }

  return test;

}

void GetSolutionNorm(MultiLevelSolution& mlSol, const unsigned & group, std::vector <double> &data)
{

  int  iproc, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  NumericVector* p2;
  NumericVector* v2;
  NumericVector* vol;
  NumericVector* vol0;
  p2 = NumericVector::build().release();
  v2 = NumericVector::build().release();
  vol = NumericVector::build().release();
  vol0 = NumericVector::build().release();

  if(nprocs == 1) {
    p2->init(nprocs, 1, false, SERIAL);
    v2->init(nprocs, 1, false, SERIAL);
    vol->init(nprocs, 1, false, SERIAL);
    vol0->init(nprocs, 1, false, SERIAL);
  }
  else {
    p2->init(nprocs, 1, false, PARALLEL);
    v2->init(nprocs, 1, false, PARALLEL);
    vol->init(nprocs, 1, false, PARALLEL);
    vol0->init(nprocs, 1, false, PARALLEL);
  }

  p2->zero();
  v2->zero();
  vol->zero();
  vol0->zero();

  unsigned level = mlSol._mlMesh->GetNumberOfLevels() - 1;

  Solution* solution  = mlSol.GetSolutionLevel(level);
  Mesh* msh = mlSol._mlMesh->GetLevel(level);


  const unsigned dim = msh->GetDimension();


  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

  vector< double > solP;
  vector< vector < double> >  solV(dim);
  vector< vector < double> > x0(dim);
  vector< vector < double> > x(dim);

  solP.reserve(max_size);
  for(unsigned d = 0; d < dim; d++) {
    solV[d].reserve(max_size);
    x0[d].reserve(max_size);
    x[d].reserve(max_size);
  }
  double weight;
  double weight0;

  vector <double> phiV;
  vector <double> gradphiV;
  vector <double> nablaphiV;

  double *phiP;

  phiV.reserve(max_size);
  gradphiV.reserve(max_size * dim);
  nablaphiV.reserve(max_size * (3 * (dim - 1) + !(dim - 1)));

  vector < unsigned > solVIndex(dim);
  solVIndex[0] = mlSol.GetIndex("U");    // get the position of "U" in the ml_sol object
  solVIndex[1] = mlSol.GetIndex("V");    // get the position of "V" in the ml_sol object
  if(dim == 3) solVIndex[2] = mlSol.GetIndex("W");       // get the position of "V" in the ml_sol object

  unsigned solVType = mlSol.GetSolutionType(solVIndex[0]);    // get the finite element type for "u"

  vector < unsigned > solDIndex(dim);
  solDIndex[0] = mlSol.GetIndex("DX");    // get the position of "U" in the ml_sol object
  solDIndex[1] = mlSol.GetIndex("DY");    // get the position of "V" in the ml_sol object
  if(dim == 3) solDIndex[2] = mlSol.GetIndex("DZ");       // get the position of "V" in the ml_sol object

  unsigned solDType = mlSol.GetSolutionType(solDIndex[0]);

  unsigned solPIndex;
  solPIndex = mlSol.GetIndex("P");
  unsigned solPType = mlSol.GetSolutionType(solPIndex);

  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
    if(msh->GetElementGroup(iel) == group) {
      short unsigned ielt = msh->GetElementType(iel);
      unsigned ndofV = msh->GetElementDofNumber(iel, solVType);
      unsigned ndofP = msh->GetElementDofNumber(iel, solPType);
      unsigned ndofD = msh->GetElementDofNumber(iel, solDType);
      // resize

      phiV.resize(ndofV);
      gradphiV.resize(ndofV * dim);
      nablaphiV.resize(ndofV * (3 * (dim - 1) + !(dim - 1)));

      solP.resize(ndofP);
      for(int d = 0; d < dim; d++) {
        solV[d].resize(ndofV);
        x0[d].resize(ndofD);
        x[d].resize(ndofD);
      }
      // get local to global mappings
      for(unsigned i = 0; i < ndofD; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solDType);
        for(unsigned d = 0; d < dim; d++) {
          x0[d][i] = (*msh->_topology->_Sol[d])(idof);

          x[d][i] = (*msh->_topology->_Sol[d])(idof) +
                    (*solution->_Sol[solDIndex[d]])(idof);
        }
      }

      for(unsigned i = 0; i < ndofV; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solVType);    // global to global mapping between solution node and solution dof
        for(unsigned  d = 0; d < dim; d++) {
          solV[d][i] = (*solution->_Sol[solVIndex[d]])(idof);      // global extraction and local storage for the solution
        }
      }



      for(unsigned i = 0; i < ndofP; i++) {
        unsigned idof = msh->GetSolutionDof(i, iel, solPType);
        solP[i] = (*solution->_Sol[solPIndex])(idof);
      }


      for(unsigned ig = 0; ig < mlSol._mlMesh->_finiteElement[ielt][solVType]->GetGaussPointNumber(); ig++) {
        // *** get Jacobian and test function and test function derivatives ***
        msh->_finiteElement[ielt][solVType]->Jacobian(x0, ig, weight0, phiV, gradphiV, nablaphiV);
        msh->_finiteElement[ielt][solVType]->Jacobian(x, ig, weight, phiV, gradphiV, nablaphiV);
        phiP = msh->_finiteElement[ielt][solPType]->GetPhi(ig);

        vol0->add(iproc, weight0);
        vol->add(iproc, weight);

        std::vector < double> SolV2(dim, 0.);
        for(unsigned i = 0; i < ndofV; i++) {
          for(unsigned d = 0; d < dim; d++) {
            SolV2[d] += solV[d][i] * phiV[i];
          }
        }

        double V2 = 0.;
        for(unsigned d = 0; d < dim; d++) {
          V2 += SolV2[d] * SolV2[d];
        }
        v2->add(iproc, V2 * weight);

        double P2 = 0;
        for(unsigned i = 0; i < ndofP; i++) {
          P2 += solP[i] * phiP[i];
        }
        P2 *= P2;
        p2->add(iproc, P2 * weight);
      }
    }
  }

  p2->close();
  v2->close();
  vol0->close();
  vol->close();

  double p2_l2 = p2->l1_norm();
  double v2_l2 = v2->l1_norm();
  double VOL0 = vol0->l1_norm();
  double VOL = vol->l1_norm();

  std::cout.precision(14);
  std::scientific;
  std::cout << " vol0 = " << VOL0 << std::endl;
  std::cout << " vol = " << VOL << std::endl;
  std::cout << " (vol-vol0)/vol0 = " << (VOL - VOL0) / VOL0 << std::endl;
  std::cout << " p_l2 norm / vol = " << sqrt(p2_l2 / VOL)  << std::endl;
  std::cout << " v_l2 norm / vol = " << sqrt(v2_l2 / VOL)  << std::endl;

  data[1] = (VOL - VOL0) / VOL0;
  data[2] = VOL;
  data[3] = sqrt(p2_l2 / VOL);
  data[4] = sqrt(v2_l2 / VOL);

  delete p2;
  delete v2;
  delete vol;

}


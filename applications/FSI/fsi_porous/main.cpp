// library includes
#include "FemusInit.hpp"
#include "ElemType.hpp"
#include "NumericVector.hpp"
#include "PetscVector.hpp"
#include "LinearEquationSolver.hpp" //linear_solver
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Parameter.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelMesh.hpp"
#include "FemusInputParser.hpp"

// application includes and prototypes
#include "main.hpp"
#include "MyMultigrid.hpp"

using namespace femus;

//===============================
int AssembleMatrixResP(MultiLevelProblem &mg2, unsigned level, /*const elem_type *type_elem[6][5], vector <vector <double> > &vt,*/ const unsigned &gridn, const unsigned &ipde, const bool &assemble_matrix);
int AssembleMatrixResD(MultiLevelProblem &mg2, unsigned level, /*const elem_type *type_elem[6][5], vector <vector <double> > &vt,*/ const unsigned &gridn, const unsigned &ipde, const bool &assemble_matrix);
int AssembleMatrixResVel(MultiLevelProblem &mg2, unsigned level, /*const elem_type *type_elem[6][5], vector <vector <double> > &vt,*/ const unsigned &gridn, const unsigned &ipde, const bool &assemble_matrix);
bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level);
bool BoundaryND(/*MultiLevelProblem& mg_in,*/const double &x, const double &y, const double &z,const char name[], double &value, const int FaceName,const double time);


FemusInputParser<double> * runtime_double; //per ora devo usarla cosi' global... dovrebbe appartenere a tutte le classi...
                                     //e' chiaro che lasciarla qui e' un problema perche' gli argomenti del costruttore...

int main(int argc,char **args) {

// // // // //   bool linear=1;
// // // // //   bool vanka=1;
// // // // //   if(argc == 2) {
// // // // //     if( strcmp("vanka",args[1])) vanka=0;
// // // // //   }
// // // // //   else {
// // // // //     cout << "No input arguments!" << endl;
// // // // //     exit(0);
// // // // //   }
// // // // //   
// // // // //   /// Init Petsc-MPI communicator
// // // // //   FemusInit mpinit(argc,args,MPI_COMM_WORLD);
// // // // //   
// // // // //  //READ STRING of TIME FROM SHELL;
// // // // //   std::string outfolder = "output"; //was getenv(OUTFOLDER);
// // // // //     
// // // // //   // READ DOUBLES FROM FILE ======== WAS declared as GLOBAL SCOPE, now NO MORE
// // // // // //   FemusInputParser<double> * runtime_double; //this line is not needed, it works nevertheless but it's not needed //so the brutal way to make it visible everywhere is to put the declaration OUTSIDE the function and to declare it with extern in all the files where it's needed
// // // // //                                           //non serve il singleton pattern per questo! serve solo chiamare il costruttore QUI NEL MAIN e non OUTSIDE
// // // // //   runtime_double = new FemusInputParser<double>("Doubles","./");
// // // // //   runtime_double->read();
// // // // //   runtime_double->print();
// // // // // 
// // // // //   // READ STRINGS FROM FILE ========
// // // // //   FemusInputParser<std::string> * runtime_string = new FemusInputParser<std::string>("Strings","./");  
// // // // //   runtime_string->read();
// // // // //   runtime_string->print();
// // // // //   
// // // // //   
// // // // //    //  OPEN BIG TABLE ========
// // // // //   std::ofstream ofs;
// // // // //   ofs.open("./output/BigTable.txt",std::fstream::app);
// // // // //   ofs /*<< outfolder << " =============="*/ << std::endl;
// // // // //   ofs  << std::left << std::setw(15) << std::setprecision(12) << outfolder << "   " << "E_frac" << "  " << "E_well" << "  "  << "K_frac" << "  " << "K_well" << "  " << "nlevs" << " *** " << "NUM_FLUX"  << "  " << "DEN_AVG_PRESS"  << "  "  << "P.I."  << "  "  <<  std::endl;
// // // // //   
// // // // //   unsigned short nm,nr;
// // // // //   std::cout << "#MULTIGRID levels? (>=1) \n";
// // // // //   nm = (int) runtime_double->get("nlevs");
// // // // // 
// // // // //   std::cout << "#MAX_REFINEMENT levels? (>=0) \n";
// // // // //   nr = (int) runtime_double->get("nrefins");
// // // // //   int tmp = nm;
// // // // //   nm += nr;
// // // // //   nr = tmp;
// // // // // 
// // // // //   cout << "nm  ==== " << nm << " nr ==== " << nr << std::endl;
// // // // //   
// // // // //   char *infile = new char [50];
// // // // //   std::ostringstream meshfile; meshfile << "./input/" << runtime_string->get("mesh_name");
// // // // //   sprintf(infile, meshfile.str().c_str() );
// // // // // 
// // // // //   //Adimensional quantity (Lref,Uref)
// // // // //   double Lref = 1.0;
// // // // //   double Uref = 1.0;
// // // // //   Parameter parameter(Lref,Uref);
// // // // // 
// // // // //   // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
// // // // //   Fluid fluid(parameter,1.,1.,"Newtonian");
// // // // //   Solid solid(parameter,runtime_double->get("young_well"),0.4,1000.,"Linear_elastic");
// // // // // //AAAAAAA: HERE you are setting 1 kg/m^3 for the FLUID DENSITY and 1000 for the SOLID DENSITY!!!
// // // // // // In this case we are not actually using those numbers, but PAY ATTENTION!  
// // // // //   
// // // // // //Steadystate MultiGrid
// // // // //   MultiLevelMesh ml_msh(nm,nr,infile,"fifth",Lref,SetRefinementFlag); 
// // // // // //   ml_msh.EraseCoarseLevels(2);
// // // // //   MyMultiGrid  mg(&ml_msh);
// // // // // 
// // // // // //i wanted to print the mesh without variables, we need to add a function for that
// // // // // //   std::vector<std::string> print_vars_tmp;
// // // // // //   print_vars_tmp.resize(1);
// // // // // //   print_vars_tmp[0] = "TMP";
// // // // // //   mg.printsol_vtu_inline("biquadratic",print_vars_tmp);
// // // // // // // //   // END MESH =================================
// // // // // 
// // // // // // PHYSICS ===========================
// // // // //   mg.AddParameters(runtime_double);  //TODO where do i need this?
// // // // //   mg.parameters.set<Fluid>("Fluid") = fluid;
// // // // //   mg.parameters.set<Solid>("Solid") = solid;
// // // // // 
// // // // // //Start System Variables;===========================
// // // // //   std::vector<std::string> varnames_p(NVAR_P);
// // // // //   varnames_p[0] = "p";
// // // // //   std::vector<std::string> varnames_d(NVAR_D);
// // // // //   varnames_d[0] = "DX";
// // // // //   varnames_d[1] = "DY";
// // // // //   varnames_d[2] = "DZ";
// // // // //   std::vector<std::string> varnames_u(NVAR_VEL);
// // // // //   varnames_u[0] = "UX";
// // // // //   varnames_u[1] = "UY";
// // // // //   varnames_u[2] = "UZ";
// // // // // // // //   //all these variables are added to the vector that is "spanned" by INDEX
// // // // // // // //   for (int i=0; i<NVAR_P; ++i)     mg.AddSolution(varnames_p[i].c_str(),"biquadratic"); //Add Solution does not exist anymore
// // // // // // // //   for (int i=0; i<NVAR_D; ++i)     mg.AddSolution(varnames_d[i].c_str(),"biquadratic");
// // // // // // // //   for (int i=0; i<NVAR_VEL; ++i)   mg.AddSolution(varnames_u[i].c_str(),"biquadratic");
// // // // // // // // 
// // // // // // // //   for (int i=0; i<NVAR_P; ++i)     mg.Initialize(varnames_p[i].c_str());
// // // // // // // //   for (int i=0; i<NVAR_D; ++i)     mg.Initialize(varnames_d[i].c_str());
// // // // // // // //   for (int i=0; i<NVAR_VEL; ++i)   mg.Initialize(varnames_u[i].c_str());
// // // // // // // // //  mg.Initialize("All"); //TODO would this "All" still be working?
// // // // // // // // 
// // // // // // // //   mg.AttachSetBoundaryConditionFunction(BoundaryND);
// // // // // // // //   
// // // // // // // //   //Set Boundary (update Dirichlet(...) function)
// // // // // // // //   for (int i=0; i<NVAR_P; ++i)   mg.GenerateBdc(varnames_p[i].c_str());
// // // // // // // //   for (int i=0; i<NVAR_D; ++i)   mg.GenerateBdc(varnames_d[i].c_str());
// // // // // // // //   for (int i=0; i<NVAR_VEL; ++i) mg.GenerateBdc(varnames_u[i].c_str());
// // // // //       //TODO do i have to generate them also for the velocity? by default they do nothing?
// // // // //       //it seems like it was necessary... so the default is not neumann "do nothing"?!
// // // // //    
// // // // //   //End System Variables; ==============================
// // // // // 
// // // // //   // START EQUATIONS =================================
// // // // //   //TODO use the new system-based-framework
// // // // // //   mg.AddPde("EQN_P");
// // // // // //   mg.AddPde("EQN_D");
// // // // // //   mg.AddPde("EQN_VEL");
// // // // //   
// // // // // //   mg.ClearSolPdeIndex();
// // // // // //   for (int i=0; i<NVAR_P; ++i)       mg.AddSolutionToSolPdeIndex("EQN_P",varnames_p[i].c_str());
// // // // // //   for (int i=0; i<NVAR_D; ++i)       mg.AddSolutionToSolPdeIndex("EQN_D",varnames_d[i].c_str());
// // // // // //   for (int i=0; i<NVAR_VEL; ++i)     mg.AddSolutionToSolPdeIndex("EQN_VEL",varnames_u[i].c_str());
// // // // //   
// // // // // //   mg.CreatePdeStructure();
// // // // // // // //   mg.BuildSparsityPattern();  //TODO this will not be needed now
// // // // //   
// // // // //  for (unsigned nonlin = 0; nonlin < runtime_double->get("nonlin_iter"); nonlin++) {
// // // // //    
// // // // // // // //   mg.ClearVankaIndex();  // create index of solutions to be to used in the Vanka Smoother
// // // // // // // //   mg.SetAssembleFunction(AssembleMatrixResP);
// // // // // // // //   for (int i=0; i<NVAR_P; ++i)       mg.AddToVankaIndex(varnames_p[i].c_str());
// // // // // // // //   mg.SetVankaSchurOptions(false);//(true,true,1);
// // // // // // // //   mg.SetSolverFineGrids("GMRES");
// // // // // // // //   mg.SetPreconditionerFineGrids("LU");
// // // // // // // //   mg.SetTolerances(1.e-12,1.e-20,1.e+50,20);
// // // // // // // //   mg.SetDimVankaBlock("All"/*3*/);                //2^lev 1D 4^lev 2D 8^lev 3D
// // // // // // // //   mg.FullMultiGrid(7,1,1);
// // // // // // // // 
// // // // // // // //   
// // // // // // // //   mg.ClearVankaIndex();  // create index of solutions to be to used in the Vanka Smoother
// // // // // // // //   mg.SetAssembleFunction(AssembleMatrixResD);
// // // // // // // //   mg.AddToVankaIndex("DX");
// // // // // // // //   mg.AddToVankaIndex("DY");
// // // // // // // //   mg.AddToVankaIndex("DZ");  //remember it must be in MGIndex first //i dont solve 3d comp in 2D
// // // // // // // //   mg.SetVankaSchurOptions(false);
// // // // // // // //   mg.SetSolverFineGrids("GMRES");
// // // // // // // //   mg.SetPreconditionerFineGrids("LU");
// // // // // // // //   mg.SetTolerances(1.e-12,1.e-20,1.e+50,20);
// // // // // // // //   mg.SetDimVankaBlock("All"/*3*/);                //2^lev 1D 4^lev 2D 8^lev 3D     //No DD:pass "All"
// // // // // // // //   mg.FullMultiGrid(7,1,1);
// // // // // // // //   
// // // // // // // // 
// // // // // // // //   mg.ClearVankaIndex();  // create index of solutions to be to used in the Vanka Smoother
// // // // // // // //   mg.SetAssembleFunction(AssembleMatrixResVel);
// // // // // // // //   mg.AddToVankaIndex("UX");
// // // // // // // //   mg.AddToVankaIndex("UY");
// // // // // // // //   mg.AddToVankaIndex("UZ");  //remember it must be in MGIndex first //i dont solve 3d comp in 2D
// // // // // // // //   mg.SetVankaSchurOptions(false);
// // // // // // // //   //Solver configuration for Temperature problem
// // // // // // // //   mg.SetSolverFineGrids("GMRES");
// // // // // // // //   mg.SetPreconditionerFineGrids("LU");
// // // // // // // //   mg.SetTolerances(1.e-12,1.e-20,1.e+50,20);
// // // // // // // //   mg.SetDimVankaBlock("All"/*3*/);                //2^lev 1D 4^lev 2D 8^lev 3D     //No DD:pass "All"
// // // // // // // //   mg.FullMultiGrid(7,1,1);
// // // // // // // // 
// // // // // // // //  
// // // // // // // //   double Num = mg.ComputeProductivityIndexNum( (int) runtime_double->get("inner") );  //face 1
// // // // // // // //   std::vector<double> integrals(mg.ComputeProductivityIndexDen(2));
// // // // // // // // 
// // // // // // // //   double Den = (integrals[0]/integrals[1] - runtime_double->get("p_in")  );  
// // // // // // // //   
// // // // // // // //   std::cout << " The productivity index is " << Num / Den << std::endl;
// // // // // // // //   
// // // // // // // // //move mesh according to displacement
// // // // // // // //   std::vector<std::string> mov_vars(varnames_d);  //copy constructor
// // // // // // // //   mg.SetMovingMesh(mov_vars);
// // // // // // // //   mg.IncreaseStep(nonlin);
// // // // // // // //   mg.printsol_vtu_binary(outfolder,"biquadratic");
// // // // // // // //
// // // // //    
// // // // // // // // //   if (nonlin == runtime_double.get("nonlin_iter") -1 ) {
// // // // // // // //   ofs  << std::left << std::setw(15) << std::setprecision(12) << outfolder << "   " << runtime_double->get("young_frac") << "  " << runtime_double->get("young_well") << "  " << runtime_double->get("kappa_frac") << "  "  << runtime_double->get("kappa_well") << "  "  << runtime_double->get("nlevs") << " *** " << Num << "   " << Den << "  " <<  Num / Den <<  "   " << std::endl;
// // // // // // // // //   }
// // // // //   
// // // // // }
// // // // //    
// // // // //    ofs << std::endl;
// // // // //    
// // // // // 
// // // // //   mg.clear();
// // // // // //   mg.FreeMultigrid(); //it does not exit anymore
// // // // //   delete [] infile;
// // // // // //   delete [] runtime_double; //check seg fault on these two... the global variable may not help
// // // // // //   delete [] runtime_string;
  
  return 0;
}


bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level) {
  bool refine=0;
  //refinement based on center element coordinate and level
  //   if( (y<(1.5-0.25 *level) && (x<0.25 || x >0.75)) ||
  // 	 (y<(1.5-0.25 *(level+1)) && x>0.25 && x <0.75) ||
  // 	  x<(0.5-0.25*(level+1)) || x>(0.5+0.25*(level+1))){
  //     refine=0;
  //   }
  //refinemenet based on Elemen Group Number
  const unsigned nlevs = runtime_double->get("nlevs");
  
  if (ElemGroupNumber== (int) runtime_double->get("frac_reg_idx") && level< 3+nlevs-1) refine=1;
  if (ElemGroupNumber== (int) runtime_double->get("layer1")       && level< 2+nlevs-1) refine=1;
  if (ElemGroupNumber== (int) runtime_double->get("layer2")       && level< 1+nlevs-1) refine=1;
  if (ElemGroupNumber== (int) runtime_double->get("layer3")       && level< 0+nlevs-1) refine=1;

  return refine;
}

/// 1 Dirichlet 0 Neumann - value: any nonhomogeneous BC
bool BoundaryND(/*MultiLevelProblem& mg_in,*/const double &x, const double &y, const double &z,const char name[], double &value, const int FaceName,const double time) {

  bool test=1; //Dirichlet
  bool DIR=1;
  bool NEU=0;
  value=0.;

//   MyMultiGrid* mg = static_cast<MyMultiGrid*>(&mg_in); //we use the "global" runtime map so far
  
  if (!strcmp(name,"p")) {

    if ( (int) runtime_double->get("inner") == FaceName) {
      test=DIR;
      value = runtime_double->get("p_in"); //well
    }
    else if ( (int) runtime_double->get("outer") == FaceName) {
      test=DIR;
      value = runtime_double->get("p_out");
    }

  }

  
  else if (!strcmp(name,"DX")) {

    if ( (int) runtime_double->get("inner") == FaceName) {
      test=DIR;
      value=0.;
    }
    else if ( (int) runtime_double->get("outer") == FaceName) {
      test=DIR;
      value=0.;
    }

  }

  else if (!strcmp(name,"DY")) {

    if ( (int) runtime_double->get("inner") == FaceName) {
      test=DIR;
      value=0.;
    }
    else if ( (int) runtime_double->get("outer") == FaceName) {
      test=DIR;
      value=0.;
    }

  }
  else if (!strcmp(name,"DZ")) {

    if ( (int) runtime_double->get("inner") == FaceName) {
      test=NEU;
      value=0.;
    }
    else if ( (int) runtime_double->get("outer") == FaceName) {
      test=NEU;
      value=0.;
    }

  }
  
   else if (!strcmp(name,"UX")) {

    if ( (int) runtime_double->get("inner") == FaceName) {
      test=NEU;
      value=0.;
    }
    else if ( (int) runtime_double->get("outer") == FaceName) {
      test=NEU;
      value=0.;
    }

  }

  else if (!strcmp(name,"UY")) {

    if ( (int) runtime_double->get("inner") == FaceName) {
      test=NEU;
      value=0.;
    }
    else if ( (int) runtime_double->get("outer") == FaceName) {
      test=NEU;
      value=0.;
    }

  }
  else if (!strcmp(name,"UZ")) {

    if ( (int) runtime_double->get("inner") == FaceName) {
      test=NEU;
      value=0.;
    }
    else if ( (int) runtime_double->get("outer") == FaceName) {
      test=NEU;
      value=0.;
    }

  } 

  return test;
}


int AssembleMatrixResP(MultiLevelProblem &mg2, unsigned level, /*const elem_type *type_elem[6][5], vector <vector <double> > &vt,*/ const unsigned &gridn, const unsigned &ipde, const bool &assemble_matrix) {

// not time dependent
  MyMultiGrid& mg = static_cast<MyMultiGrid&>(mg2);
  
  unsigned eqn_idx_current = IDX_P;

  
// // //   LinearSolverM*    lsyspde_lev     = mg.Lin_Solver_[level];
// // //   LinearSolverM*    lsyspdemesh_lev = mg.Lin_Solver_[level];
// // // 
// // //   PetscErrorCode ierr;
// // //   
// // // //global lin algebra
// // //   Mat&           myKK      = lsyspde_lev->KK;
// // //   Vec&           myRES     = lsyspde_lev->RES;
// // //   vector <int>& myKKIndex  = lsyspde_lev->KKIndex;
// // // 
// // // //geometry
// // //   PetscInt node_geom_el[MAX_EL_NODES];
// // //   double geom_el_coords[SPACEDIM][MAX_EL_NODES];
// // //   elem*    myel  = lsyspdemesh_lev->el;
// // //   unsigned nel   = lsyspdemesh_lev->GetElementNumber();
// // //   unsigned igrid = lsyspdemesh_lev->GetLevel();
// // // 
// // // // variables
// // //   PetscInt nodeP[MAX_EL_NODES];
// // //   unsigned indexP    = mg.GetMGIndex("p"); //block //corresponds to position of P in mg (set in main)
// // //   unsigned indexSolP = mg.GetIndex("p");
// // //   unsigned order_ind2 = lsyspde_lev->SolType[indexSolP];  //quadratic
// // //   unsigned end_ind2   = lsyspde_lev->END_IND[order_ind2];
// // // 
// // //   double B[NVAR_P][NVAR_P][MAX_EL_NODES*MAX_EL_NODES];
// // //   double F[NVAR_P][MAX_EL_NODES];
// // // 
// // // 
// // // // FE order
// // //   double phi2[MAX_EL_NODES],gradphi2[MAX_EL_NODES][SPACEDIM];
// // // 
// // // 
// // // // quadrature
// // //   double WeightXJac;
// // // 
// // // // Physics, and Math (operators)
// // //   double    _ILambda = .0e-3;
// // //   double    _IRe     = mg._fluid->get_IReynolds_number();
// // //   bool _NavierStokes = 1;
// // //   bool       _newton = 0*_NavierStokes;
// // //   bool      _penalty = 0;
// // // 
// // //   int indexSolD[NVAR_D];
// // //   std::string varnamesD[NVAR_D];
// // //   varnamesD[0] = "DX";
// // //   varnamesD[1] = "DY";
// // //   varnamesD[2] = "DZ";
// // // 
// // //    for (int i=0; i<NVAR_D; ++i)    indexSolD[i] = mg.GetIndex(varnamesD[i].c_str());  
// // // //   unsigned varindx=mg.GetIndex("DX");
// // //  
// // //  //==== Parameters ======
// // //  int moving_dom_p = mg._runtime_double->get("moving_dom_p");  
// // // 
// // //  //the following is a loop over all the system variables, because all of them may be involved
// // //   
// // // // vectors for getting dofs   ====================
// // //   int sol_size = lsyspde_lev->Sol_.size();
// // //   PetscScalar **vec_sol     = new PetscScalar*[sol_size];
// // //   vector<PetscVectorM*> petsc_vec_sol;
// // //   petsc_vec_sol.resize(sol_size);
// // // 
// // //   for (int k=0; k<sol_size; k++) {
// // //     petsc_vec_sol[k]  = static_cast<PetscVectorM*>(lsyspde_lev->Sol_[k]);
// // //     ierr = VecGetArray(petsc_vec_sol[k]->vec(),&vec_sol[k]);    CHKERRQ(ierr);
// // //   }
// // // // vectors for getting dofs   ====================
// // // 
// // //   ierr = MatZeroEntries(myKK);  CHKERRQ(ierr);
// // //   ierr = VecZeroEntries(myRES);  CHKERRQ(ierr);
// // //    
// // // 
// // //   clock_t AssemblyTime=0;
// // //   clock_t start_time, end_time;
// // //   start_time=clock();
// // // 
// // // 
// // //   for (unsigned kel=0;kel<nel;kel++) {
// // // 
// // //     short unsigned kelt = myel->GetElementType(kel);
// // //     unsigned       nve2 = myel->GetElementDofNumber(kel,end_ind2);
// // //     int        myregion = myel->GetElementGroup(kel);
// // //     
// // //     double kappa = mg._runtime_double-> get("kappa_well");
// // //     if ( myregion == (int) mg._runtime_double-> get("frac_reg_idx") ) kappa =  mg._runtime_double-> get("kappa_frac");
// // // 
// // //     memset(F[indexP],0,nve2*sizeof(double));
// // //     memset(B[indexP][indexP],0,nve2*nve2*sizeof(double));
// // // 
// // // // geometry and dof indices
// // //     for (unsigned i=0;i<nve2;i++) {
// // //       unsigned inode=myel->GetElementDofIndex(kel,i)-1u;
// // //       node_geom_el[i]=inode;
// // //     for (unsigned j=0; j<SPACEDIM; j++)    {  geom_el_coords[j][i] = vt[j][inode] + moving_dom_p*vec_sol[ indexSolD[j] ][inode]; }
// // //       nodeP[i] = node_geom_el[i] + myKKIndex[indexP];
// // //     }
// // // 
// // // 
// // //     if (igrid==gridn || !myel->GetRefinedElementIndex(kel)) {   //
// // //       
// // //       for (unsigned ig=0; ig < type_elem[kelt][order_ind2]->GetGaussPointNumber(); ig++) {
// // //         // *** get Jacobian and test function and test function derivatives ***
// // //        
// // //         (type_elem[kelt][order_ind2]->*type_elem[kelt][order_ind2]->Jacobian_ptr)(geom_el_coords,ig,WeightXJac,phi2,gradphi2);
// // // 
// // //         double SolP = 0;
// // //         double gradSolP[3]={0.,0.,0.};
// // //         for (unsigned i=0; i<nve2; i++) {
// // //           unsigned k=mg.GetIndex("p");
// // //           double soli=vec_sol[k][node_geom_el[i]];
// // //           SolP += phi2[i]*soli;
// // // 	  for (unsigned j=0; j<SPACEDIM; j++) {  gradSolP[j] += gradphi2[i][j]*soli;   }
// // //         }
// // // 
// // //         const double *gradfi=gradphi2[0];
// // //         const double *fi=phi2;
// // //         // *** phi_i loop ***
// // //         for (unsigned i=0; i<nve2; i++,gradfi+=3,fi++) {
// // // 
// // //           //BEGIN RESIDUALS A block ===========================
// // // 
// // //             double Lap_rhs =gradphi2[i][0]*gradSolP[0]+gradphi2[i][1]*gradSolP[1]+gradphi2[i][2]*gradSolP[2];
// // //             // residual equation U
// // //           F[indexP][i] += WeightXJac*(-0*phi2[i]
// // //                           -kappa*Lap_rhs
// // // // 			-_NavierStokes*(SolP*gradSolP[0] + SolV*gradSolP[1] + SolW*gradSolP[2])*phi2[i]
// // // // 			+SolP*gradphi2[i][0]
// // //                           );
// // // 
// // //           //END RESIDUALS A block ===========================
// // // 
// // //           const double *gradfj=gradphi2[0];
// // //           const double *fj=phi2;
// // //           // *** phi_j loop ***
// // //           for (unsigned j=0; j<nve2; j++,gradfj+=3,fj++) {
// // // 
// // //             // Laplacian
// // //             double Lap = ((*(gradfi+0))*(*(gradfj+0))+(*(gradfi+1))*(*(gradfj+1))+(*(gradfi+2))*(*(gradfj+2)));
// // //             // advection term I
// // // // 	    double Adv1 = ( SolP*(*(gradfj+0))*(*(fi))+ SolV*(*(gradfj+1))*(*(fi))+ SolW*(*(gradfj+2))*(*(fi)) );
// // // 
// // //             B[indexP][indexP][i*nve2+j] += WeightXJac*kappa*Lap;
// // // 
// // //           }   //end phij loop
// // //         } //end phii loop
// // // 
// // // 
// // //       }  // end gauss point loop
// // //     }  // endif single element not refined or fine grid loop
// // // 
// // // 
// // //     // U-equation
// // //     ierr = VecSetValues(myRES,nve2,nodeP,F[indexP],ADD_VALUES);    CHKERRQ(ierr);
// // //     ierr = MatSetValuesBlocked(myKK,nve2,nodeP,nve2,nodeP,B[indexP][indexP],ADD_VALUES);    CHKERRQ(ierr);
// // // 
// // //   } //end list of elements loop
// // // 
// // //   //BEGIN MATRIX ASSEMBLY ============
// // //   ierr = MatAssemblyBegin(myKK,MAT_FINAL_ASSEMBLY);  CHKERRQ(ierr);
// // //   ierr = MatAssemblyEnd(myKK,MAT_FINAL_ASSEMBLY);  CHKERRQ(ierr);
// // //   //END MATRIX ASSEMBLY   ============
// // // 
// // //   //BEGIN RESIDUAL ASSEMBLY ============
// // //   ierr = VecAssemblyBegin(myRES);  CHKERRQ(ierr);
// // //   ierr = VecAssemblyEnd(myRES);  CHKERRQ(ierr);
// // //   //END RESIDUAL ASSEMBLY ============
// // // 
// // // //  ierr= MatView(myKK, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
// // // 
// // // //   if(igrid==30){
// // // //     PetscViewer viewer;
// // // //     ierr=PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);CHKERRQ(ierr);
// // // //     ierr= MatView(myKK,viewer);CHKERRQ(ierr);
// // // //         ierr= VecView(myRES,viewer);CHKERRQ(ierr);
// // // //     double ff;
// // // //     std::cin>>ff;
// // //     
// // // //     system("sleep 5");
// // //     
// // // //   }
// // // 
// // //   // *************************************
// // //   end_time=clock();
// // //   AssemblyTime+=(end_time-start_time);
// // //   // ***************** END ASSEMBLY *******************
// // // 
// // //   // *** Computational info ***
// // // //  cout<<"Grid="<< lsyspdemesh_lev->GetLevel()<<endl;
// // // //  cout<<"ASSEMBLY + RESIDUAL time="<<static_cast<double>(AssemblyTime)/CLOCKS_PER_SEC<<endl;
// // // 
// // //   //free memory
// // //   for (int k=0; k<lsyspde_lev->Sol_.size(); k++) {
// // //     ierr = VecRestoreArray(petsc_vec_sol[k]->vec(),&vec_sol[k]);
// // //     CHKERRQ(ierr);
// // //     delete [] vec_sol[k];
// // //   }
// // //   delete [] vec_sol;

  return 1;

}


//ok, in the pressure part i don't have anything and then i fill, without preparing any sparsity pattern;
// then in the following equation i only filled the pressure part, and the rest was never filled,
//and i want to fill the displacement blocks, but for some reason it doesn't let me do that...
//ok it says new nonzero causes a malloc... of course, the sparsity pattern was never put into that!
//why nothing bad happened before in the pressure equation?! it seems like we need to call something for the matrix
// that we didn't call...

//new nonzero caused a malloc.. when can this happen?! if the sparsity pattern is built on-the-fly,
//of course i will always have mallocs!!
//why do certain things give a malloc and other things don't?!?
//what did i say about the matrix when i first initialized it? did i give at least a maximum number of columns,
// (without saying the positions of the nonzero columns) ?


//Physically speaking, we must take into account the couplings,
//and they must be somewhat IN ACCORDANCE: the permeabilities must be in agreement with the elastic constants...

int AssembleMatrixResD(MultiLevelProblem &mg2, unsigned level, /*const elem_type *type_elem[6][5], vector <vector <double> > &vt,*/ const unsigned &gridn, const unsigned &ipde, const bool &assemble_matrix) {

// not time dependent
// // //   MyMultiGrid& mg = static_cast<MyMultiGrid&>(mg2);
// // //   
// // //   unsigned eqn_idx_current = IDX_D;
// // // 
// // //    unsigned eqn_offset = 0;
// // //    for (unsigned i=0; i< eqn_idx_current; i++) { eqn_offset += mg._qty_ncomps[i];  }
// // //     
// // //   PetscErrorCode ierr;
// // // 
// // //   LinearSolverM*     lsyspde_lev = mg.Lin_Solver_[level];
// // //   LinearSolverM* lsyspdemesh_lev = mg.Lin_Solver_[level];
// // //   
// // //   Mat&           myKK      = lsyspde_lev->KK;
// // //   Vec&           myRES     = lsyspde_lev->RES;
// // //   vector <int>& myKKIndex  = lsyspde_lev->KKIndex;   //offset of the variables in the global matrix and rhs 
// // // 
// // //  
// // // // variables  
// // //   std::string varnamesD[NVAR_D];
// // //   varnamesD[0] = "DX";
// // //   varnamesD[1] = "DY";
// // //   varnamesD[2] = "DZ";
// // //   
// // //   unsigned order_ind2 = lsyspde_lev->SolType[mg.GetIndex( varnamesD[0].c_str() )];
// // //   unsigned end_ind2   = lsyspde_lev->END_IND[order_ind2];
// // //   
// // //   PetscInt nodeD[NVAR_D][MAX_EL_NODES];  //node+offset for variables
// // //   unsigned indexD[NVAR_D];
// // //   unsigned indexD_LocBlocks[NVAR_D];
// // //   
// // //   for (int i=0; i<NVAR_D; ++i) {   
// // //     indexD[i]           = mg.GetMGIndex(varnamesD[i].c_str());
// // //     indexD_LocBlocks[i] = indexD[i] - eqn_offset;
// // //   }
// // //   double B[NVAR_D][NVAR_D][MAX_EL_NODES*MAX_EL_NODES];
// // //   double Rhs[NVAR_D][MAX_EL_NODES];
// // // 
// // //   int indexSolD[NVAR_D];
// // //    for (int i=0; i<NVAR_D; ++i)    indexSolD[i] = mg.GetIndex(varnamesD[i].c_str());  
// // //   
// // //    //Parameters
// // //    int moving_dom_disp = mg._runtime_double->get("moving_dom_disp");
// // //    int Lap_flag = mg._runtime_double->get("Lap");
// // //    int Graddiv_flag = mg._runtime_double->get("Graddiv");
// // //   
// // //   // FE
// // //   double phi[MAX_EL_NODES],gradphi[MAX_EL_NODES][SPACEDIM];
// // //   
// // //   // quadrature
// // //   double WeightXJac,Weight_nojac;
// // // 
// // //   // Element, geometry
// // //   elem*           myel     = lsyspdemesh_lev->el;
// // //   double elem_coords[SPACEDIM][MAX_EL_NODES];
// // //   PetscInt node_geom_el[MAX_EL_NODES]; 
// // // 
// // //   // Physics
// // //   double g[SPACEDIM] = {0.,0.,0.};
// // //   double dt = 1.; //mg.GetTimeStep();
// // //   double eps_pen = 1.e40;
// // // 
// // //   double _rhof = mg._fluid->get_density();
// // //   double _IRe  = mg._fluid->get_IReynolds_number();
// // //   double _rhos = mg._solid->get_density();
// // //   
// // //   double     _mu_lame = mg._solid->get_lame_shear_modulus();
// // //   double _lambda_lame = mg._solid->get_lame_lambda();
// // //   
// // //   double _mus = _mu_lame/_rhof;
// // //   double _lambda  = _lambda_lame / _rhof;
// // //   double _betafsi = _rhos; //_rhos / _rhof;
// // //   double _betans   = 1.;
// // //   double _theta_ns = 1.0;
// // //   double _theta_sm = 1.0; //0.5;
// // //   double _lambda_map = 0.;
// // //   int    _solid_model = mg._solid->get_physical_model();
// // //   bool   _newton=0;
// // // 
// // //   bool   penalty = lsyspde_lev->GetStabilization();
// // // 
// // //   double C_mat[3][3][3][3];
// // //   double Cauchy[3][3];
// // //   double Cauchy_old[3][3];
// // //   double Jnp1_hat;
// // //   double Jn_hat;
// // //   double I_bleft;
// // //   double F[3][3];
// // //   double b_left[3][3];
// // //   double e[3][3];
// // //   double e_old[3][3];
// // //   double tg_stiff_matrix[3][3];
// // //   const double Id2th[3][3] = {
// // //     { 1.,0.,0.},
// // //     { 0.,1.,0.},
// // //     { 0.,0.,1.}
// // //   };
// // // 
// // //   //initialization: Saint-Venaint Kirchoff model
// // //   for (int I=0; I<3; ++I) {
// // //     for (int J=0; J<3; ++J) {
// // //       for (int K=0; K<3; ++K) {
// // //         for (int L=0; L<3; ++L) {
// // // // 	 C_mat[I][J][K][L] = 2.*_mus*Id2th[I][K]*Id2th[J][L]; //incompressible
// // //           C_mat[I][J][K][L] = _lambda*Id2th[I][J]*Id2th[K][L] + 2.*_mus*Id2th[I][K]*Id2th[J][L]; //compressible
// // // // 	 cout << C_mat[I][J][K][L] << endl;
// // //         }
// // //       }
// // //     }
// // //   }
// // // 
// // //   
// // // // take the dof vector ============================  
// // //   int sol_size = lsyspde_lev->Sol_.size(); //number of variables in sol
// // //   PetscScalar **vec_sol     = new PetscScalar*[sol_size];
// // //   vector<PetscVectorM*> petsc_vec_sol;
// // //   petsc_vec_sol.resize(sol_size);
// // // 
// // //   for (int k=0; k< sol_size; k++) {
// // //     petsc_vec_sol[k]  = static_cast<PetscVectorM*>(lsyspde_lev->Sol_[k]);
// // //     ierr = VecGetArray(petsc_vec_sol[k]->vec(),&vec_sol[k]);    CHKERRQ(ierr);  //here the input is the first, the output to be used later is the last
// // //   }
// // // // end take the dof vector ============================  
// // // 
// // //   clock_t AssemblyTime=0;
// // //   clock_t start_time, end_time;
// // // 
// // //   unsigned nel   = lsyspdemesh_lev->GetElementNumber();
// // //   unsigned igrid = lsyspdemesh_lev->GetLevel();
// // // 
// // //   start_time=clock();
// // //   
// // //   ierr = MatZeroEntries(myKK);  CHKERRQ(ierr);
// // //   ierr = VecZeroEntries(myRES);  CHKERRQ(ierr);
// // // 
// // //   
// // // //         PetscViewer viewer;
// // // //         ierr=PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);CHKERRQ(ierr);
// // // //         ierr= MatView(myKK,viewer);CHKERRQ(ierr);
// // // // //         ierr= VecView(myRES,viewer);CHKERRQ(ierr);
// // // //         double ff;
// // // //         std::cin>>ff;
// // // 
// // //  
// // //   for (unsigned kel=0;kel<nel;kel++) {
// // // 
// // //     short unsigned kelt = myel->GetElementType(kel);
// // //     unsigned nve        = myel->GetElementDofNumber(kel,end_ind2);
// // //     int myregion        = myel->GetElementGroup(kel);  //choose the material
// // // 
// // // // SECOND SOLID =======
// // //     //i could instantiate another solid here, but then i should change the _mus and _lambda in the equation...i'll do it here
// // //         double YoungFrac = mg._runtime_double->get("young_frac");
// // // 	double poisson = mg._solid->get_poisson_coeff();
// // //    if (myregion== (int) mg._runtime_double-> get("frac_reg_idx") ) {
// // //   _lambda_lame = (YoungFrac*poisson)/((1.+poisson)*(1.-2.*poisson));
// // //   _mu_lame     = YoungFrac/(2.*(1.+poisson));
// // //   _mus         = _mu_lame/_rhof;
// // //   _lambda      = _lambda_lame / _rhof;  //the nondimensional ones
// // //    }
// // // // END SECOND SOLID =======
// // //     
// // //     for (unsigned i=0; i<NVAR_D;i++) {
// // //       memset(Rhs[indexD_LocBlocks[i]],0,nve*sizeof(double));
// // //       for (unsigned j=0; j<NVAR_D; j++) {
// // //          memset(B[indexD_LocBlocks[i]][indexD_LocBlocks[j]],0,nve*nve*sizeof(double));
// // //       }
// // //     }
// // // 
// // //     for (unsigned i=0;i<nve;i++) {
// // //       unsigned inode = myel->GetElementDofIndex(kel,i)-1u;
// // //       for (unsigned j=0; j<SPACEDIM; j++)    elem_coords[j][i] = vt[j][inode] + moving_dom_disp*vec_sol[ indexSolD[j] ][inode];
// // //       node_geom_el[i] = inode;
// // //       for (unsigned j=0; j<NVAR_D; j++)       nodeD[ indexD_LocBlocks[j] ][i] = node_geom_el[i] + myKKIndex[indexD[j]];   //global
// // //     }
// // // 
// // // 
// // //     if (igrid==gridn || !myel->GetRefinedElementIndex(kel)) {
// // //  
// // //       for (unsigned igauss=0;igauss < type_elem[kelt][order_ind2]->GetGaussPointNumber(); igauss++) {
// // // 
// // //         (type_elem[kelt][order_ind2]->*type_elem[kelt][order_ind2]->Jacobian_ptr)(elem_coords,igauss,WeightXJac,phi,gradphi);         // *** get Jacobian and test function and test function derivatives ***
// // // 
// // // 	
// // //         double SolU[NVAR_VEL]={0.,0.,0.};
// // //         for (unsigned i=0; i<nve; i++) {
// // //           unsigned  k = mg.GetIndex("p");
// // //           double soli = vec_sol[k][node_geom_el[i]];
// // //           for (unsigned j=0; j<SPACEDIM; j++)  SolU[j] += - gradphi[i][j]*soli;  // drag term (- alpha*1/alpha grad P)
// // //          }
// // //          
// // //         double SolD[NVAR_D]={0.,0.,0.};
// // //         double gradSolD[NVAR_D][SPACEDIM] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
// // // 
// // // 	for (unsigned i=0; i<nve; i++) {
// // //              for (unsigned d=0; d<NVAR_D; d++) {
// // // 	       unsigned k = mg.GetIndex( varnamesD[d].c_str() );
// // // 	       
// // //           double soli = vec_sol[k][node_geom_el[i]];
// // //           SolD[d]+=phi[i]*soli;
// // // 	  for (unsigned j=0; j<SPACEDIM; j++)  gradSolD[d][j] += gradphi[i][j]*soli;
// // // 
// // // 	    }
// // // 	}
// // // 
// // // 
// // //         //--------------------------------------------------------------------
// // //         if (_solid_model==0) {
// // // 
// // //           //computation of the stress tensor
// // //           e[0][0] = gradSolD[0][0];
// // //           e[0][1] = 0.5*(gradSolD[0][1] + gradSolD[1][0]);
// // //           e[0][2] = 0.;
// // //           e[0][2] = 0.*0.5*(gradSolD[0][2] + gradSolD[2][0]);
// // // 
// // //           e[1][0] = 0.5*(gradSolD[1][0] + gradSolD[0][1]);
// // //           e[1][1] = gradSolD[1][1];
// // //           e[1][2] = 0.;
// // //           e[1][2] = 0.*0.5*(gradSolD[1][2] + gradSolD[2][1]);
// // // 
// // //           e[2][0] = 0.*0.5*(gradSolD[0][2] + gradSolD[2][0]);
// // //           e[2][1] = 0.*0.5*(gradSolD[1][2] + gradSolD[2][1]);
// // //           e[2][2] = 0.*gradSolD[2][2];
// // // 
// // // 
// // //           // Cauchy stress tensor
// // //           for (int irow=0; irow<3; ++irow) {
// // //             for (int jcol=0; jcol<3; ++jcol) {
// // //               //compressible
// // // // 	   Cauchy[irow][jcol] = _lambda*I_e*Id2th[irow][jcol] + 2*_mus*e[irow][jcol];
// // //               //incompressible
// // //               Cauchy[irow][jcol] = 2*_mus*e[irow][jcol];
// // //             }
// // //           }
// // // 
// // //         }  //end solid model 0
// // // 
// // //         else if (_solid_model==1) {
// // // 
// // //           //deformation gradient
// // //           F[0][0] = 1. + gradSolD[0][0];
// // //           F[0][1] = gradSolD[0][1];
// // //           F[0][2] = gradSolD[0][2];
// // // 
// // //           F[1][0] = gradSolD[1][0];
// // //           F[1][1] = 1. + gradSolD[1][1];
// // //           F[1][2] = gradSolD[1][2];
// // // 
// // //           F[2][0] = gradSolD[2][0];
// // //           F[2][1] = gradSolD[2][1];
// // //           F[2][2] = 1. + gradSolD[2][2];
// // // 
// // //           Jnp1_hat =     F[0][0]*F[1][1]*F[2][2] + F[0][1]*F[1][2]*F[2][0] + F[0][2]*F[1][0]*F[2][1]
// // //                        - F[2][0]*F[1][1]*F[0][2] - F[2][1]*F[1][2]*F[0][0] - F[2][2]*F[1][0]*F[0][1];
// // // 
// // //           // computation of the the three deformation tensor b
// // //           for (int I=0; I<3; ++I) {
// // //             for (int J=0; J<3; ++J) {
// // //               b_left[I][J]=0.;
// // //               for (int K=0; K<3; ++K) {
// // //                 //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
// // //                 b_left[I][J] += F[I][K]*F[J][K];
// // //               }
// // // //         Cauchy compressible
// // //               Cauchy[I][J] = (_mus/Jnp1_hat)*(b_left[I][J] - Id2th[I][J])
// // //                              + (_lambda/Jnp1_hat)*log(Jnp1_hat)*Id2th[I][J];
// // // //         Cauchy incompressible
// // // //            Cauchy[I][J] = (_mus/Jnp1_hat)*(b_left[I][J] - Id2th[I][J]);
// // //             }
// // //           }
// // // 
// // //           I_bleft = b_left[0][0] + b_left[1][1] + b_left[2][2];
// // // 
// // //           for (int ii=0; ii<3; ++ii) {
// // //             for (int jj=0; jj<3; ++jj) {
// // //               for (int kk=0; kk<3; ++kk) {
// // //                 for (int ll=0; ll<3; ++ll) {
// // //                   C_mat[ii][jj][kk][ll] = (_lambda/Jnp1_hat)*Id2th[ii][jj]*Id2th[kk][ll]
// // //                                           + 2.*((_mus - _lambda*log(Jnp1_hat))/Jnp1_hat)*Id2th[ii][kk]*Id2th[jj][ll];
// // // 
// // // //                      C_mat[ii][jj][kk][ll] = 2.*_mus*pow(Jnp1_hat,-1.6666666666666)*(
// // // //  	                        0.333333333333*I_bleft*Id2th[ii][kk]*Id2th[jj][ll]              //1/3*I_c*i
// // // // // 	                        +0.111111111111*I_C*Id2th[i][j]*Id2th[k][l]             //1/9*I_b*IxI
// // // // // 				-0.333333333333*b_left[i][j]*Id2th[k][l]                //-1/3*b*I
// // // // // 				-0.333333333333*Id2th[i][j]*b_left[k][l]                //-1/3*b*I
// // // // 				)
// // // //  	                     -SolP*(Id2th[ii][jj]*Id2th[kk][ll]-2.*Id2th[ii][kk]*Id2th[jj][ll] );  // -p(IxI-2i)
// // // 
// // //                 }
// // //               }
// // //             }
// // //           }
// // // 
// // //           //Old deformation gradient
// // // //          F[0][0] = 1. + gradSolOldhatDX[0];
// // // 
// // // 
// // //         }  //end solid model 1
// // //         //----------------------------------------------------------------------------------------------------------------------------
// // // 
// // //         /////////////
// // // 
// // // 
// // //         {
// // //           const double *gradfi=gradphi[0];
// // //           const double *fi=phi;
// // // 
// // //           /// *** phi_i loop ***
// // //           for (unsigned i=0; i<nve; i++,gradfi+=3,fi++) {
// // // 
// // // 
// // //             //BEGIN RESIDUALS A + Bt block ==========================
// // // 
// // //               for (int r=0; r<NVAR_D; ++r) {
// // // 
// // // 		double Div_d = 0.;
// // // 		             for (int j=0; j<3; ++j) { Div_d += gradSolD[ j ][ j ]  ;}
// // //             Rhs[ indexD_LocBlocks[r] ][i] += 
// // //                                      WeightXJac*(
// // // 				       +        SolU[ indexD_LocBlocks[r] ]*phi[i] /*DRAG TERM*/ 
// // // 				       +  _betafsi*g[ indexD_LocBlocks[r] ]*phi[i] 
// // //                                   - Lap_flag*_mus * ( gradSolD[indexD_LocBlocks[r]][0]*gradphi[i][0] + gradSolD[indexD_LocBlocks[r]][1]*gradphi[i][1] + gradSolD[indexD_LocBlocks[r]][2]*gradphi[i][2] )  //laplace
// // //                                    - Graddiv_flag*(_mus + _lambda)*( Div_d )*gradphi[i][indexD_LocBlocks[r]]
// // //                                    + 0.*(-_theta_sm*dt*(gradphi[i][0]*Cauchy[r][0] + gradphi[i][1]*Cauchy[r][1] + gradphi[i][2]*Cauchy[r][2]) )
// // //                                  );
// // // // 		    -(1-_theta_sm)*dt*(gradphi_old[i][0]*Cauchy_old[0][0] + gradphi_old[i][1]*Cauchy_old[0][1] + gradphi_old[i][2]*Cauchy_old[0][2])*Weight_old
// // // // 	            -phi[i]*_betafsi*(SolU - SolOldU)*Weight
// // // // 		    +dt*SolP*gradphi[i][0]*Weight
// // // // 		    +(1-_theta_sm)*dt*SolOldP*gradphi_old[i][0]*Weight_old
// // // 			 }
// // //             //-------------------------------------------------------------------------------------------------------------------------------
// // // 
// // // 
// // //             //END RESIDUALS A + Bt block ===========================
// // // 
// // //             const double *gradfj=gradphi[0];
// // //             const double *fj=phi;
// // //             // *** phi_j loop ***
// // //             for (unsigned j=0; j<nve; j++,gradfj+=3,fj++) {
// // // 
// // //               //Laplacian
// // //               double Lap = ((*(gradfi+0))*(*(gradfj+0))+(*(gradfi+1))*(*(gradfj+1))+(*(gradfi+2))*(*(gradfj+2)));
// // // 
// // //               //now only 2D
// // //               for (int icount=0; icount<2; ++icount) {
// // //                 for (int jcount=0; jcount<2; ++jcount) {
// // //                   tg_stiff_matrix[icount][jcount] = 0.;
// // //                   for (int kcount=0; kcount<2; ++kcount) {
// // //                     for (int lcount=0; lcount<2; ++lcount) {
// // //                       tg_stiff_matrix[icount][jcount] += (*(gradfi+kcount))*0.25*(
// // //                                                            C_mat[icount][kcount][jcount][lcount]+C_mat[icount][kcount][lcount][jcount]
// // //                                                            +C_mat[kcount][icount][jcount][lcount]+C_mat[kcount][icount][lcount][jcount]
// // //                                                          )*(*(gradfj+lcount));
// // //                     }
// // //                   }
// // //                 }
// // //               }
// // // 
// // //               /// Stiffness operator -- Elasticity equation (Linear or not)
// // //              for (int r=0; r<2; ++r) {
// // //                 for (int c=0; c<2; ++c) {
// // // 		  if (c==r){
// // //               B[ indexD_LocBlocks[r] ][ indexD_LocBlocks[c] ][i*nve+j] += WeightXJac*(
// // // 		                                                                     + Lap_flag*_mus*Lap                                  //_theta_sm*dt*tg_stiff_matrix[r][c]*Weight;
// // // 	                                                                              + Graddiv_flag*(_mus + _lambda)*gradphi[i][r]*gradphi[j][c]  //so check that it is ADDING to the LAPLACIAN
// // // 										     );                          
// // // 		  }
// // // 		    
// // // 		  }
// // // 	     }
// // // 
// // //               B[ indexD_LocBlocks[2] ][ indexD_LocBlocks[2] ][i*nve+j] += WeightXJac*(
// // // 		                                                                  Lap_flag*_mus*Lap
// // //                                                                             + Graddiv_flag*(_mus + _lambda)*gradphi[i][2]*gradphi[j][2]
// // // 										     );
// // // 
// // //             }
// // //           }
// // //         }
// // //         ////////////
// // // 
// // //       }
// // //     }
// // // 
// // //     for (unsigned i=0;i<NVAR_D;i++) {
// // //       ierr = VecSetValues(myRES,nve,nodeD[ indexD_LocBlocks[i] ],Rhs[ indexD_LocBlocks[i] ],ADD_VALUES);      CHKERRQ(ierr);
// // //       for (unsigned j=0;j<NVAR_D;j++) {
// // //         ierr = MatSetValuesBlocked(myKK,nve,nodeD[ indexD_LocBlocks[i] ],nve,nodeD[ indexD_LocBlocks[j] ],B[ indexD_LocBlocks[i] ][ indexD_LocBlocks[j] ],ADD_VALUES);        CHKERRQ(ierr);
// // //       }
// // //     }
// // // 
// // //    } //end list of elements loop
// // // 
// // // //----------------------------------------------------------------------
// // // 
// // //   //BEGIN MATRIX ASSEMBLY ==========
// // //   ierr = MatAssemblyBegin(myKK,MAT_FINAL_ASSEMBLY);  CHKERRQ(ierr);
// // //   ierr = MatAssemblyEnd(myKK,MAT_FINAL_ASSEMBLY);  CHKERRQ(ierr);
// // //   //END MATRIX ASSEMBLY ============
// // // 
// // //   //BEGIN RESIDUAL ASSEMBLY ==========
// // //   ierr = VecAssemblyBegin(myRES);  CHKERRQ(ierr);
// // //   ierr = VecAssemblyEnd(myRES);  CHKERRQ(ierr);
// // //   //END RESIDUAL ASSEMBLY ============
// // //   
// // //   
// // // //         PetscViewer viewer;
// // // //         ierr=PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);CHKERRQ(ierr);
// // // //         ierr= MatView(myKK,viewer);CHKERRQ(ierr);
// // // // //         ierr= VecView(myRES,viewer);CHKERRQ(ierr);
// // // //         double ff;
// // // //         std::cin>>ff;
// // // 
// // // 
// // //   //free memory
// // //   for (int k=0; k<lsyspde_lev->Sol_.size(); k++) {
// // //     ierr = VecRestoreArray(petsc_vec_sol[k]->vec(),&vec_sol[k]);    CHKERRQ(ierr);
// // //     delete [] vec_sol[k];
// // //   }
// // //   delete [] vec_sol;
// // // 
// // //   // *************************************
// // //   end_time=clock();
// // //   AssemblyTime+=(end_time-start_time);
// // //   // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

//----------------------------------------------------------------------------------------------

  return 1;

}


int AssembleMatrixResVel(MultiLevelProblem &mg2, unsigned level, /*const elem_type *type_elem[6][5], vector <vector <double> > &vt,*/ const unsigned &gridn, const unsigned &ipde, const bool &assemble_matrix) {

// not time dependent
// // //   MyMultiGrid& mg = static_cast<MyMultiGrid&>(mg2);
// // //   
// // //   unsigned eqn_idx_current = IDX_VEL;
// // // 
// // //   unsigned eqn_offset = 0;
// // //   for (unsigned i=0; i< eqn_idx_current; i++) { eqn_offset += mg._qty_ncomps[i];  }
// // // 
// // //   
// // //   LinearSolverM*    lsyspde_lev     = mg.Lin_Solver_[level];
// // //   LinearSolverM*    lsyspdemesh_lev = mg.Lin_Solver_[level];
// // // 
// // //   PetscErrorCode ierr;
// // //   
// // // //global lin algebra
// // //   Mat&           myKK      = lsyspde_lev->KK;
// // //   Vec&           myRES     = lsyspde_lev->RES;
// // //   vector <int>& myKKIndex  = lsyspde_lev->KKIndex;
// // // 
// // // //geometry
// // //   double vx[SPACEDIM][MAX_EL_NODES];
// // //   PetscInt node_geom[MAX_EL_NODES];
// // //   elem*    myel  = lsyspdemesh_lev->el;
// // //   unsigned nel   = lsyspdemesh_lev->GetElementNumber();
// // //   unsigned igrid = lsyspdemesh_lev->GetLevel();
// // // 
// // // // variables  
// // //   std::string varnames[NVAR_VEL];
// // //   varnames[0] = "UX";
// // //   varnames[1] = "UY";
// // //   varnames[2] = "UZ";
// // // 
// // //   unsigned order_ind2 = lsyspde_lev->SolType[mg.GetIndex( varnames[0].c_str() )];
// // //   unsigned end_ind2   = lsyspde_lev->END_IND[order_ind2];
// // //   
// // //   PetscInt nodeU[NVAR_VEL][MAX_EL_NODES];  //node+offset for variables (with offset)
// // //   unsigned indexU[NVAR_VEL];
// // //   unsigned indexU_localBlocks[NVAR_VEL];
// // // 
// // //   int indexSolD[NVAR_D];
// // //   std::string varnamesD[NVAR_D];
// // //   varnamesD[0] = "DX";
// // //   varnamesD[1] = "DY";
// // //   varnamesD[2] = "DZ";
// // // 
// // //    for (int i=0; i<NVAR_D; ++i)    indexSolD[i] = mg.GetIndex(varnamesD[i].c_str());  
// // //   
// // // 
// // //   for (int i=0; i<NVAR_VEL; ++i) {
// // //     indexU[i] = mg.GetMGIndex(varnames[i].c_str());
// // //     indexU_localBlocks[i] = indexU[i] - eqn_offset;
// // //   }
// // //     
// // //   double B[NVAR_VEL][NVAR_VEL][MAX_EL_NODES*MAX_EL_NODES];
// // //   double F[NVAR_VEL][MAX_EL_NODES];
// // // 
// // // // FE
// // //   double phi2[MAX_EL_NODES],gradphi2[MAX_EL_NODES][SPACEDIM];
// // // 
// // // // quadrature
// // //   double Weight2XJac;
// // // 
// // // //Parameters
// // //   int moving_dom_vel = mg._runtime_double->get("moving_dom_vel");
// // //   
// // // //the following is a loop over all the system variables, because all of them may be involved  
// // //   // vectors for getting dofs   ====================
// // //   int sol_size = lsyspde_lev->Sol_.size();
// // //   PetscScalar **vec_sol     = new PetscScalar*[sol_size];
// // //   vector<PetscVectorM*> petsc_vec_sol;
// // //   petsc_vec_sol.resize(sol_size);
// // // 
// // //   for (int k=0; k<sol_size; k++) {
// // //     petsc_vec_sol[k]  = static_cast<PetscVectorM*>(lsyspde_lev->Sol_[k]);
// // //     ierr = VecGetArray(petsc_vec_sol[k]->vec(),&vec_sol[k]);    CHKERRQ(ierr);
// // //   }
// // // // vectors for getting dofs   ====================
// // // 
// // //   ierr = MatZeroEntries(myKK);  CHKERRQ(ierr);
// // //   ierr = VecZeroEntries(myRES);  CHKERRQ(ierr);
// // // 
// // //   clock_t AssemblyTime=0;
// // //   clock_t start_time, end_time;
// // //   start_time=clock();
// // // 
// // // 
// // //   for (unsigned kel=0;kel<nel;kel++) {
// // // 
// // //     short unsigned kelt = myel->GetElementType(kel);
// // //     unsigned       nve2 = myel->GetElementDofNumber(kel,end_ind2);
// // //     int        myregion = myel->GetElementGroup(kel);
// // //     
// // //     double kappa = mg._runtime_double->get("kappa_well");
// // //     if (myregion == (int) mg._runtime_double-> get("frac_reg_idx")) kappa = mg._runtime_double->get("kappa_frac");
// // // 
// // //     for (unsigned i=0; i<NVAR_VEL;i++) {
// // //       memset(F[indexU_localBlocks[i]],0,nve2*sizeof(double));
// // //       for (unsigned j=0; j<NVAR_VEL; j++) {
// // //         memset(B[indexU_localBlocks[i]][indexU_localBlocks[j]],0,nve2*nve2*sizeof(double));
// // //       }
// // //     }    
// // // 
// // //     for (unsigned i=0;i<nve2;i++) {
// // //       unsigned inode = myel->GetElementDofIndex(kel,i)-1u;
// // //       for (unsigned j=0; j<SPACEDIM; j++)    vx[j][i] = vt[j][inode] + moving_dom_vel*vec_sol[ indexSolD[j] ][inode];
// // //       node_geom[i] = inode;
// // //       for (unsigned j=0; j<NVAR_VEL; j++)       nodeU[ indexU_localBlocks[j] ][i] = node_geom[i] + myKKIndex[indexU[j]];   //global
// // //     }    
// // //     
// // // 
// // //     if (igrid==gridn || !myel->GetRefinedElementIndex(kel)) {
// // //       
// // //       for (unsigned ig=0; ig < type_elem[kelt][order_ind2]->GetGaussPointNumber(); ig++) {
// // // 	
// // //         (type_elem[kelt][order_ind2]->*type_elem[kelt][order_ind2]->Jacobian_ptr)(vx,ig,Weight2XJac,phi2,gradphi2);        // *** get Jacobian and test function and test function derivatives ***
// // // 
// // // 
// // //         double SolP=0;
// // //         double gradSolP[SPACEDIM]={0.,0.,0.};
// // //         for (unsigned i=0; i<nve2; i++) {
// // //           unsigned k=mg.GetIndex("p");
// // //           double soli=vec_sol[k][node_geom[i]];
// // //           SolP+=phi2[i]*soli;
// // // 	  for (unsigned j=0; j<SPACEDIM; j++) {  gradSolP[j]+=gradphi2[i][j]*soli;   }
// // //         }
// // // 
// // //          double SolU[NVAR_VEL]={0.,0.,0.};
// // //         for (unsigned i=0; i<nve2; i++) {
// // // 	  for (unsigned j=0; j<NVAR_VEL; j++) {
// // // 	  unsigned k=mg.GetIndex(varnames[j].c_str());
// // //           double soli=vec_sol[k][node_geom[i]];
// // //           SolU[j]+=phi2[i]*soli;
// // // 	  }
// // //         }       
// // //         
// // // // 	  for (unsigned j=0; j<SPACEDIM; j++) {  std::cout << SolU[j] << " == "<< std::endl; }
// // // 	    
// // // 	const double *fi=phi2;
// // //         const double *gradfi=gradphi2[0];
// // // 	
// // //         for (unsigned i=0; i<nve2; i++,gradfi+=3,fi++) {
// // // 	  
// // //               for (int r=0; r<NVAR_VEL; ++r) { F[ indexU_localBlocks[r] ][i]+=( (-kappa*gradSolP[r]- SolU[r])*phi2[i] )*Weight2XJac; }
// // // 
// // //           const double *fj=phi2;
// // //           const double *gradfj=gradphi2[0];
// // // 
// // // 	  for (unsigned j=0; j<nve2; j++,gradfj+=3,fj++) {
// // //              for (int r=0; r<NVAR_VEL; ++r) {
// // //                B[ indexU_localBlocks[r] ][ indexU_localBlocks[r] ][i*nve2+j] += phi2[i]*phi2[j]*Weight2XJac;
// // // 		}
// // //             }   //end phij loop
// // //          } //end phii loop
// // // 
// // //       }  // end gauss point loop
// // //     }  // endif single element not refined or fine grid loop
// // // 
// // //    for (unsigned i=0;i<NVAR_VEL;i++) {
// // //       ierr = VecSetValues(myRES,nve2,nodeU[indexU_localBlocks[i]],F[indexU_localBlocks[i]],ADD_VALUES);      CHKERRQ(ierr);
// // //       for (unsigned j=0;j<NVAR_VEL;j++) {
// // //         ierr = MatSetValuesBlocked(myKK,nve2,nodeU[indexU_localBlocks[i]],nve2,nodeU[indexU_localBlocks[j]],B[indexU_localBlocks[i]][indexU_localBlocks[j]],ADD_VALUES);        CHKERRQ(ierr);
// // //       }
// // //     }
// // //  
// // //   } //end list of elements loop
// // // 
// // //   // MATRIX ASSEMBLY ============
// // //   ierr = MatAssemblyBegin(myKK,MAT_FINAL_ASSEMBLY);  CHKERRQ(ierr);
// // //   ierr = MatAssemblyEnd(myKK,MAT_FINAL_ASSEMBLY);  CHKERRQ(ierr);
// // //   // RESIDUAL ASSEMBLY ============
// // //   ierr = VecAssemblyBegin(myRES);  CHKERRQ(ierr);
// // //   ierr = VecAssemblyEnd(myRES);  CHKERRQ(ierr);
// // // 
// // // //   if(igrid==30){
// // // //     PetscViewer viewer;
// // // //     ierr=PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);CHKERRQ(ierr);
// // // //     ierr= MatView(myKK,viewer);CHKERRQ(ierr);
// // // //     double ff;
// // // //     std::cin>>ff;
// // // //   }
// // // 
// // //   // *************************************
// // //   end_time=clock();
// // //   AssemblyTime+=(end_time-start_time);
// // //   // ***************** END ASSEMBLY *******************
// // // 
// // //   // *** Computational info ***
// // // //  cout<<"Grid="<< lsyspdemesh_lev->GetLevel()<<endl;
// // // //  cout<<"ASSEMBLY + RESIDUAL time="<<static_cast<double>(AssemblyTime)/CLOCKS_PER_SEC<<endl;
// // // 
// // //   //free memory
// // //   for (int k=0; k<lsyspde_lev->Sol_.size(); k++) {
// // //     ierr = VecRestoreArray(petsc_vec_sol[k]->vec(),&vec_sol[k]);
// // //     CHKERRQ(ierr);
// // //     delete [] vec_sol[k];
// // //   }
// // //   delete [] vec_sol;

  
  return 1;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// COMMENTS ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


// =====================
// You can do git clone even if the original repositories is messed up, you don't take what's uncommitted

// ========== SCRIPT, LAUNCH DEBUGGER WITH IT ====================
// The script is launched in parallel with kdevelop...
// We should tell kdevelop to launch the script when doing the debugging...
// But you also need to tell kdevelop what is the EXECUTABLE to DEBUG. In fact, the focus 
// is the EXECUTABLE when doing the debugging!!!
// maybe i should launch kdbg IN THE SCRIPT

// =========== GOOD MESHING ============
// If you want to do a good mesh, try to make "pseudo square" elements, as square as possible

// WHEN YOU COPY FACES/EDGES/WHATEVER by ROTATION, pay attention that you have to REMOVE DOUBLE EDGES!!!

// In gambit, i have to learn more about the COLORS of the EDGES (yellow, pale blue, pink, ..)

// If you have to do "divergent" "non-square" faces, then it's better if you do it AT THE COARSEST LEVEL as possible

// If you have some SYMMETRIES, maybe it's better to exploit them AFTER DOING THE MESH, so that 
// you copy GEOMETRY and MESH at the same time.
// CAN WE EVEN DO the SYMMETRIES AFTER DEFINING the BOUNDARY AND VOLUME ZONES?!?

// WHEN YOU DEFINE the ZONES, THEY DON'T HAVE TO BE INTERCONNECTED

// Differenza tra COPY MESH LINKED e COPY MESH UNLINKED:
// You have to do UNLINKED, so that the mesh is not linked to a specific geometrical entity. 
// Otherwise you cannot do CONNECT EDGES

// When you do SAVE, then you don't have any UNDO, they reset them.

// PAY ATTENTION to REMOVE MULTIPLE EDGES and MULTIPLE FACES.


// ======== GIT REVISION ========
// My idea is that BEFORE ANY COMMIT you have to run the FORMATTER so that
// all the files comply to that standard...
// but astyle is for C/C++ files...
// It does not work for the other file types in a package
// Makefiles, scripts (shell, python, whatever)


// ==============
// Un fluido puo' essere considerato come un continuo PLASTICO


// ====== REDIRECT STD OUT and STD ERR =============
// I guess that to every launched executable the operating system
// associates some DEFAULT FILE DESCRIPTORS, 0, 1 and 2.
// Therefore i would say, even if the executable links to other executables,
// the 0,1,2 file descriptors should be UNIQUE, be them with printf or with std::cout...


// =========== LINUX, sad observation =====
// Sometimes even linux can be easily fixed with a RESTART...
// I had a problem with repositories... even though i refreshed them, i was still getting
// an error about a /var/tmp... file that did not exist...
// So i could not solve it, then i restarted and everything went fine...


// ========= LINUX - ENVIRONMENT VARIABLES on suspend to ram =============
// They are not cancelled, they are still there!
// They are indeed cancelled when you RESTART the system!


// =========== MESH STRATEGIES ======================
// Se fai SUBTRACT between faces devi RIMUOVERE gli EDGE DOPPI!
// Lo SPLIT di una faccia si puo' fare O CON FACCE O CON EDGE... pero' gli edge devono essere virtuali...



// ============ PROBLEMA PRESSIONE ZERO in REGIONI 5,6,7 =================
// la matrice totale e' assemblata senza problemi. Probabilmente solo il pezzo della regione 8 viene estratto
// voglio che non ci sia NE' multigriglia ne'Adaptive mesh refinement ne' vanka.

// Qual e' l'opzione per non usare il MULTIGRIGLIA?  credo nlevs = 1
// Qual e' l'opzione per non usare, nel multigriglia, il VANKA SMOOTHER?  credo SetDimVankaBlock...
// Qual e' l'opzione per non usare l'ADAPTIVE MESH REFINEMENT?  credo nrefins

// Queste tre cose sono completamente indipendenti: posso fare raffinamenti in determinate zone
// restando ad un livello multigriglia, e SENZA Vanka smoother (quando sono ad un livello che solver usa,
// quello settato dallo SMOOTHER o quello settato dal COARSE LEVEL?)

// Allora, se ho degli ZERI, le questioni possono essere:
// IL RHS e' ZERO.
// ci sono delle bc bastarde
// la stampa e' bastarda...

// tolgo il multigriglia, tolgo l'adaptive mesh refinement, tolgo tutto, e voglio vedere nonzero.

// Il sospetto e' sul MESH. Sembra che nel mesh vecchio non ci siano gli edge doppi.

// IL CODICE DOVREBBE NON GIRARE SE GLI EDGE SONO DOPPI.
// E invece gira, semplicemente hai un edge doppio e quindi in pratica
// i due elementi ad esso adiacenti sono NON-COMUNICANTI, perche' i nodi comuni sono riprodotti due volte distintamente!!!


// ========= REINSTALLAZIONE PETSC =======
// Ho reinstallato le petsc aggiungendo il link a X11, ma lasciando gli stessi nomi...
// libmesh funzionera' ancora?! credo che l'installazione di petsc sia  HARD CODED in libmesh,
// con il nuovo sistema di configurazione. Con il vecchio sistema si includevano direttamente 
// i Makefile di petsc

// ======== DIFF with VERY DIFFERENT REPOS, which once were the same... =========
// mi conviene sempre fare un MEGA-commit, 
// fare git-fetch,
// e poi fare  un MEGA-diff tra HEAD e FETCH_HEAD 
// cercando di capire che cosa e' cambiato prima di fare il merge...
// Dopo si puo' fare un merge  di prova (esiste merge con dry-run?),
// e poi eventualmente gestire a mano i conflitti.

// AUTOMATIC DESTRUCTION of OBJECTS: Smart Pointers, o Classe con destructor =======

// =============== GIT CLONE o GIT PULL ============
// Quando importo da qualche repository, o perche' lo inizio o perche' lo aggiorno,
// importo che cosa? TUTTI i BRANCHES o solo il master? Importo anche i TAGS?
// Dipende chiaramente da che cosa gli dici.
// In particolare per GIT PULL, puoi specificare il REPOSITORY e la REFSPEC 
// di quel repository, tipicamente un BRANCH di quel repository.
// Se non specifichi niente, va a leggere dai file di configurazione di git.
// Nel dubbio e' sempre meglio specificare la refspec che vado a prendere dal repository remoto.

// Una domanda e': quando faccio git CLONE, i tag vengono importati? credo di si'

// ======= GIT CHECKOUT =============
// Se ho un BARE repository e faccio "git checkout master", sembra che non si "palesi", ma mi aspettavo che si palesasse...

// ======= I MIEI SCRIPT DI CONFIGURAZIONE =========
// In sostanza quello che fanno i miei script di configurazione
// e' evitare di MODIFICARE i FILE, 
// sia i veri e propri SORGENTI (.C, .h), sia i MAKEFILE...
// pero' QUALCOSA DI HARD CODED ORA CE L'HO, in quanto uso
// git checkout con le TAG.
// Prima ho cercato di EVITARE TUTTO L'HARD-CODED usando
// essenzialmente le ENVIRONMENT VARIABLES,
// che vengono lette dai MAkefile,
// oppure passate ai sorgenti tramite la macro -D in riga di compilazione...
// Le due principali categorie di file da mantenere sono, come dicevo:
// - Sorgenti
// - Makefile
// Poi ovviamente anche
// - Script di configurazione
// - Documentazione

// Io faccio in modo da non dover riscrivere NE' i MAKEFILE, NE' i SORGENTI (.hpp, .cpp).

// Quando usero' CMake mi rassegnero' all'HARD CODING...
// ogni volta che cambia la libreria dovro' ri-runnare configure e make,
// ma in effetti questo lo faccio gia'.
// In una macchina infatti ho un certo PETSC, un certo MPI, un certo HDF5,
// in un'altra ho un altro PETSC, un altro MPI, un altro HDF5,
// e il mio configure "sempliciotto" fa lo switch tra i due set.
// Quindi la procedura di riconfigurare la faccio gia',
// l'unica differenza e' che faccio in modo che essa non comporti RISCRITTURA DI FILE.
// In realta', ora che sto facendo checkout tra due diverse versioni,
// sto proprio reintroducendo la RISCRITTURA dei FILE.

// Se nel futuro voglio introdurre modifiche che vanno bene 
// sia per  la 3.3.6 sia per la 3.4.0, 
// non dovro' piu' direi usare i tag su due commit consecutivi, 
// ma portare avanti due BRANCH PARALLELI DIREI da switchare...:

//          ----- x ---- 3.3.6
// --x--x---
//          ----- x ---- 3.4.0

// Il problema e' che, se ogni futuro commit va bene per entrambi,
// dovrei fare piuttosto un "rebase", in modo che i commit nuovi
// si inseriscano "da sotto", come un'unghia che cresce dalla radice 
// (lo spostamento e' imposto sotto).
// il rebase pero' rischia di compromettermi la storia
// se qualcuno esterno fa una pull request...

// Oppure potrei crescere "dalla punta" di uno dei due branch, 
// e mergiare sempre verso l'altro branch SOLO i COMMIT "SUPERIORI allo SWITCH".
// Bisogna vedere se e' possibile fare un merge di due pezzi che non arrivano alla biforcazione upstream, ma "volanti"...

// In ogni modo, la cosa migliore sarebbe mantenere il supporto per ENTRAMBE le librerie 
// tramite uno script di configurazione che sia capace di RISCRITTURA.
// Quindi alla fine, siccome librerie diverse possono avere chiamate diverse, 
// il punto sara' fare degli #ifdef con le macro che tengono traccia delle   VERSIONI 3.3.6 o 3.4.0.

//E' chiaro infatti che non si puo' andare avanti su due binari paralleli. Il supporto alla 3.4.0
// non si deve dimenticare del supporto alle precedenti

// =========== SUPPORTO per "DIVERSI COMPILATORI" simile a SUPPORTO PER DIVERSI "SET DI LIBRERIE" ===========
// Libmesh ora gestisce la compilazione con diversi compilatori
// come tante sottocartelle di build.
// E' simile a quello che fa PETSC nell'avere tante sottocartelle
// di build corrispondenti a diverse PETSC_ARCH.
// Possiamo considerare ogni cartella come corrispondente 
// ad un diverso SET DI LIBRERIE. Il "compilatore" e' una diversa libreria tra le tante.
// Per esempio potrei fare una compilazione con {OpenMPI, HDF5},
// e un'altra con {MPICH, VTK}.
// Il compilatore e' una "LIBRERIA SPECIALE" nel senso che e' la "PRIMA LIBRERIA",
// in quanto le altre, se non precompilate, verranno compilate con quel compilatore scelto inizialmente
// (Gli external packages di petsc sono gia' precompilati oppure li compila sul momento?)
// [ Non so bene quali problemi sorgano quando si linkano librerie compilate con compilatori diversi ]




// ======== KEEP TRACK OF LINEAR ALGEBRA CLASSES ==========
// The best would maybe be to see the EVOLUTION of LIBMESH VERSIONS,
// and see the diffs in the libmesh commits (now in git)
// Another thing would be to do the diff between femTTU and libmesh

// NOW that we moved  to petsc 3.4.0, we must pay attention to the fact that
// some things are HARD CODED in the files, and this hard coding is made by
// the configure script
// if you link to two different versions of the same library, you cannot think 
// of switching automatically.
// You would need to COMPILE your code TWICE, once with ONE VERSION of the library
// and once with ANOTHER VERSION.
// You also need to compile several times your library if it is linked with 
// the same library but different compiling modes, can you imagine...?
// So, the point is that your library must be written in such a way that
// it may support SEVERAL VERSIONS of EACH EXTERNAL LIBRARY USED.
// In each installation, ONE AND ONLY ONE version will be used.
// For instance, in my case, I would like to be able to use  petsc 3.3.6 in one computer,
// and petsc 3.4.0 in another, just by setting different variables in the CONFIGURATION SCRIPT.

// You have to remember that after running a ./configure script everything is HARD CODED,
// so there is no way one could switch automatically between one library and another.

// In my case, I basically did two separate commits in the library,
// one associated to 3.3.6, another associated to 3.4.0.
// Therefore, SINCE I AM RECOMPILING ALL THE LIBRARY EVERY TIME,
//  I could just checkout one commit or the other and work in different computers.
// I could put two tags. What is the difference between switching between TAGS and between BRANCHES?
// Arent tags another type of branch?


// how to shrink two commits into one:    git rebase -i HEAD~2


// one thing that is important in git is: ONCE YOU SHARE, YOU CANNOT MODIFY LATER, so think twice before a pull request!
// - learn the TYPES of OBJECTS in git: heads, tags, replaces, whatever...
// http://gitready.com/advanced/2009/02/10/squashing-commits-with-rebase.html
// git rebase -i HEAD~2

// In git, pay attention to all the changes that would MODIFY the HISTORY: once you share, you cannot change
// because it would be hard for others to do merges and other stuff like that afterwards

// how to modify a commit message: git commit --amend -m "New commit message"

// how to place a tag
// git tag tag_name <commit>   //this is a lightweight tag, the easiest one

// OBJECTS IN GIT: blobs, trees, tags and commits

// can i checkout a tag which is not the tip of a branch, or will it tell me "detached head"? Certo, sarai in detached head,
// se quel commit non e' la tip di alcun branch

// Nell'altra  macchina, in cui ho un'altra installazione di petsc, 
// per usarla faccio semplicemente, nel repo della libreria (non dell'applicazione):
// git checkout petsc-3.3.6
// source myconf.sh

// Nella nuova macchina, rimango sul petsc-3.4.0, poi faccio
// export FM_MYMACHINE=dellnew
// source myconf.sh

// ============ PETSC configuration at runtime ==========
// E' possibile configurare a runtime tante caratteristiche di petsc perche' 
// c'e' PetscInitialize !! Quella fa la differenza... 
// Fara' qualcosa di simile al GetPot di LibmeshInit, che legge gli argomenti della riga di comando...
// la domanda pero' e': in PETSC, e in LIBMESH, DOVE VENGONO STORATI I DATI DELLA RIGA DI COMANDO?
// IN UNA VARIABILE CON GLOBAL SCOPE?!? 
// l'obiettivo e' sempre che sia accessibile a tutte le routine...

// Posso esimermi dall'usare LibmeshInit per le chiamate a Libmesh ?!?
// se SIA LIBMESH SIA FEMUS sono installate CON le PETSC,
// ci puo' essere il classico DOPPIO PetscInitialize
// IL PUNTO PRINCIPALE STA PROPRIO QUI

// ========== settare il COMUNICATORE MPI ===============

// ============= GIT FETCH, deeper ==============
// Se voglio pigliare le TAGS, faccio git fetch -t.
// DOPO NON HO BISOGNO DI FARE NESSUN git-merge per le tags!!!
// il git merge evidentemente riguarda i commits, o come si fanno chiamare...


// ========= GIT ===============
// the big difference with other revision control systems is that 
// git tracks content and not files.
// If two files have different names but the exact same content,
// their content will be represented by the exact same SHA1.


// ======= DEBUGGING STRATEGIES and FUNCTIONS ===========
// the calls to typical debug functions must be done in such a way that 
// the user can add or remove the calls  with an "if" in front, possibly at RUNTIME,
// but without losing PERFORMANCE.

// Well, actually RECOMPILING for DEBUGGING PURPOSES is QUITE NATURAL in the concept 
// of "dbg mode"... but here we want some "lighter" kind of debugging...

// ========= HOW TO REDIRECT X WINDOWS TO FILE (PETSC OUTPUT) ============

// ========= REDIRECTING OUTPUTS for every run ============
// Ora lo faccio con uno script ESTERNO, che sembra avere alcuni vantaggi essenziali.
// Il punto pero' e': COME POSSO REDIRIGERE gli output se ho un esecuzione PARALLELA?!
// Vorrei un file separato per ogni processo, con _NUMERO, pero' quel numero e' noto solo 
// DENTRO AL MAIN, DOPO CHE il main e' stato lanciato. DIREI CHE NON C'E' MODO DI SAPERLO PRIMA...
// A MENO CHE UNO NON RIESCA A SAPERE "AL VOLO", cioe' nel momento del lancio del main,
// il PID del processo (se non sbaglio ogni processo MPI ha un PID diverso)...
// Pero' sono quasi sicuro che il PID e' noto UN ATTIMO DOPO, APPENA DOPO il lancio dell'eseguibile...

//  Quindi, per ora, se voglio tenere il loop di generazione della cartella di output ESTERNO,
// devo accontentarmi di redirigere in UNICO FILE.
// La cosa interessante e' vedere che cosa succede se uso insieme la redirezione SERIALE e quella PARALLELA...
// Se fuori dal main dico di redirigere a fileA, ma dentro il main dico di redirigere a fileB, CHI VINCE?
// Immagino che l'ULTIMO SOVRASCRIVA prendendo il sopravvento su tutti.
// Dopodiche' il redirect esterno puo' fungere da "COLLETTORE DI SPORCIZIA", per far si' che tutto sia scritto da qualche parte...

// Per avere il pid di un programma IN ESECUZIONE, basta usare "pidof" e il nome dell'eseguibile.
// Pero' non posso fare "./main.out; pidof main.out" perche' dopo il punto e virgola l'eseguibile muore...

// Niente da fare, tutto avviene per forza UN ISTANTE DOPO il  lancio del main, quindi se voglio fare qualcosa 
// lo posso fare DENTRO AL MAIN.

// Per quanto riguarda il parallelo, se voglio fare dei run totalmente scollegati, cioe' con MPI_COMM_SELF, 
// allora posso generare la stringa di istante temporale ESTERNAMENTE, prima di chiamare mpiexec, 
// in modo che sia sicuramente UNICA e non ci possano essere dei ritardi tra un thread e l'altro. 
// Poi ALL'INTERNO DI CIASCUN PROCESSO, aggiungero' "_PROC_RANK" a quella stringa, E SOLO IN QUEL MOMENTO, 
// dentro il main cioe', CREERO' LA CARTELLA DEL PROCESSO.
// Questa non e' una cattiva idea, ma a quel punto e' piu' compatto fare tutto dentro al main...

// Ci potrebbe essere un'altra alternativa. 
// "Se tu non sai bene che cosa succedera', ma sai che quello che succedera' avra' almeno un certo confine, 
// allora puoi provare ad anticipare le mosse".
// Se per esempio so che lancero' con 4 processi, allora posso CREARE FIN DA SUBITO 4 CARTELLE 
// con un istante temporale e _0, _1, _2 e _3. Lo posso fare, certo!
// Dopodiche' il processo 0 andra' a sputare dentro 0, e cosi' via...
// No, aspetta, siamo da capo... IL PUNTO E' SAPER DIRE A MPIEXEC che OGNI PROCESSO SCRIVA IN UNA CARTELLA
// data da una certa stringa piu' il rank, del tipo "output_0".
// DEVO POTER SAPER DIRE QUESTA COSA ad MPIEXEC, NON PIU' TARDI!

// MA ASPETTA: siamo sicuri che quando redirigo lo std::cout in realta' non lo redirigo anche per TUTTE LE LIBRERIE LINKATE,
// anche se queste hanno printf o chesso' io?
// Magari se ALL'INTERNO DEL MAIN faccio un REDIRECT di std::cout e di std::cerr, forse i messaggi di errore 
// di PETSc o HDF5 finiscono TUTTI DENTRO IL FILE, e NON VANNO ALL'ESTERNO COME PENSO IO...
// Non escluderei questa prospettiva a priori...

// ====== GENERATION OF THE OUTPUT FOLDER STRING ==========
// Se decido di farla nel loop ESTERNO, devo trovare un modo per comunicarla al MAIN.
// Qualunque modo scelgo, e' importante che non ci sia il rischio di SOVRAPPOSIZIONI.
// Sto pensando in particolare al caso in cui io facessi un "mpiexec -n NP" con PROCESSI NON-COMUNICANTI
// (per fare questo uso il comunicatore MPI_COMM_SELF).
// Questi sono discorsi di THREAD SAFETY.
// Allora, bisogna confrontare i metodi di comunicazione al main, cioe':
// - lettura da file
// - lettura da shell environment
// - lettura da riga di comando

// E' chiaro che la lettura da riga di comando non va bene, perche' io passo UNA stringa del nome della cartella come argomento a "mpiexec -n 4"
// Anche la lettura da file puo' essere rischiosa perche' io avro' UN FILE UNICO che verra' letto da quattro processi.
// Per la variabile d'ambiente, dovrei essere sicuro che LE SHELL SIANO INDIPENDENTI,
// ma credo che la shell di "mpiexec -n NP" sia unica... non credo ci siano quattro sottoshell...

// Se le cose stanno cosi', la generazione dell'output folder puo' essere thread-safe solo ALL'INTERNO del thread, chiaramente...

// In questo modo pero' io ridirigerero' lo std::cout del mio programma, MA NON L'OUTPUT di HDF5, MPI, PETSC, eccetera... STA LI' IL PUNTO della questione.
// Quegli output sono ad es dei printf, non c'entrano niente con il mio std::cout!

// Giusto per avere un'idea, se voglio sapere qual e' la lista delle ENVIRONMENT VARIABLES di un dato processo,
// posso fare "cat /proc/<PID>/environ"


// ======== FUNCTION POINTERS vs. VIRTUAL FUNCTIONS ===========
// Virtual functions are good because they belong to the class, so they can access 
// the data of the class.
// Function pointers are MORE FLEXIBLE in the sense that you can change them AT RUNTIME,
// while virtual functions cannot be changed once compiled.


// ============ GOOD MESH GENERATOR =========
// - open source
// - possibilita' di fare SCRIPT, o avere il TRACING dei comandi grafici
// - supporto per HEX, e SURFACE MESHES
// - libmesh puo' leggerlo
// - COMSOL puo' leggerlo
// - che sia abbastanza diffuso
// - che dia le FLAG per i BOUNDARY e per le REGIONI DI VOLUME


// Salome
// OpenCascade
// GMSH
// Cubit (a pagamento)



// in gambit le ZONE che definisco devono essere delle facce, o edge, o volumi
// ok, if a mesh does not work, i'll do separate faces, do the mesh on each, and then UNITE THEM...
// when you UNITE faces, the mesh DISAPPEARS
// OK, that's NO PROBLEM, because a ZONE can be the UNION of FACES, so it's ok!
// AAA: devo ricordarmi di come sono ORIENTATI gli EDGE!!!
// Se mesho, partendo dal cerchio piu' piccolo, con 4 elementi sul boundary, e faccio map,
// mi dice "The intervals would force the mesh to cross itself"
// Invece con quattro elementi sul boundary funziona
// gli edge sono numerati doppi, perche' appartengono a due facce...
// quindi, quando faccio il mesh devo stare attento a questo...

// poi, non si possono cancellare gli EDGE se prima non si cancellano le FACE
// In giallo fa vedere gli edge che NON APPARTENGONO ad alcuna faccia.
// AAA: il tipo di solver lo devi selezionare PRIMA di decidere le ZONE di BOUNDARY!!!

// allora, la sola equazione della pressione non funziona, perche' non dovrebbe risolvermi 
// solo nella regione 8 e darmi zero in tutto il resto
// se cambio le condizioni al contorno della pressione, mi modifica il boundary, 
// ma nel resto e' come se NON RISOLVESSE... quindi credo che NON ASSEMBLI le regioni 5,6, e 7.
// Come fare per capirlo? Ad es. STAMPARE GRAFICAMENTE la matrice!

// i dont understand if i am SOLVING only on a subregion, or i am PRINTING only on a subregion.
// Since the bc's in the well are correct, i guess it is a problem about SOLVING, i.e.,
// about FILLING the MATRIX and RIGHT HAND SIDE...
// maybe the problem is about the ADAPTIVE REFINEMENT...
// if i remove the adaptive refinement, still i get nonzeros only in region 8... WHAT IS WRONG WITH THAT?!?
// let me just do one region with gambit, and that's it.

// ================
// Allora, con un mesh cosi' il flusso non viene proprio bellissimo...
// adesso faccio il mesh SENZA FRATTURA. sempre 4 zone, 4 cerchi stavolta, 500 200 50

// ================
// How can the numerators be THE SAME, if I changed the SCALE!!!


//========== TO SPECIFY BOUNDARY TYPES ==========
// You must use the ZONES part in GAMBIT.
// I don't remember if the flags for the volume parts must start at 5...


// ============ DIMENSIONAL ANALYSIS =================
// I need to check whether the numbers are correct or not
// GEOMETRICAL SETUP
//well radius: 15cm
//fracture longitudinal length: 100 m
//field radius: 500 m
//fracture diameter: 5mm

// our mesh has outer radius equal to ONE.
// but the field radius is 1/2 km, let us say 1 km,
// so there is a 10^6 factor in the mesh dimensionalization

// ======== PRODUCTIVITY INDEX =========
// What are its DIMENSIONS? It does not seem to be nondimensional!!!

// ============== continuo in questo computer per ora ========
// ok, ora vorrei vedere come cambia il productivity index con o senza frattura.
// posso farlo senza riaprire il mesh? vediamo.



// =====zypper see  file list of a package ===================
// // Fetch the files without installing:
// // zypper in --download-only <package>
// // 
// // Find the file:
// // find /var/cache/zypp -iname "package*rpm"
// // 
// // List the files in an uninstalled package:
// // rpm -qlp /var/cache/zypp/packages/<repo_alias>/suse/<arch>/<package-file-name>
// // 
// // That's --query --list --package <file>.

// ============
// I don't think it makes a lot of sense to do the simulation of the whole well...
//but, if you want to compute the P.I. you need to integrate over the whole volume

// ======= TO DETERMINE THE TYPE AT RUN-TIME
// can we find a way to determine the type of a variable 
// by reading a STRING at runtime which says "int" or "double" or "std::string" or whatever?
// well, that is not possible because, once a program is compiled, it has the type...
// unless maybe you start with a void* which may become a int* or a double* or a std::string*...
// since we use type templates, we could do something like T myvar = (T) runtime_map.get("pippo")...
// I think we could do a macro 
// How does PETSC do to keep track of the NAMES of the FUNCTIONS? They use some MACRO mechanism.
// Then we could do the same for the TYPES... mmmh i don't know, i think maybe the void* thing works

// ========= TYPE INFERENCE =========
// STATICALLY TYPED = a language that defines the type of a variable at compile time
// DYNAMICALLY TYPED = the type is known only at runtime

// ======= PHYSICAL DEPENDENCIES ==========
// Dependency on Es and Ks
// Where the dependency is the strongest, see for different LEVELS.
// It seems like the P.I. is a LINEAR FUNCTION of k... because the surface flux depends linearly on K

// From the results of the computations it seems like the productivity depends mostly on the WELL PERMEABILITY.
// It seems like it does not depend on BOTH ELASTIC CONSTANTS, but mainly on the PERMEABILITIES!!!
// So let us focus on the permeabilities first

// Sembra che la siano solo i KAPPA ad incidere sul P.I.
// Considerando invece gli E, l'unico da cui sembri dipendere il P.I.
// (seppure alla sesta cifra decimale) e' E_frac principalmente.

// La grossa variazione e' data da k_well. Quello che varia soprattutto con k_well 
// e' il NUMERATORE del P.I.  u = -k gradP
// gradP e' circa costante e quindi u varia con k, anche al bordo del well, e quindi il P.I.


// ====================== RUNTIME dependencies ================
//ok, i decided to make it runtime but then i want to change.
// do i have to do the loop outside, from shell, or can i do it inside here?
// Actually, you can do it all inside. 
// The only problem about doing from inside is to REDIRECT ALL the OUTPUT.
// If you do it from the shell, you could also redirect possible error messages,
// or other messages, that come FROM OTHER LIBRARIES than yours.
//In this way you would collect it all.
// Does the program need to know the string of time?
// yes, in order to create the directory (well, this can be better done outside)
// and in order to redirect NON JUST THE TERMINAL OUTPUTS, 
// but also the FILE PRINTING.
// now, how do I do that? I could write from shell to a file, and then read from file in the main.
// or, i could pass that shell variable to command line and read it from "argv" in the main.

//but, the point is that i want to have a different TIME for ANY VECTOR of PARAMETERS [E_1, E_2, k_1, k_2],
//so i must create the time as an INNER LOOP.
//So, in this way I would have to do the LOOPS on the vector of parameters IN THE SHELL AS WELL.
//Well, I may change the way of generating the FOLDER NAME in such a way that I can do the loops on the parameters INSIDE.
// I  can just do "initial time" outside + COUNTER OF THAT SPECIFIC RUN.
//No, I can't do that, otherwise I should generate LATER the folder in which I redirect all the INPUTS and OUTPUTS!!!
//So the loop on parameters must be either ALL INSIDE or ALL OUTSIDE.
// If the loop is ALL OUTSIDE, then you do not use C++ Math functions... on ther other hand,
//you are sure to collect ALL the OUTPUTS from ALL the CALLED LIBRARIES, which is good
//if you want to get possible errors from petsc,hdf5,libmesh or what else.

//So i will do like this: choose a vector of parameters, 
//create the TIME FOLDER,
//copy the template parameter file to the folder,
//modify the parameter file IN THE FOLDER with sed,
// PUT THE OUTPUT FOLDER NAME IN THE PARAMETER FILE?!? (it is a STRING, not a DOUBLE),
//copy the parameters file in it,
// run the main that reads from the parameter file IN THE FOLDER?!?, and send the OUTPUTS to the given folder, as well as 
//some QoI value to the 


//In any case, if the time is generated outside, i have to pass it to the MAIN.


// VARIOUS WAYS TO READ PARAMETERS:
// either COMPILE TIME
// or RUN TIME.

//If COMPILE TIME, you substitute with SED in some SOURCE FILE (for instance, with #define, or even in
// the variable assignment in the .C file), then you compile and run.
// So, you change the executable for every parameter.

//If RUN TIME, you have JUST ONE EXECUTABLE, and you read the parameters at run-time, in these ways:
// - either from command line
// - or from the shell environment
// - or from file.

// I decide to read from SHELL, so I just use "getenv" and i do not handle command line arguments
// it is clear that whether you read from SHELL, or from COMMAND LINE, you need to modify the EXECUTION COMMAND,
// and you do not need that for file. 

// The safest thing is TO CREATE A FOLDER FOR AN EXECUTION and DO EVERYTHING INSIDE THAT FOLDER.
// In this way you don't have any SUPERPOSITION of FILES

//Now the thing is this: how do I make these executions EMBARRASSINGLY PARALLEL?

//Of course you can make things ZERO INTRUSIVE by doing RUNS and then COPIES.
//But that means that you cannot do EMBARRASSINGLY PARALLEL RUNS


//The important thing is that I don't want to RECOMPILE at least the library, 
// and possibly also the main

//The point of the time is that no folder is ever superimposed
//The point about generating the time outside is that you "collect all the outputs, FROM ALL THE LINKED LIBRARIES",
//and not just by REDIRECTING the std::cout and std::cerr OF MY OWN LIBRARY.

//So the procedure is:
// create folder:
//launch executable by redirecting the output to that folder
// change the vector 


//The runtime map must be available to ALL THE SYSTEM. It must be accessible FROM ANYWHERE.
//Even if it concerns APPLICATION SPECIFIC NUMBERS, these are called in FUNCTIONS WHICH ARE CALLED BY THE LIBRARY,
//even though IMPLEMENTED BY THE USER (this happens through FUNCTION POINTERS, or through VIRTUALITY in INHERITANCE).
//Therefore those function must have the data AVAILABLE TO THEM!

// Adesso il problema, se lo voglio fare in parallelo, e' che l'istante potrebbe essere lo stesso,
// quindi devo fare il controllo all'atto della CREAZIONE DELLA DIRECTORY. Se c'e' gia' fermati e non andare avanti.



//How do i make an instantiation of a class be available
// to all the system?
// Either I make a global variable (can I do a global class?),
// or I do a Singleton (something that is instantiated only ONCE),
// or I pass it EXPLICITLY to all the classes/functions I need

//Will a change in model put us CLOSER or not to "one-day-available" experimental models?

//What is the point from a NUMERICAL POINT OF VIEW?
////Why should one use MULTILEVEL DOMAIN DECOMPOSITION rather than NOTHING?
// So now i am doing A PRIORI ADAPTIVE.This has nothing to deal with domain decomposition.
// I may not even mention it.
//The domain decomposition that i am using is NOT for the parallel computing thing...

// AM I USING IT NOW OR NOT ?!

//===============
// I must pass the RUNTIME MAP to the MULTIGRID CONSTRUCTOR, perche' ne puo' avere bisogno 
// GIA' in COSTRUZIONE.

//VANTAGGI del loop esterno di shell:
// puoi raccogliere TUTTI gli STDOUT e STDERR di TUTTE LE LIBRERIE in modo organico

// SVANTAGGI del loop esterno di shell
// le operazioni sui parametri vengono fatte con le funzioni aritmetiche di shell e non di C++;
// c'e' una variabile che devi passare dalla shell al main con getenv, il nome temporale della cartella,
// quindi se vuoi fare il main senza lo script devi prima fare almeno export di quel nome di cartella.

//REDO MESH WITH MORE REALISTIC MEASURES


//what is wrong? With a moving domain it seems like the equation 
// is not converging... why?!?

//remember one thing: ALL THE EQUATIONS MUST BE WRITTEN ON A MOVING DOMAIN!
// and also the routines that compute the integrals in the productivity index!

// How can we check that an element has not "superimposed itself"?
// It is not about CONVEXITY...
// The element must always be COUNTERCLOCKWISE like by convention.
// In questi movimenti puo' succedere di tutto, il quadrato puo' diventare come un triangolo e due lati coincidenti esterni... or whatever...
// Forse un check che si puo' far abbastanza facilmente e' il check di CONVESSITA'... prendi due lati consecutivi 
// orientati e controlli che l'angolo sia minore di 180...

//after one nonlinear iteration, nothing seems to change a lot, so it looks like one or very few nonlinear iterations are enough...

//Even though the displacement is negligible,
//this doesnt mean that it doesnt count, since there are DIFFERENT LENGTH SCALES involved

//The model is not well refined (models can be well refined only if they are matched against EXPERIMENTAL RESULTS.
// We do not have any refinement data available like that, so we see how we could use those models better...)


//variation of properties to study phenomena which would typically occur:
// - CHANGE OF PRODUCTIVITY INDEX
// - TANGENCY CONDITION: the two boundaries touch themselves,
// due to the SWELLING
// limit: VARIATIONAL STUDY, PARAMETER IDENTIFICATION PROBLEM (finite-dimensional parameters)

//NUMERICS: plot the RESIDUALS per iteration

//Of course, when studying the dependency on the material properties,
//even if the equations are THE SAME and only the COEFFICIENTS change,
//what is important are not the DIFFERENCES but the ABSOLUTE VALUES
//(typically, absolute values are described by RATIOS, not differences
// but it's also clear that the solution does not only depend on the RATIO,
// but also on the SINGLE ABSOLUTE VALUES)

// MOTIVATION

//MODELING
// NONLINEARITIES in FLUID


//NONLINEARITIES in SOLID

//NUMERICAL ANALYSIS


// GEOMETRICAL SETUP
//well radius: 15cm
//fracture longitudinal length: 100 m
//field radius: 500 m
//fracture diameter: 5mm

//it seems like the vector has NOT ZEROs
// an ARRAY of PetscScalar seems to be different from a PetscVectorM
// it seems like the Sol vector is not set to zero from the very beginning
//does the sol vector have to be set to zero only at the beginning or at the beginning of each iteration...
//only at the beginning i'd say
//Il vettore Sol contiene la condizione iniziale, che non e' detto che sia zero.
// Quindi,se vuoi simulare un'equazione sola SENZA ACCOPPIAMENTI, devi anche mettere ZERO nelle condizioi iniziail!

// laplacian of pressure with different coefficients:
// different permeabilities, much higher in fracture
// Figure with man-made fracturing
//http://www.ireservoir.com/case_fractured.html


// FRACTURE MODELING
// inside fracture we can expect even navier stokes
// in porous media the size of the fracture is very small
// Fracture NS
// solid: 

// Solid properties for rock and sand
// compute the productivity index
// what are tipical values for the Productivity Index?!



//the displacement in the fluid is given by the ale equation
//the displacement in the solid is given by the kinematic equation
//the velocity in the fluid is given by ale navier stokes
//the velocity in the solid is given by the equation of structural mechanics
//of course in both cases we need a time discretization for both displacement and velocity in both materials

// Mesh
// Dof,FE
// Assemble Eqns
//LoopSolver

//INDICES
// Index
// MGIndex
// VankaIndex

// what is the difference between level and gridn?!

//do the mesh with the ZONES for both BOUNDARY and CONTINUUM

//what is the DEFAULT for the BOUNDARY CONDITIONS?!
//do I have to put any REFERENCE VALUES?!

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////


//Here, if you want to avoid singleton pattern problems, you must pass the runtime map 
// to the MULTIGRID, which passes it to linearsolver, which passes it to lsyspde, which passes it to mesh
//therefore, the runtime map must belong to all objects. In particular it must belong to the mesh

//The SINGLETON pattern does not necessarily mean to have a PRIVATE constructor.
// The ones in libmesh for instance do not have private constructors...

// PROVIAMO A VEDERE SE QUESTO SINGLETON FUNZIONA. E' un barbatrucco tipo variabile globale
// che in genere non e' molto consigliato.
// Pero' vediamo se ci evita di fare i cambi delle interfacce.
// Oppure quella cosa tipo "anonymous namespace" funziona?

//I want to put STATIC data members to 

// static functions can only act on static data, clearly

//static stuff must be in a c++ file, while template stuff must be in the .h file, how do you concile the two...
//Also, static member does not mean CONST. static means that IT IS THE SAME FOR ALL INSTANTIATIONS.
// A singleton is something that is instantiated ONLY ONCE, so in that case it doesn't make a big difference
//to distinguish between what is instantiation-dependent and what is independent...

//in any case, either with singleton or changing all the interfaces,
//it means to CHANGE ALSO THE LIBRARY because the library must be aware of 
//the singleton which is now added for the specific application...
      
//anyway, let us try the singleton pattern for runtimemap.
//this class is made of just a MAP, a couple of strings, and few functions to read/write

//when i have a static datum, when and where do i fill it? do i have to fill it 
//before calling the constructor (the initializer list), or even after?
//can i declare static classes, not just C primitive datatypes?

//Ok, if you start from a class and you try to "MAKE ALL THE MEMBERS STATIC step by step",
//then it happens that the compiler complains, of course, because he wants you to initialize that member somehow.
// "A static data member can only be initialized at its definition".
// does it mean that the VALUE of a static data member must be known at COMPILE time, or not necessarily?!?
// maybe initializing does not necessarily mean put the value... maybe it can be changed later... let me try

//also, in-class initialization of non-const static member is prohibited, so i must do externally of the class definition.
// instead, "static const" member can be initialized in-class. 
// This is of course according to ISO C++.
// EXTERNALLY means in the .h file or even in a separate .C file?

// =============================
//What is the difference between static and const?

// these keywords already exist in C, but in C++ they assume further meanings, in particular "static"

// "static" means that each instantiation has the SAME VALUE of that data member.

// "const" means that it MUST be initialized. You cannot do const int pippo; you must put = something.
// "const" does NOT mean that it must be initialized AT COMPILE TIME! it just means that,
// once the variable is created in memory, a value for it must be set which will never be modified
// in the rest of the code!!! So basically, the variable is "created and filled once and for all" at the same time!
 
// const doesnt mean COMPILE TIME: it can be initialized at runtime, AND THEN REMAIN UNCHANGED!
// static doesnt mean COMPILE TIME: it can be filled at runtime, AND EVEN CHANGED!
    
// static doesnt mean const
// "static" basically means "GLOBAL variable"

// A "static" data member or member function basically means "IT DOESN'T BELONG TO A PARTICULAR CLASS INSTANCE".
// It MUST be initialized, but it doesn't mean IT CANNOT CHANGE.

// Also, notice that "static" concerns the VARIABLES, not the TYPES! 
// You cannot declare a "static class" or "static struct"


// ========== Declaration vs Definition: In Summary ===================

// We must have clear in mind the difference between DECLARATION and DEFINITION...
// well, actually according to the C++ standard, a definition IS a declaration...

// 
// A declaration provides basic attributes of a symbol: its type and its name. 
// 
// A definition provides all of the details of that symbol
// --if it's a function, what it does;
// if it's a class, what fields and methods it has;
// if it's a variable, where that variable is stored (you are telling the compiler where to create the STORAGE for that variable)
//   Often, the compiler only needs to have a declaration for something
//   in order to compile a file into an object file,
//   expecting that the linker can find the definition from another file.
//   If no source file ever defines a symbol, but it is declared,
//   you will get errors at link time complaining about undefined symbols.

//  When you use the EXTERN keyword, it means you are DECLARING a variable but you ARE NOT DEFINING IT:
// it is defined SOMEWHERE ELSE which will be found by the linker.
// When do I use extern? When i have GLOBAL variables?!?
//=======================


//So the goal is to make a GLOBAL CLASS
//i can make a static map, fill it in the main and then call it from elsewhere

//you cannot put the IMPLEMENTATIONS (DEFINITIONS) in HEADER files, you would have "MULTIPLE DEFINITION" in linking!!!
// the same applies to STATIC variables!!!

//now since a singleton can have a unique instantiation, how can i have it in several points of the code?
//maybe because i instantiate POINTERS...? No i don't think so

//for me it's ok even to have multiple instantiations but ALL static, so that i get THE SAME NUMBERS 
// in every part of the code (maybe if they are also CONST static i am sure that it is not changed...)

//if i instantiate a class WITHOUT new/delete, can i destruct the variable BEFORE the END OF SCOPE?

// ============ CONST ==========
// What does it mean the INITIALIZATION of a variable when it is CONST?!
// When a variable is declared as CONST, it must be initialized at the time it is allocated.
// Now, why
//const int pippo = "some value",
// but i can do
// const std::string giggio; 
// WITHOUT ANY "value"?!? 
  
//It means that SOME CONSTRUCTOR for the INITIALIZATION is called in any case?!
// So, what does it mean INITIALIZATION?! It means that the allocated memory 
// is accessed for the first time and "filled with something"? In that case it means that what is given
// by the default constructor will not change at all! Therefore, the only way to initialize 
// a const class is in the initializer list of the constructor!
// In fact, 
  
// non-integral type means a type which is not an integer

// ========= INITIALIZATION of a CONST data member ================
// When you have some "const" member datum, you need to initialize it
// If i want to initialize, i cannot do the IN-CLASS DEFINITION (.h file),
//   because i would have "multiple definitions", it is forbidden by the ISO standard. 
//   IN-CLASS is forbidden for ALL TYPES, const or non-const.
// So, I can do it OUTSIDE.
// I cannot do it outside like i do for STATIC members  WITH THE CLASS SCOPE ::.
// So, i can initialize in the INITIALIZER LIST in the CONSTRUCTORS.
// Notice that the compiler gives an error at the point of the CONSTRUCTORS! If you have several constructors,
// you need to initialize the CONST data members in the initializer list of EACH OF THEM!!!
// Also, in the initializer list you can initialize other class data which are NOT CONST, it's good to initialize more...

// ========== INITIALIZATION of a STATIC data member ====================
// if you declare a data member as STATIC, then you need to provide a DEFINITION of it.
// now i dont understand why if i have a STATIC INT or STATIC DOUBLE, i have to provide a DEFINITION
// in some cpp file, while this doesnt hold for STATIC STD::STRING... HOW COME?
// OOOOH HO CAPITO, non e' un problema di INT o DOUBLE vs STD::STRING!
// E' semplicemente che SE NON LE USI il LINKER NON TI DA' un "undefined reference"!!!
// E' per questo che uno POTREBBE NON INIZIALIZZARE le STATIC VARIABLES... 
// banalmente perche' non le usa, e quindi sarebbero totalmente inutili...
// magari il compiler e' cosi' furbo che non le considera di striscio...
// Siccome si spera che uno USI almeno una volta una static variable, allora DEVI DEFINIRLA,
// nel senso che devi dire al compilatore dove storarla in memoria!
// NOTA quindi che il compilatore e' piu' STRINGENTE con le const non usate rispetto alle static non usate:
// per le prime da' sempre un errore in compilazione, per le seconde non dice nulla.
// "Static" infatti e' qualcosa che tipicamente rimanda al LINKER.

// Now, NOTICE that "a static data member can only be initialized at its definition".
// Therefore, it cannot be initialized in the initializer list of the constructors!
// So, you must provide a DEFINITION, in a .cpp file, with the scope operator :: .

   
// Ora, per come e' fatto il Singleton, l'istanziazione e' unica. Ma il singleton e' GLOBAL SCOPE, o no?
// credo di no... Il singleton ha lo scope del main... io non sto usando una variabile con global scope
//Lo instanzio nel main, dopodiche' in ogni altro punto faccio getInstance e ottengo la stessa cosa.
//sarebbe meglio fare le cose come const, in modo che nessuno le possa toccare in seguito.

//La questione cruciale e': come faccio ad inizializzare una classe nell'inizializer list?
// la mia inizializzazione sarebbe la funzione read... come faccio a chiamare la funzione nell'inizializer?
// devo metterla nel COSTRUTTORE e sono a posto...

//Ok, so now the point is to initialize a static data member at runtime.
//Of course, if to initialize a datum i need some runtime operations in advance, then I cannot put
//all of them in front...
// of course i could just put variables with GLOBAL SCOPE, BUT I CANNOT EXECUTE ANY CODE "OUT OF SCOPE"!
// So I will initialize with something "fake" (default constructor, stuff like that), 
// but i won't make the variable const of course... this is not too good because every instantiation
//might in principle CHANGE the value, but since you make the variables PRIVATE then YOU KNOW HOW THE VALUE
// will be changed, because only you (the Class) are allowed to change the value.

//So, if you do cannot do static const, at least do STATIC, NON-CONST but at least PRIVATE.

//Now, this is a key point: if you are in a templated class, how do you INITIALIZE STATIC TEMPLATED data members?!?

// in conclusion... it's not that i have only one instance...i have more than one... each of them has a part of the data
// that are STATIC (instantiation independent), while other parts are instantiation-dependent (e.g., the counter)
// The point is that wherever i have to create an instance, i have to pass the constructor parameters, which might be 
// generated at runtime. So wherever i have i need to have those parameters available

// so, truly speaking, this is not exactly what i wanted to do... i wanted something COMPLETELY GLOBAL SCOPE...
// I think i could do this with some anonymous namespace... The problem is always this: if i put something
// under a GLOBAL SCOPE, but this thing needs the EXECUTION of some INSTRUCTIONS before its constructions,
//then of course i cannot do them outside of any function call...!
// I can make this trick!!! wrapping a DATUM into a FUNCTION that RETURNS THAT DATUM!!!

// ============ extern ======================
// "extern" is a way of communicating data between function without passing through the ARGUMENTS of the functions.
// "extern" in C is a variable defined OUTSIDE any function block
// To understand how external variables relate to the extern keyword,
// it is necessary to know the difference between defining and declaring a variable.
// When a variable is defined, the compiler allocates memory for that variable 
// and possibly also initializes its contents to some value.
// When a variable is declared, the compiler requires that the variable be defined elsewhere. 
// The declaration informs the compiler that a variable by that name and type exists,
// but the compiler need not allocate memory for it since it is allocated elsewhere.
// 
// The extern keyword means "declare without defining"

// extern is very DANGEROUS, because you define a variable in global scope, but in the various points 
// where you use it you don't actually know the status of that variable, you don't know if it was filled or not!!!

//E' QUESTO IL MODO PIU' PORCOSO PER EVITARE DI CAMBIARE GLI ARGOMENTI DELLE FUNZIONI.



// ========== WAYS TO COMMUNICATE between the EXECUTABLE and "THE EXTERNAL WORLD" ==============
// you have two alternatives: at COMPILE TIME, or at RUNTIME
// At RUNTIME you can:
// 1 - read from file
// 2 - read from command line
// 3 - read from shell environment

// At COMPILE TIME you could 
// 4 - use #define 

// 1 - you must do some PARSING routine; the file might be text, or binary; you might use some established file format
// 2 - these data are accessible only to main()
// 3 - these data can be read from anywhere, with getenv(...)

// 4 - you must do some PARSING that substitutes the strings with what is needed,
//     then you COMPILE and RUN. 
//     If you do runtime it is clear that you save the compile step all the time

//in all the 4 cases, it must be known a priori what is the variable TYPE:
// int, double, whatever,
// so there is some information that always needs to be known at COMPILE time.


// another tip: only constructors take member initializers...


// ====== PERFORMANCE of C++ GET CALLS ================
// what is the difference between accessing an object public datum 
// and getting it with a get function where the datum is private.

// so, what i mean is the difference between int foo =  myobject._myvariable; 
// and int foo = myobject.getmyvariable();
// this might make some difference if you are performing an operation a tremendous amount of times
// because you would save time
// to understand that, one should do a small program and analyze the result of the ASSEMBLER CODE.
// To understand the performances one should repeat the same operation an enormous amount of times
// and see if one highlights some "time increase"...

// also one should see what happens with or without "inline"
// i would say a priori that "inline" is the key (it belongs to both C and C++)

// ================== CONSTRUCTOR CALL ============
// if you declare a pointer to a class, then you call the constructor with the "new" operation.
// if you declare it as a class object, then you call a constructor, and then you cannot recall another constructor afterwards.


// ========== EXECUTE INSTRUCTIONS OUT OF ANY FUNCTION BLOCK ===========
// GLOBAL: when i declare a variable as global scope, so outside of any function block, the problem is that
// in order to manipulate it i need to EXECUTE INSTRUCTIONS and these must be executed within SOME FUNCTION BLOCK necessarily.
// So basically i am separating the declaration from the set of instructions needed to define the function.
// the constructor would be a way to put declaration and definition together...
// so, that is a way to EXECUTE INSTRUCTIONS OUT OF ANY FUNCTION BLOCK...
// another way would be 

// ========== ANONYMOUS NAMESPACE ===========
// it has FILE SCOPE. A global variable want to have GLOBAL scope, not just file scope, so that doesnt work.



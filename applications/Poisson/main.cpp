#include "MultiLevelProblem.hpp"
#include "MultiLevelMesh.hpp"
#include "TransientSystem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "FemTTUInit.hpp"
#include "SparseMatrix.hpp"
#include "VTKOutput.hpp"
#include "GMVOutput.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "SolvertypeEnum.hpp"
#include <json/json.h>
#include <json/value.h>

using std::cout;
using std::endl;

using namespace femus;


void AssembleMatrixResPoisson(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix);


// double InitVariableU(const double &x, const double &y, const double &z);


bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[], 
			  double &value, const int FaceName, const double time);

bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level);


static std::string
readInputTestFile( const char *path )
{
   FILE *file = fopen( path, "rb" );
   if ( !file )
      return std::string("");
   fseek( file, 0, SEEK_END );
   long size = ftell( file );
   fseek( file, 0, SEEK_SET );
   std::string text;
   char *buffer = new char[size+1];
   buffer[size] = 0;
   if ( fread( buffer, 1, size, file ) == (unsigned long)size )
      text = buffer;
   fclose( file );
   delete[] buffer;
   return text;
}


int main(int argc,char **args) {
  
    std::string path = "/home/simone/Software/Devel/femus/applications/Poisson/input/input.json";
//    Json::Features features;
//    bool parseOnly;
//    int exitCode = parseCommandLine( argc, argv, features, path, parseOnly );
//    if ( exitCode != 0 )
//    {
//       return exitCode;
//    }

   std::string input = readInputTestFile( path.c_str() );
   if ( input.empty() )
   {
      printf( "Failed to read input or empty input: %s\n", path.c_str() );
      return 3;
   }

  
  Json::Value root;   // will contains the root value after parsing.
  Json::Reader reader;
  bool parsingSuccessful = reader.parse(input, root );
  if ( !parsingSuccessful )
  {
     // report to the user the failure and their locations in the document.
     std::cout  << "Failed to parse configuration\n" << reader.getFormatedErrorMessages();
     return 1;
  }

  // Get the value of the member of root named 'encoding', return 'UTF-8' if there is no
  // such member.
  std::string encoding = root.get("encoding", "UTF-8" ).asString(); 
  // Get the value of the member of root named 'encoding', return a 'null' value if
  // there is no such member.
  const Json::Value plugins = root["plug-ins"];
  for ( int index = 0; index < plugins.size(); ++index )  // Iterates over the sequence elements.
    std::cout << plugins[index].asString() << std::endl;
     
  std::cout << root["indent"].get("length", 3).asInt() << std::endl;
  std::cout << root["indent"].get("use_space", true).asBool() << std::endl;
 
  std::string filename = root["mesh"].get("filename", "").asString();
  
  int numelemx;
  int numelemy;
  int numelemz;
  double xa, xb, ya, yb, za, zb;
  ElemType elemtype;
  
  bool isBox = root["mesh"].get("box", false).asBool();
  if(isBox) {
    numelemx = root["mesh"].get("box","").get("nx", 2).asUInt();
    numelemy = root["mesh"].get("box","").get("ny", 2).asUInt();
    numelemz = root["mesh"].get("box","").get("nz", 2).asUInt();
    xa = root["mesh"].get("box","").get("xa", 0.).asDouble();
    xb = root["mesh"].get("box","").get("xb", 1.).asDouble();
    ya = root["mesh"].get("box","").get("ya", 0.).asDouble();
    yb = root["mesh"].get("box","").get("yb", 1.).asDouble();
    za = root["mesh"].get("box","").get("za", 0.).asDouble();
    zb = root["mesh"].get("box","").get("zb", 0.).asDouble();
    std::string elemtypestr = root["mesh"].get("box","").get("elemtype", "Quad9").asString();
    if(elemtypestr == "Quad9") 
    {
      elemtype = QUAD9;
    }
  }
  
  std::string variableName = root["variable"].get("name", "Q").asString();
  std::cout << variableName << std::endl;
 
  std::string fe_order = root["variable"].get("fe_order", "biquadratic").asString();
  std::cout << fe_order << std::endl;
  
  unsigned int nlevels = root["mgsolver"].get("nlevels", 1).asInt();
  unsigned int npresmoothing = root["mgsolver"].get("npresmoothing", 1).asUInt();
  unsigned int npostmoothing = root["mgsolver"].get("npostsmoothing", 1).asUInt();
  std::string smoother_type  = root["mgsolver"].get("smoother_type", "gmres").asString();
  std::string mg_type        = root["mgsolver"].get("mg_type", "V_cycle").asString();
  unsigned int max_number_linear_iteration = root["mgsolver"].get("max_number_linear_iteration", 6).asUInt();
  double abs_conv_tol        = root["mgsolver"].get("abs_conv_tol", 1.e-09).asDouble();

  MgType mgtype;
  if (!strcmp("V_cycle",mg_type.c_str())) 
  {
    mgtype = V_CYCLE;
  }
  else if(!strcmp("F_cycle",mg_type.c_str())) 
  {
    mgtype = F_CYCLE;
  }
  else if(!strcmp("F_cycle",mg_type.c_str())) 
  {
    mgtype = F_CYCLE;
  }
  else {
    cout << "The selected MG cycle does not exist!" << endl;
    exit(1);
  }
  
  bool Vanka=0, Gmres=0, Asm=0;

  if( !strcmp("vanka",smoother_type.c_str()))          Vanka=1;
  else if( !strcmp("gmres",smoother_type.c_str()))     Gmres=1;
  else if( !strcmp("asm",smoother_type.c_str()))       Asm=1;
    
  if(Vanka+Gmres+Asm==0) {
    cout << "The selected MG smoother does not exist!" << endl;
    exit(1);
  }
  
  // end reading input from file
  //-----------------------------------------------------------------------------------------------
  
  /// Init Petsc-MPI communicator
  FemTTUInit mpinit(argc,args,MPI_COMM_WORLD);
  
  /// INIT MESH =================================  
  
  unsigned short nm,nr;
  nm=nlevels;
  std::cout<<"MULTIGRID levels: "<< nm << endl;

  nr=0;
  std::cout<<"MAX_REFINEMENT levels: " << nr << endl<< endl;
  
  int tmp=nm;  nm+=nr;  nr=tmp;
  
  //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;
  
  //Steadystate NonLinearMultiLevelProblem  
  MultiLevelMesh ml_msh;
  
  if(filename != "") 
  {
    ml_msh.ReadCoarseMesh(filename.c_str(),"seventh",Lref);
  }
  else
  {
    ml_msh.BuildBrickCoarseMesh(numelemx,numelemy,numelemz,xa,xb,ya,yb,za,zb,elemtype,"seventh");
  }
  ml_msh.RefineMesh(nm,nr,SetRefinementFlag);
  
  // ml_msh.EraseCoarseLevels(2);
  
  MultiLevelSolution ml_sol(&ml_msh);
  
  // generate solution vector
  ml_sol.AddSolution("Sol", fe_order.c_str());
 
  //Initialize (update Init(...) function)
  ml_sol.Initialize("Sol");
  
  //Set Boundary (update Dirichlet(...) function)
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  ml_sol.GenerateBdc("Sol");
  
  MultiLevelProblem ml_prob(&ml_msh,&ml_sol);
  
  // add fluid material
  Parameter parameter(Lref,Uref);
  
  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,0.001,1,"Newtonian",0.001,1.);
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;
  
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
   
  
  //BEGIN Poisson MultiLevel Problem
  std::cout << std::endl;
  std::cout << " *********** Poisson ************* " << std::endl;
    
  LinearImplicitSystem & system2 = ml_prob.add_system<LinearImplicitSystem>("Poisson");
  system2.AddSolutionToSytemPDE("Sol");
  
  // Set MG Options
  system2.AttachAssembleFunction(AssembleMatrixResPoisson);
  system2.SetMaxNumberOfLinearIterations(max_number_linear_iteration);
  system2.SetAbsoluteConvergenceTolerance(abs_conv_tol);  
  system2.SetMgType(mgtype);
  system2.SetNumberPreSmoothingStep(npresmoothing);
  system2.SetNumberPostSmoothingStep(npostmoothing);
   
  //Set Smoother Options
  if(Gmres) 		system2.SetMgSmoother(GMRES_SMOOTHER);
  else if(Asm) 		system2.SetMgSmoother(ASM_SMOOTHER);
  else if(Vanka)	system2.SetMgSmoother(VANKA_SMOOTHER);
  
  system2.init(); 
  //common smoother option
  system2.SetSolverFineGrids(GMRES); 
  system2.SetTolerances(1.e-12,1.e-20,1.e+50,4);
  system2.SetPreconditionerFineGrids(ILU_PRECOND);
  //for Vanka and ASM smoothers
  system2.ClearVariablesToBeSolved();
  system2.AddVariableToBeSolved("All");
  system2.SetNumberOfSchurVariables(0);
  system2.SetElementBlockNumber(4);                
  //for Gmres smoother
  system2.SetDirichletBCsHandling(PENALTY); 
  
  // Solve Temperature system
  ml_prob.get_system("Poisson").solve();
  //END Temperature Multilevel Problem
    
  /// Print all solutions
  std::vector<std::string> print_vars;
  print_vars.push_back("Sol");
      
  VTKOutput vtkio(ml_sol);
  vtkio.write_system_solutions("biquadratic",print_vars);
  
  //GMVOutput gmvio(ml_sol);
  //gmvio.write_system_solutions("biquadratic",print_vars);
  
  //Destroy all the new systems
  ml_prob.clear();
  
  //delete [] infile;
  return 0;
}

//-----------------------------------------------------------------------------------------------------------------

bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber, const int &level) {
  bool refine=0;
  // refinemenet based on Elemen Group Number
  if(ElemGroupNumber==5 ) {
    refine=1;
  }
  if(ElemGroupNumber==6 && level<2) {
    refine=1;
  }
  if(ElemGroupNumber==7 ) {
    refine=0;
  }

  return refine;
}

//--------------------------------------------------------------------------------------------------------------

// double InitVariableU(const double &x, const double &y, const double &z) { 
//    double um = 0.2;
//    double  value=1.5*um*(4.0/(0.1681))*y*(0.41-y); 
//    return value;
// }

//-------------------------------------------------------------------------------------------------------------------

bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[], 
			  double &value, const int FaceName, const double time){
  bool test=1; //Dirichlet
  value=0.;

  if(!strcmp(name,"Sol")) {
    if(1==FaceName){        // bottom face
      test=1;
      value=1;
    }  
    else if(2==FaceName ){  // right face
      test=0;
      value=0;
    }
    else if(3==FaceName ){  // top face
      test=0;
      value=10;
    }
    else if(4==FaceName ){  // left face
      test=1;
      value=1.;
    }
  }  
  
  return test;
}

// //------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------
void AssembleMatrixResPoisson(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix){
  
  //pointers and references
  Solution*      mysolution	       = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearImplicitSystem& mylin_impl_sys = ml_prob.get_system<LinearImplicitSystem>("Poisson");
  LinearEquationSolver*  mylsyspde     = mylin_impl_sys._LinSolver[level];   
  mesh*          mymsh		       = ml_prob._ml_msh->GetLevel(level);
  elem*          myel		       = mymsh->el;
  SparseMatrix*  myKK		       = mylsyspde->_KK;
  NumericVector* myRES		       = mylsyspde->_RES;
  MultiLevelSolution* ml_sol           = ml_prob._ml_sol;
  
  //data
  const unsigned	dim	= mymsh->GetDimension();
  unsigned 		nel	= mymsh->GetElementNumber();
  unsigned 		igrid	= mymsh->GetGridNumber();
  unsigned 		iproc	= mymsh->GetProcID();
  double		IPe	= 1./(ml_prob.parameters.get<Fluid>("Fluid").get_Peclet_number());  
  
  //solution variable
  unsigned SolIndex;  
  unsigned SolPdeIndex;
  SolIndex=ml_sol->GetIndex("Sol");
  SolPdeIndex=mylin_impl_sys.GetSolPdeIndex("Sol");
  //solution order
  unsigned order_ind = ml_sol->GetSolutionType(SolIndex);
  unsigned end_ind   = mymsh->GetEndIndex(order_ind);
  
  //coordinates
  vector< vector < double> > coordinates(dim); 
  
  // declare 
  vector< int > metis_node;
  vector< int > KK_dof;
  vector <double> phi;
  vector <double> gradphi;  
  double weight;
  vector< double > F;
  vector< double > B;
 
  // reserve 
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  metis_node.reserve(max_size);
  KK_dof.reserve(max_size);
  for(int i=0;i<dim;i++) 
    coordinates[i].reserve(max_size);
  phi.reserve(max_size);
  gradphi.reserve(max_size*dim);
  F.reserve(max_size);
  B.reserve(max_size*max_size);
  
  // Set to zeto all the entries of the Global Matrix
  if(assembe_matrix) myKK->zero();
  
  // *** element loop ***
 
  for (int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

    unsigned kel = mymsh->IS_Mts2Gmt_elem[iel];
    short unsigned kelt=myel->GetElementType(kel);
    unsigned nve=myel->GetElementDofNumber(kel,end_ind);
    
    // resize
    metis_node.resize(nve);
    KK_dof.resize(nve);
    phi.resize(nve);
    gradphi.resize(nve*dim);
    for(int i=0;i<dim;i++){
      coordinates[i].resize(nve);
    }
    
    // set to zero all the entries of the FE matrices
    F.resize(nve);
    memset(&F[0],0,nve*sizeof(double));
    if(assembe_matrix){
      B.resize(nve*nve);
      memset(&B[0],0,nve*nve*sizeof(double));
    }
    
    // get local to global mappings
    for( unsigned i=0;i<nve;i++){
      unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;
      unsigned inode_coord_metis=mymsh->GetMetisDof(inode,2);
      metis_node[i]=mymsh->GetMetisDof(inode,order_ind);
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coordinates[ivar][i]=(*mymsh->_coordinate->_Sol[ivar])(inode_coord_metis);
      }
      KK_dof[i]=mylsyspde->GetKKDof(SolIndex,SolPdeIndex,inode);
    }
        
    if(igrid==gridn || !myel->GetRefinedElementIndex(kel)) {
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < ml_prob._ml_msh->_type_elem[kelt][order_ind]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	(ml_prob._ml_msh->_type_elem[kelt][order_ind]->*(ml_prob._ml_msh->_type_elem[kelt][order_ind])->Jacobian_ptr)(coordinates,ig,weight,phi,gradphi);
	//current solution
	double SolT=0;
	vector < double > gradSolT(dim,0.);
	for(unsigned ivar=0; ivar<dim; ivar++){
	  gradSolT[ivar]=0; 
	}
  
	unsigned SolType=ml_sol->GetSolutionType("Sol");
	for(unsigned i=0; i<nve; i++) {
	  double soli = (*mysolution->_Sol[SolIndex])(metis_node[i]);
	  SolT+=phi[i]*soli;
	  for(unsigned ivar2=0; ivar2<dim; ivar2++) gradSolT[ivar2] += gradphi[i*dim+ivar2]*soli; 
	}
	// *** phi_i loop ***
	for(unsigned i=0; i<nve; i++){
	  //BEGIN RESIDUALS A block ===========================
	  double Adv_rhs=0;
	  double Lap_rhs=0;
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    Lap_rhs += gradphi[i*dim+ivar]*gradSolT[ivar];
	  }
	  
	  F[i]+= (-IPe*Lap_rhs + 1.*phi[i] )*weight;     
	  
	  //END RESIDUALS A block ===========================
	  if(assembe_matrix){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve; j++) {
	      double Lap=0;
	      double Adv1=0;
	      for(unsigned ivar=0; ivar<dim; ivar++) {
	      // Laplacian
		Lap  += gradphi[i*dim+ivar]*gradphi[j*dim+ivar]*weight;
	      }
	      B[i*nve+j] += IPe*Lap;
	    } // end phij loop
	  } // end phii loop
	} // endif assembe_matrix
      } // end gauss point loop
    } // endif single element not refined or fine grid loop
    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the global Matrix/vector
      
    myRES->add_vector_blocked(F,KK_dof);
    if(assembe_matrix) myKK->add_matrix_blocked(B,KK_dof,KK_dof);  
  } //end list of elements loop for each subdomain
    
  myRES->close();
  if(assembe_matrix) myKK->close();
  
   // ***************** END ASSEMBLY *******************
  
}


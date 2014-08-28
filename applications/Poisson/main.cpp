#include "MultiLevelProblem.hpp"
#include "MultiLevelMesh.hpp"
#include "TransientSystem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "FemTTUInit.hpp"
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "XDMFWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "SolvertypeEnum.hpp"
#include "FElemTypeEnum.hpp"
#include <json/json.h>
#include <json/value.h>
#include "ParsedFunction.hpp"

using std::cout;
using std::endl;

using namespace femus;


void AssemblePoissonMatrixandRhs(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix);


// double InitVariableU(const double &x, const double &y, const double &z);


// bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[],
//                           double &value, const int FaceName, const double time);

// bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level);



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

static void show_usage()
{
    std::cout << "Use --inputfile variable to set the input file" << endl;
    std::cout << "e.g.: ./Poisson --inputfile ./input/input.json" << std::endl;
}

ParsedFunction fpsource;


int main(int argc,char **argv) {

    std::string path;

    if(argc < 2)
    {
        std::cout << "Error: no input file specified!" << std::endl;
        show_usage();
        return 1;
    }

    for (int count = 1; count < argc; ++count)
    {
        std::string arg = argv[count];

        if ((arg == "-h") || (arg == "--help")) {
            show_usage();
            return 0;
        }
        else if ((arg == "-i") || (arg == "--inputfile"))
        {
            if (count + 1 < argc) {
                path = argv[++count];
            }
            else
            {
                std::cerr << "--input file option requires one argument." << std::endl;
                return 1;
            }
        }
    }

    // start reading input from file
    //-----------------------------------------------------------------------------------------------

    std::string input = readInputTestFile( path.c_str() );
    if ( input.empty() )
    {
        printf( "Failed to read input or empty input: %s\n", path.c_str() );
        return 1;
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
        std::string elemtypestr = root["mesh"].get("box","").get("elem_type", "Quad9").asString();
        if(elemtypestr == "Quad9")
        {
            elemtype = QUAD9;
        }
        else if(elemtypestr == "Tri6")
        {
            elemtype = TRI6;
        }
        else if(elemtypestr == "Edge3")
        {
            elemtype = EDGE3;
        }
        else if(elemtypestr == "Hex27")
        {
            elemtype = HEX27;
        }
        else
	{
	    elemtype = INVALID_ELEM;
	}
        
    }

    std::string variableName = root["variable"].get("name", "Q").asString();

    FEOrder fe_order;
    std::string fe_order_str = root["variable"].get("fe_order", "first").asString();
    if (!strcmp(fe_order_str.c_str(),"first"))
    {
        fe_order = FIRST;
    }
    else if (!strcmp(fe_order_str.c_str(),"serendipity"))
    {
        fe_order = SERENDIPITY;
    }
    else if (!strcmp(fe_order_str.c_str(),"second"))
    {
      fe_order = SECOND;
    }
    else
    {
      std::cerr << " Error: Lagrange finite element order not supported!" << std::endl;
      exit(1);
    }
    
    
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

    
    // reading function
     std::string function;
     function = root["variable"].get("func_source", "0.").asString();
     std::string variables = "x";
     variables += ",y";
     variables += ",z";
     variables += ",t";
     
#ifdef HAVE_FPARSER       
       fpsource.SetExpression(function);
       fpsource.SetIndependentVariables(variables);
       fpsource.Parse();
#endif
       
       std::vector<std::string> facenamearray;
       std::vector<bool> ishomogeneousarray;
       std::vector<ParsedFunction> parsedfunctionarray;
       std::vector<BDCType> bdctypearray;
       
       const Json::Value boundary_conditions = root["boundary_conditions"];
       for(unsigned int index=0; index<boundary_conditions.size(); ++index) {
	 
	 std::string facename = boundary_conditions[index].get("facename","to").asString();
         facenamearray.push_back(facename);
	 
	 std::string bdctypestr = boundary_conditions[index].get("bdc_type","to").asString(); 
         BDCType bdctype = DIRICHLET;
	 if (bdctypestr.compare("dirichlet") == 0) {
	     bdctype = DIRICHLET;
	 }
	 else if(bdctypestr.compare("neumann") == 0) {
	      bdctype = NEUMANN;
	 }
         else {
	       std::cout << "Boundary condition not implemented! Default one (Dirichlet) is set." << std::endl;
         }
         bdctypearray.push_back(bdctype);

	 bool ishomo = boundary_conditions[index].get("is_homogeneous",true).asBool(); 
	 ishomogeneousarray.push_back(ishomo);
	 
	 std::string bdcfuncstr = boundary_conditions[index].get("bdc_func","0.").asString();
	 ParsedFunction pfunc(bdcfuncstr, "x,y,z,t");
	 parsedfunctionarray.push_back(pfunc);
       }
       
    //---------------------------------------------------------------------------

    
    /// Init Petsc-MPI communicator
    FemTTUInit mpinit(argc,argv,MPI_COMM_WORLD);

    /// INIT MESH =================================

    unsigned short nm,nr;
    nm=nlevels;

    nr=0;

    int tmp=nm;
    nm+=nr;
    nr=tmp;

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
    ml_msh.RefineMesh(nm,nr, NULL);
    
    ml_msh.print_info();

    MultiLevelSolution ml_sol(&ml_msh);

    // generate solution vector
    ml_sol.AddSolution("Sol", LAGRANGE, fe_order);

    //Initialize (update Init(...) function)
    ml_sol.Initialize("Sol");

    //Set Boundary (update Dirichlet(...) function)
    ml_sol.InitializeBdc();
    
//     ml_sol.SetBoundaryCondition("Sol","right", NEUMANN, false, false, &bdcfunc);
//     ml_sol.SetBoundaryCondition("Sol","top", NEUMANN);
    
    for(int i=0; i<boundary_conditions.size(); ++i) {
      ml_sol.SetBoundaryCondition("Sol",facenamearray[i],bdctypearray[i],ishomogeneousarray[i],false,&parsedfunctionarray[i]);
    }
    
    ml_sol.GenerateBdc();
    
    
    //ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
    //ml_sol.GenerateBdc("Sol");
    
    

    MultiLevelProblem ml_prob(&ml_msh,&ml_sol);
    
    
//     ml_prob.parameters.set<func>("func_source") = fpsource;

    // add fluid material
    Parameter parameter(Lref,Uref);

    //BEGIN Poisson MultiLevel Problem
    std::cout << std::endl;
    std::cout << " PDE problem to solve: Poisson " << std::endl;

    LinearImplicitSystem & system2 = ml_prob.add_system<LinearImplicitSystem>("Poisson");
    system2.AddSolutionToSytemPDE("Sol");

    // Set MG Options
    system2.AttachAssembleFunction(AssemblePoissonMatrixandRhs);
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
    if(Asm || Vanka)
    {
      system2.SetNumberOfSchurVariables(0);
      system2.SetElementBlockNumber(4);
    }
    //for Gmres smoother
    system2.SetDirichletBCsHandling(PENALTY);

    // Solve Temperature system
    ml_prob.get_system("Poisson").solve();
    //END Temperature Multilevel Problem

    /// Print all solutions
    std::vector<std::string> print_vars;
    print_vars.push_back("Sol");

    VTKWriter vtkio(ml_sol);
    vtkio.write_system_solutions("biquadratic",print_vars);

    GMVWriter gmvio(ml_sol);
    gmvio.write_system_solutions("biquadratic",print_vars);
// 
//     XDMFWriter xdmfio(ml_sol);
//     xdmfio.write_system_solutions("biquadratic",print_vars);
    
    //Destroy all the new systems
    ml_prob.clear();

    return 0;
}

//-----------------------------------------------------------------------------------------------------------------

// bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber, const int &level) {
//     bool refine=0;
//     // refinemenet based on Elemen Group Number
//     if(ElemGroupNumber==5 ) {
//         refine=1;
//     }
//     if(ElemGroupNumber==6 && level<2) {
//         refine=1;
//     }
//     if(ElemGroupNumber==7 ) {
//         refine=0;
//     }
// 
//     return refine;
// }

//--------------------------------------------------------------------------------------------------------------

// double InitVariableU(const double &x, const double &y, const double &z) {
//    double um = 0.2;
//    double  value=1.5*um*(4.0/(0.1681))*y*(0.41-y);
//    return value;
// }

// 2D benchmark Full Dirichlet solution = sin^2(pi*x)*sin^2(pi*y)
// double Source(const double* xyz) {
//     const double pi = 3.1415926535897932;
//     double value = -2.*pi*pi*( cos(2.*pi*xyz[0])*sin(pi*xyz[1])*sin(pi*xyz[1]) + sin(pi*xyz[0])*sin(pi*xyz[0])*cos(2.*pi*xyz[1]) ) ;
//     return value;
// }

// 3D benchmark sin^2(pi*x)*sin^2(pi*y)*sin^2(pi*z)
//  double Source(const double* xyzt) {
//      const double pi = 3.1415926535897932;
//      double value = -2.*pi*pi*(  cos(2.*pi*xyzt[0])*sin(pi*xyzt[1])*sin(pi*xyzt[1])*sin(pi*xyzt[2])*sin(pi*xyzt[2]) 
//                                + sin(pi*xyzt[0])*sin(pi*xyzt[0])*cos(2.*pi*xyzt[1])*sin(pi*xyzt[2])*sin(pi*xyzt[2]) 
// 			       + sin(pi*xyzt[0])*sin(pi*xyzt[0])*sin(pi*xyzt[1])*sin(pi*xyzt[1])*cos(2.*pi*xyzt[2]) );
//      return value;
// }

// 1D benchmark 
// double Source(const double* xyz) {
//     double value = 1. - 2.*xyz[0]*xyz[0] ;
//     return value;
// }

//-------------------------------------------------------------------------------------------------------------------

// bool SetBoundaryCondition(const double &x, const double &y, const double &z,const char name[],
//                           double &value, const int FaceName, const double time) {
//     bool test=1; //Dirichlet
//     value=0.;
//     
//     std::cout << FaceName << std:: endl;
// 
//     if(!strcmp(name,"Sol")) {
//         if(1==FaceName) {       // bottom face
//             test=1;
//             value=0;
//         }
//         else if(2==FaceName ) { // right face
//              test=0;
//              value=0.;
//         }
//         else if(3==FaceName ) { // top face
//             test=1;
//             value=0.;
//         }
//         else if(4==FaceName ) { // left face
//             test=1;
//             value=0.;
//         }
//     }
// 
//     return test;
// }

// //------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------
void AssemblePoissonMatrixandRhs(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix) {

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
    unsigned 		iproc	= mymsh->processor_id();

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
    double src_term = 0.;

    // reserve
    const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
    metis_node.reserve(max_size);
    KK_dof.reserve(max_size);
    for(int i=0; i<dim; i++)
        coordinates[i].reserve(max_size);
    phi.reserve(max_size);
    gradphi.reserve(max_size*dim);
    F.reserve(max_size);
    B.reserve(max_size*max_size);

    // Set to zeto all the entries of the Global Matrix
    if(assembe_matrix) 
      myKK->zero();

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
        for(int i=0; i<dim; i++) {
            coordinates[i].resize(nve);
        }

        // set to zero all the entries of the FE matrices
        F.resize(nve);
        memset(&F[0],0,nve*sizeof(double));
        if(assembe_matrix) {
            B.resize(nve*nve);
            memset(&B[0],0,nve*nve*sizeof(double));
        }

        // get local to global mappings
        for( unsigned i=0; i<nve; i++) {
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
            for(unsigned ig=0; ig < ml_prob._ml_msh->_type_elem[kelt][order_ind]->GetGaussPointNumber(); ig++) {
                // *** get Jacobian and test function and test function derivatives ***
                (ml_prob._ml_msh->_type_elem[kelt][order_ind]->*(ml_prob._ml_msh->_type_elem[kelt][order_ind])->Jacobian_ptr)(coordinates,ig,weight,phi,gradphi);
                //current solution
                double SolT=0;
                vector < double > gradSolT(dim,0.);
                for(unsigned ivar=0; ivar<dim; ivar++) {
                    gradSolT[ivar]=0;
                }

                double xyzt[4] = {0.,0.,0.,0.};
                unsigned SolType=ml_sol->GetSolutionType("Sol");
                for(unsigned i=0; i<nve; i++) {
                    double soli = (*mysolution->_Sol[SolIndex])(metis_node[i]);
		    for(unsigned ivar=0; ivar<dim; ivar++) {
		      xyzt[ivar] += coordinates[ivar][i]*phi[i]; 
		    }
                    SolT+=phi[i]*soli;
                    for(unsigned ivar2=0; ivar2<dim; ivar2++) gradSolT[ivar2] += gradphi[i*dim+ivar2]*soli;
                }
                  
                // *** phi_i loop ***
                for(unsigned i=0; i<nve; i++) {
                    //BEGIN RESIDUALS A block ===========================
                    double Adv_rhs=0;
                    double Lap_rhs=0;
                    for(unsigned ivar=0; ivar<dim; ivar++) {
                        Lap_rhs += gradphi[i*dim+ivar]*gradSolT[ivar];
                    }
                    
             
                   //src_term = Source(xyzt);
#ifdef HAVE_FPARSER
                   // src_term = fpsource.Eval(xyzt);
                   src_term = fpsource(xyzt);
#endif
                    
                    F[i]+= (-Lap_rhs + src_term*phi[i] )*weight;
		    
                    //END RESIDUALS A block ===========================
                    if(assembe_matrix) {
                        // *** phi_j loop ***
                        for(unsigned j=0; j<nve; j++) {
                            double Lap=0;
                            double Adv1=0;
                            for(unsigned ivar=0; ivar<dim; ivar++) {
                                // Laplacian
                                Lap  += gradphi[i*dim+ivar]*gradphi[j*dim+ivar]*weight;
                            }
                            B[i*nve+j] += Lap;
                        } // end phij loop
                    } // end phii loop
                } // endif assembe_matrix
            } // end gauss point loop

            //number of faces for each type of element
            unsigned nfaces = myel->GetElementFaceNumber(kel);

            // loop on faces
            for(unsigned jface=0; jface<nfaces; jface++) {

                // look for boundary faces
                if(myel->GetFaceElementIndex(kel,jface)<0) {

                    double tau = 0.;
                    double time = 0.;
                    int node[27];
                    double vx[3][27];
                    double phi[27],gradphi[27][3],Weight;
                    double normal[3];
		    double xyzt[4] = {0.,0.,0.,0.};
		    ParsedFunction* bdcfunc;
		    
		    
		    unsigned int face = -(mymsh->el->GetFaceElementIndex(kel,jface)+1) - 1;
		    
		    if(ml_sol->GetBoundaryCondition("Sol",face) == NEUMANN && !ml_sol->Ishomogeneous("Sol",face)) {

		        bdcfunc = (ParsedFunction* )(ml_sol->GetBdcFunction("Sol", face));
			unsigned nve = mymsh->el->GetElementFaceDofNumber(kel,jface,order_ind);
                        const unsigned FELT[6][2]= {{3,3},{4,4},{3,4},{5,5},{5,5},{6,6}};
                        unsigned felt = FELT[kelt][jface<mymsh->el->GetElementFaceNumber(kel,0)];

                        for(unsigned i=0; i<nve; i++) {
                            unsigned inode=mymsh->el->GetFaceVertexIndex(kel,jface,i)-1u;
                            node[i] = inode + mylsyspde->KKIndex[0];
                            unsigned inode_Metis=mymsh->GetMetisDof(inode,2);

                            vx[0][i]=(*mymsh->_coordinate->_Sol[0])(inode_Metis);
                            vx[1][i]=(*mymsh->_coordinate->_Sol[1])(inode_Metis);
                            vx[2][i]=(*mymsh->_coordinate->_Sol[2])(inode_Metis);
                        }

                        if(felt != 6) 
			{
                          for(unsigned igs=0; igs < ml_prob._ml_msh->_type_elem[felt][order_ind]->GetGaussPointNumber(); igs++) {
                            (ml_prob._ml_msh->_type_elem[felt][order_ind]->*ml_prob._ml_msh->_type_elem[felt][order_ind]->Jacobian_sur_ptr)(vx,igs,Weight,phi,gradphi,normal);

                            for(unsigned i=0; i<nve; i++) {
			        xyzt[0] = vx[0][i];
				xyzt[1] = vx[1][i];
				xyzt[2] = vx[2][i];
				xyzt[3] = time;

 				double tau = (*bdcfunc)(xyzt);
                                //SetBoundaryCondition(vx[0][i],vx[1][i],vx[2][i],"Sol",tau,-(mymsh->el->GetFaceElementIndex(kel,jface)+1),time);
                                double bdintegral = phi[i]*tau*Weight;
                                unsigned int ilocalnode = mymsh->el->GetLocalFaceVertexIndex(kel, jface, i);
                                F[ilocalnode] += bdintegral;
                            }
                          }
			}
			else // 1D : the side elems are points and does not still exist the point elem
			{
			  xyzt[0] = vx[0][0];
		          xyzt[1] = vx[1][0];
			  xyzt[2] = vx[2][0];
			  xyzt[3] = time;
				
			  //SetBoundaryCondition(vx[0][0],vx[1][0],vx[2][0],"Sol",tau,-(mymsh->el->GetFaceElementIndex(kel,jface)+1),time);
			  double bdintegral = (*bdcfunc)(xyzt);;
			  unsigned int ilocalnode = mymsh->el->GetLocalFaceVertexIndex(kel, jface, 0);
                          F[ilocalnode] += bdintegral;
			}
                    }
                }
            }
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



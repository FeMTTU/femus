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
#include "ParsedFunction.hpp"
#include "InputParser.hpp"

using namespace femus;


void AssemblePoissonMatrixandRhs(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix);

void show_usage()
{
    std::cout << "Use --inputfile variable to set the input file" << std::endl;
    std::cout << "e.g.: ./Poisson --inputfile $FEMUS_DIR/applications/Poisson/input/input.json" << std::endl;
}

ParsedFunction fpsource;


int main(int argc,char **argv) {

    std::string path;

    if(argc < 2)
    {
        std::cout << argv[0] << ": You must specify the input file" << std::endl;
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
//         else {
// 	  std::cerr << argv[count] << " : command line argument not recognized" << std::endl;
// 	  show_usage();
// 	  return 1;
// 	}
    }

    /// Init Petsc-MPI communicator
    FemTTUInit mpinit(argc,argv,MPI_COMM_WORLD);
    
    // input parser pointer
    std::auto_ptr<InputParser> inputparser = InputParser::build(path);

    /// INIT MESH =================================

    unsigned short nm,nr;
    unsigned int nlevels = inputparser->getValue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.nlevels",1);   
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

    if(inputparser->isTrue("multilevel_mesh.first.type","filename"))
    {
      std::string filename = inputparser->getValue("multilevel_mesh.first.type.box.filename", "./input/input.neu"); 
      ml_msh.ReadCoarseMesh(filename.c_str(),"seventh",Lref);
    }
    else if(inputparser->isTrue("multilevel_mesh.first.type","box"))
    {
      int numelemx = inputparser->getValue("multilevel_mesh.first.type.box.nx", 2);
      int numelemy = inputparser->getValue("multilevel_mesh.first.type.box.ny", 2);
      int numelemz = inputparser->getValue("multilevel_mesh.first.type.box.nz", 0);
      double xa = inputparser->getValue("multilevel_mesh.first.type.box.xa", 0.);
      double xb = inputparser->getValue("multilevel_mesh.first.type.box.xb", 1.);
      double ya = inputparser->getValue("multilevel_mesh.first.type.box.ya", 0.);
      double yb = inputparser->getValue("multilevel_mesh.first.type.box.yb", 1.);
      double za = inputparser->getValue("multilevel_mesh.first.type.box.za", 0.);
      double zb = inputparser->getValue("multilevel_mesh.first.type.box.zb", 0.);
      ElemType elemtype = inputparser->getValue("multilevel_mesh.first.type.box.elem_type", QUAD9);
      ml_msh.BuildBrickCoarseMesh(numelemx,numelemy,numelemz,xa,xb,ya,yb,za,zb,elemtype,"seventh");
    }
    else
    {
        std::cerr << "Error: no input mesh specified. Please check to have added the keyword mesh in the input json file! " << std::endl;
        return 1;
    }
    ml_msh.RefineMesh(nm,nr, NULL);
    
    ml_msh.print_info();

    MultiLevelSolution ml_sol(&ml_msh);

    // generate solution vector
    FEOrder fe_order = inputparser->getValue("multilevel_solution.multilevel_mesh.first.variable.first.fe_order",FIRST);
    ml_sol.AddSolution("Sol", LAGRANGE, fe_order);

    //Initialize (update Init(...) function)
    ml_sol.Initialize("Sol");

    //Set Boundary (update Dirichlet(...) function)
    ml_sol.InitializeBdc();
    
    std::vector<std::string> facenamearray;
    std::vector<ParsedFunction> parsedfunctionarray;
    std::vector<BDCType> bdctypearray;
       
    unsigned int bdcsize = inputparser->getSize("multilevel_solution.multilevel_mesh.first.variable.first.boundary_conditions");
    for(unsigned int index=0; index<bdcsize; ++index) {
      std::string facename = inputparser->getValueFromArray("multilevel_solution.multilevel_mesh.first.variable.first.boundary_conditions", index, "facename", "top");
      facenamearray.push_back(facename);
 
      BDCType bdctype = inputparser->getValueFromArray("multilevel_solution.multilevel_mesh.first.variable.first.boundary_conditions", index, "bdc_type", DIRICHLET);
      bdctypearray.push_back(bdctype);

      std::string bdcfuncstr = inputparser->getValueFromArray("multilevel_solution.multilevel_mesh.first.variable.first.boundary_conditions", index, "bdc_func", "0.");
      ParsedFunction pfunc(bdcfuncstr, "x,y,z,t");
      parsedfunctionarray.push_back(pfunc);
    }
      
    for(int i=0; i<bdcsize; ++i) {
      ml_sol.SetBoundaryCondition("Sol",facenamearray[i],bdctypearray[i],false,&parsedfunctionarray[i]);
    }
    
    ml_sol.GenerateBdc();
    
    
    MultiLevelProblem ml_prob(&ml_msh,&ml_sol);
    
    // add fluid material
    Parameter parameter(Lref,Uref);
    
    //BEGIN Poisson MultiLevel Problem
    std::cout << std::endl;
    std::cout << " PDE problem to solve: Poisson " << std::endl;

    LinearImplicitSystem & system2 = ml_prob.add_system<LinearImplicitSystem>("Poisson");
    system2.AddSolutionToSytemPDE("Sol");
    
    // reading source function
    std::string function;
    function = inputparser->getValue("func_source", "0.");
    std::string variables = "x";
    variables += ",y";
    variables += ",z";
    variables += ",t";
     
#ifdef HAVE_FPARSER       
    fpsource.SetExpression(function);
    fpsource.SetIndependentVariables(variables);
    fpsource.Parse();
#endif

    // Set MG Options
    system2.AttachAssembleFunction(AssemblePoissonMatrixandRhs);
    
    unsigned int max_number_linear_iteration = inputparser->getValue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.max_number_linear_iteration",6);
    system2.SetMaxNumberOfLinearIterations(max_number_linear_iteration);
    
    double abs_conv_tol = inputparser->getValue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.abs_conv_tol",1.e-08);
    system2.SetAbsoluteConvergenceTolerance(abs_conv_tol);
    
    MgType mgtype = inputparser->getValue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.mgtype",V_CYCLE);
    system2.SetMgType(mgtype);
    
    unsigned int npresmoothing = inputparser->getValue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.npresmoothing",1); 
    system2.SetNumberPreSmoothingStep(npresmoothing);
    
    unsigned int npostmoothing = inputparser->getValue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.npostmoothing",1); 
    system2.SetNumberPostSmoothingStep(npostmoothing);

    if(inputparser->isTrue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.smoother.type","gmres")) {
      system2.SetMgSmoother(GMRES_SMOOTHER);
    }
    else if(inputparser->isTrue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.smoother.type","asm")) {
      system2.SetMgSmoother(ASM_SMOOTHER);
    }
    else if(inputparser->isTrue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.smoother.type","vanka")) {
      system2.SetMgSmoother(VANKA_SMOOTHER);
    }
    
    system2.init();
    //common smoother option
    system2.SetSolverFineGrids(RICHARDSON);
    system2.SetTolerances(1.e-12,1.e-20,1.e+50,4);
    system2.SetPreconditionerFineGrids(SOR_PRECOND);
    //for Vanka and ASM smoothers
    system2.ClearVariablesToBeSolved();
    system2.AddVariableToBeSolved("All");
    if(inputparser->isTrue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.smoother.type","asm") ||
       inputparser->isTrue("multilevel_problem.multilevel_mesh.first.system.poisson.linear_solver.type.multigrid.smoother.type","vanka")) {
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
    vector<double> normal(3.0);
    double src_term = 0.;
    vector<double> xyzt(4,0.);
    ParsedFunction* bdcfunc = NULL;

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
            // *** Gauss point loop ***
            for(unsigned ig=0; ig < ml_prob._ml_msh->_type_elem[kelt][order_ind]->GetGaussPointNumber(); ig++) {
                // *** get Jacobian and test function and test function derivatives ***
                (ml_prob._ml_msh->_type_elem[kelt][order_ind]->*(ml_prob._ml_msh->_type_elem[kelt][order_ind])->Jacobian_ptr)(coordinates,ig,weight,phi,gradphi);
                //current solution
                double SolT=0;
                vector < double > gradSolT(dim,0.);
                for(unsigned ivar=0; ivar<dim; ivar++) {
                    gradSolT[ivar]=0;
                }

                xyzt.assign(4,0.);
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
                    
                   src_term = fpsource(&xyzt[0]);
                    
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
    
		    unsigned int face = -(mymsh->el->GetFaceElementIndex(kel,jface)+1) - 1;
		    
		    if(ml_sol->GetBoundaryCondition("Sol",face) == NEUMANN && !ml_sol->Ishomogeneous("Sol",face)) {
		      
		        bdcfunc = (ParsedFunction* )(ml_sol->GetBdcFunction("Sol", face));
			unsigned nve = mymsh->el->GetElementFaceDofNumber(kel,jface,order_ind);
                        const unsigned FELT[6][2]= {{3,3},{4,4},{3,4},{5,5},{5,5},{6,6}};
                        unsigned felt = FELT[kelt][jface<mymsh->el->GetElementFaceNumber(kel,0)];

                        for(unsigned i=0; i<nve; i++) {
                            unsigned inode=mymsh->el->GetFaceVertexIndex(kel,jface,i)-1u;
                            unsigned inode_coord_metis=mymsh->GetMetisDof(inode,2);

			    for(unsigned ivar=0; ivar<dim; ivar++) {
                              coordinates[ivar][i]=(*mymsh->_coordinate->_Sol[ivar])(inode_coord_metis);
                            }
                        }

                        if(felt != 6) 
			{
                          for(unsigned igs=0; igs < ml_prob._ml_msh->_type_elem[felt][order_ind]->GetGaussPointNumber(); igs++) {
                            (ml_prob._ml_msh->_type_elem[felt][order_ind]->*ml_prob._ml_msh->_type_elem[felt][order_ind]->Jacobian_sur_ptr)(coordinates,igs,weight,phi,gradphi,normal);

			    xyzt.assign(4,0.);
                            for(unsigned i=0; i<nve; i++) {
		              for(unsigned ivar=0; ivar<dim; ivar++) {
		               xyzt[ivar] += coordinates[ivar][i]*phi[i]; 
		              }
                            }
			    
			    // *** phi_i loop ***
                            for(unsigned i=0; i<nve; i++) {
                              double surfterm_g = (*bdcfunc)(&xyzt[0]);
                              double bdintegral = phi[i]*surfterm_g*weight;
                              unsigned int ilocalnode = mymsh->el->GetLocalFaceVertexIndex(kel, jface, i);
                              F[ilocalnode] += bdintegral;
                            }
                          }
			}
			else // 1D : the side elems are points and does not still exist the point elem
			{
			  // in 1D it is only one point
 			  xyzt[0] = coordinates[0][0];
 		          xyzt[1] = 0.;
 			  xyzt[2] = 0.;
 			  xyzt[3] = 0.;

			  double bdintegral = (*bdcfunc)(&xyzt[0]);
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



#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "NumericVector.hpp"
#include "Assemble_jacobian.hpp"

#define FACE_FOR_CONTROL 2  //we do control on the right (=2) face
#define AXIS_DIRECTION_CONTROL_SIDE  1  //change this accordingly to the other variable above
#include "../elliptic_param.hpp"

#define GROUP_INTERNAL  12
#define GROUP_EXTERNAL  13


using namespace femus;


double Solution_set_initial_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {

    double value = 0.;

    if(!strcmp(name,"state")) {
        value = 0.;
    }
    else if(!strcmp(name,"control")) {
        value = 0.;
    }
    else if(!strcmp(name,"adjoint")) {
        value = 0.;
    }
    else if(!strcmp(name,"adjoint_ext")) {
        value = 0.;
    }
    else if(!strcmp(name,"mu")) {
        value = 0.;
    }
    else if(!strcmp(name,"TargReg")) {
        value = ElementTargetFlag(x);
    }
    else if(!strcmp(name,"ContReg")) {
        value = ControlDomainFlag_external_restriction(x);
    }
    else if(!strcmp(name,"act_flag")) {
        value = 0.;
    }


    return value;
}


bool Solution_set_boundary_conditions(const std::vector < double >& x, const char name[], double& value, const int faceName, const double time) {

    bool dirichlet = true; //dirichlet
    value = 0.;

//   if(!strcmp(name,"control")) {
//       value = 0.;
//   if (faceName == 3)
//     dirichlet = false;
//
//   }

    if(!strcmp(name,"mu")) {
//       value = 0.;
//   if (faceName == 3)
        dirichlet = false;
    }

    return dirichlet;
}


void ComputeIntegral(const MultiLevelProblem& ml_prob);

void AssembleLiftExternalProblem(MultiLevelProblem& ml_prob);


int main(int argc, char** args) {

    // init Petsc-MPI communicator
    FemusInit mpinit(argc, args, MPI_COMM_WORLD);

    // ======= Files ========================
    Files files;
    files.CheckIODirectories();
    files.RedirectCout();

    // ======= Quad Rule ========================
    std::string fe_quad_rule("seventh");
    /* "seventh" is the order of accuracy that is used in the gauss integration scheme
         probably in the furure it is not going to be an argument of this function   */

    // ======= Mesh  ==================
    MultiLevelMesh ml_mesh;

    double scalingFactor = 1.;

    // read coarse level mesh and generate finers level meshes
    ml_mesh.ReadCoarseMesh("./input/ext_box.neu", fe_quad_rule.c_str(), scalingFactor);

    //ml_mesh.GenerateCoarseBoxMesh(NSUB_X,NSUB_Y,0,0.,1.,0.,1.,0.,0.,QUAD9,"seventh");
    unsigned numberOfUniformLevels = 4;
    unsigned numberOfSelectiveLevels = 0;
    ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
    ml_mesh.PrintInfo();
    ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);

    // ======= Solution  ==================
    MultiLevelSolution ml_sol(&ml_mesh);

    // add variables to ml_sol
    ml_sol.AddSolution("state", LAGRANGE, FIRST);
    ml_sol.AddSolution("control", LAGRANGE, FIRST);
    ml_sol.AddSolution("adjoint", LAGRANGE, FIRST);
    ml_sol.AddSolution("adjoint_ext", LAGRANGE, FIRST);
    ml_sol.AddSolution("mu", LAGRANGE, FIRST);
    ml_sol.AddSolution("TargReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
    ml_sol.AddSolution("ContReg",  DISCONTINUOUS_POLYNOMIAL, ZERO); //this variable is not solution of any eqn, it's just a given field
    const unsigned int fake_time_dep_flag = 2;  //this is needed to be able to use _SolOld
    const std::string act_set_flag_name = "act_flag";
    ml_sol.AddSolution(act_set_flag_name.c_str(), LAGRANGE, FIRST,fake_time_dep_flag);               //this variable is not solution of any eqn, it's just a given field

    // ======= Problem ========================
    MultiLevelProblem ml_prob(&ml_sol);

    ml_prob.SetFilesHandler(&files);

    ml_sol.Initialize("All");    // initialize all varaibles to zero

    ml_sol.Initialize("state",       Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize("control",     Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize("adjoint",     Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize("adjoint_ext", Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize("mu",          Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize("TargReg",     Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize("ContReg",     Solution_set_initial_conditions, & ml_prob);
    ml_sol.Initialize(act_set_flag_name.c_str(), Solution_set_initial_conditions, & ml_prob);

    // attach the boundary condition function and generate boundary data
    ml_sol.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions);
    ml_sol.GenerateBdc("state");
    ml_sol.GenerateBdc("control");
    ml_sol.GenerateBdc("adjoint");
    ml_sol.GenerateBdc("adjoint_ext");
    ml_sol.GenerateBdc("mu");  //we need add this to make the matrix iterations work...

    // ======= System ========================
    NonLinearImplicitSystemWithPrimalDualActiveSetMethod& system = ml_prob.add_system < NonLinearImplicitSystemWithPrimalDualActiveSetMethod > ("LiftRestr");

    system.SetActiveSetFlagName(act_set_flag_name);

    system.AddSolutionToSystemPDE("state");
    system.AddSolutionToSystemPDE("control");
    system.AddSolutionToSystemPDE("adjoint");
    system.AddSolutionToSystemPDE("adjoint_ext");
    system.AddSolutionToSystemPDE("mu");

    // attach the assembling function to system
    system.SetAssembleFunction(AssembleLiftExternalProblem);

    ml_sol.SetWriter(VTK);
    ml_sol.GetWriter()->SetDebugOutput(true);

    system.SetDebugNonlinear(true);
    system.SetDebugFunction(ComputeIntegral);
    // system.SetMaxNumberOfNonLinearIterations(4);

    // initilaize and solve the system
    system.init();
    system.MGsolve();

    // ======= Print ========================
    std::vector < std::string > variablesToBePrinted;
    variablesToBePrinted.push_back("all");
    ml_sol.GetWriter()->Write(files.GetOutputPath()/*DEFAULT_OUTPUTDIR*/, "biquadratic", variablesToBePrinted);

    return 0;

}


void AssembleLiftExternalProblem(MultiLevelProblem& ml_prob) {
    //  ml_prob is the global object from/to where get/set all the data

    //  level is the level of the PDE system to be assembled
    //  levelMax is the Maximum level of the MultiLevelProblem
    //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

    //  extract pointers to the several objects that we are going to use

    NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
    const unsigned level = mlPdeSys->GetLevelToAssemble();
    const bool assembleMatrix = mlPdeSys->GetAssembleMatrix();

    Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);

    MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
    Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

    LinearEquationSolver* pdeSys = mlPdeSys->_LinSolver[level];
    SparseMatrix*             KK = pdeSys->_KK;
    NumericVector*           RES = pdeSys->_RES;

    const unsigned  dim = msh->GetDimension();
    const unsigned dim2 = (3 * (dim - 1) + !(dim - 1));
    const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));

    const unsigned    iproc = msh->processor_id();

//***************************************************

    const int coords_fe_type = BIQUADR_FE;  //biquadratic

//************** geometry (at dofs) *************************************
    vector < vector < double > > coords_at_dofs(dim);
    vector < vector < double > > coords_at_dofs_bdry(dim);
    for (unsigned idim = 0; idim < dim; idim++) {
        coords_at_dofs[idim].reserve(max_size);
        coords_at_dofs_bdry[idim].reserve(max_size);
    }

//************** geometry (at quadrature points) *************************************
    vector < double > coord_at_qp(dim);


//************* shape functions (at dofs and quadrature points) **************************************
    double weight_qp = 0.;      // gauss point weight
    double weight_qp_bdry = 0.; // gauss point weight on the boundary

    vector < vector < double > > phi_fe_qp(NFE_FAMS);
    vector < vector < double > > phi_x_fe_qp(NFE_FAMS);
    vector < vector < double > > phi_xx_fe_qp(NFE_FAMS);

    for(int fe=0; fe < NFE_FAMS; fe++) {
        phi_fe_qp[fe].reserve(max_size);
        phi_x_fe_qp[fe].reserve(max_size*dim);
        phi_xx_fe_qp[fe].reserve(max_size*(3*(dim-1)));
    }

//************* bdry shape functions (at dofs and quadrature points) **************************************
    vector < vector < double > > phi_fe_qp_bdry(NFE_FAMS);
    vector < vector < double > > phi_x_fe_qp_bdry(NFE_FAMS);

    for(int fe=0; fe < NFE_FAMS; fe++) {
        phi_fe_qp_bdry[fe].reserve(max_size);
        phi_x_fe_qp_bdry[fe].reserve(max_size * dim);
    }

//********************* vol-at-bdry adjoint *******************
    vector <double> phi_adj_vol_at_bdry;
    phi_adj_vol_at_bdry.reserve(max_size);
    vector <double> phi_adj_x_vol_at_bdry;
    phi_adj_x_vol_at_bdry.reserve(max_size * dim);
    vector <double> sol_adj_x_vol_at_bdry_gss(dim);
//***************************************************
    
//********************* vol-at-bdry adjoint_ext *******************
    vector <double> phi_adj_ext_vol_at_bdry;
    phi_adj_ext_vol_at_bdry.reserve(max_size);
    vector <double> phi_adj_ext_x_vol_at_bdry;
    phi_adj_ext_x_vol_at_bdry.reserve(max_size * dim);
    vector <double> sol_adj_ext_x_vol_at_bdry_gss(dim);
//***************************************************

    //************** act flag ****************************
    std::string act_flag_name = "act_flag";
    unsigned int solIndex_act_flag = ml_sol->GetIndex(act_flag_name.c_str());
    unsigned int solFEType_act_flag = ml_sol->GetSolutionType(solIndex_act_flag);
    if(sol->GetSolutionTimeOrder(solIndex_act_flag) == 2) {
        *(sol->_SolOld[solIndex_act_flag]) = *(sol->_Sol[solIndex_act_flag]);
    }

    //********* variables for ineq constraints *****************
    const int ineq_flag = INEQ_FLAG;
    const double c_compl = C_COMPL;
    vector < double/*int*/ >  sol_actflag;
    sol_actflag.reserve(max_size); //flag for active set
    vector < double >  ctrl_lower;
    ctrl_lower.reserve(max_size);
    vector < double >  ctrl_upper;
    ctrl_upper.reserve(max_size);
    //***************************************************

//***************************************************
//********* WHOLE SET OF VARIABLES ******************
    const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();

    enum Sol_pos {pos_state=0, pos_ctrl, pos_adj, pos_adj_ext, pos_mu}; //these are known at compile-time

    assert(pos_state   == mlPdeSys->GetSolPdeIndex("state"));
    assert(pos_ctrl    == mlPdeSys->GetSolPdeIndex("control"));
    assert(pos_adj     == mlPdeSys->GetSolPdeIndex("adjoint"));
    assert(pos_adj_ext == mlPdeSys->GetSolPdeIndex("adjoint_ext"));
    assert(pos_mu      == mlPdeSys->GetSolPdeIndex("mu"));


    vector < std::string > Solname(n_unknowns);
    Solname[0] = "state";
    Solname[1] = "control";
    Solname[2] = "adjoint";
    Solname[3] = "adjoint_ext";
    Solname[4] = "mu";

    vector < unsigned > SolPdeIndex(n_unknowns);
    vector < unsigned > SolIndex(n_unknowns);
    vector < unsigned > SolFEType(n_unknowns);


    for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
        SolPdeIndex[ivar] = mlPdeSys->GetSolPdeIndex(Solname[ivar].c_str());
        SolIndex[ivar]    = ml_sol->GetIndex        (Solname[ivar].c_str());
        SolFEType[ivar]   = ml_sol->GetSolutionType(SolIndex[ivar]);
    }

    vector < unsigned > Sol_n_el_dofs(n_unknowns);

//***************************************************
    //----------- quantities (at dof objects) ------------------------------
    vector< int >       L2G_dofmap_AllVars;
    L2G_dofmap_AllVars.reserve( n_unknowns * max_size );
    vector < vector < int > >     L2G_dofmap(n_unknowns);
    for(int i = 0; i < n_unknowns; i++) {
        L2G_dofmap[i].reserve(max_size);
    }
    vector < vector < double > >  sol_eldofs(n_unknowns);
    for(int k = 0; k < n_unknowns; k++) {
        sol_eldofs[k].reserve(max_size);
    }
    vector< double > Res;
    Res.reserve( n_unknowns * max_size);
    vector < double > Jac;
    Jac.reserve( n_unknowns * max_size * n_unknowns * max_size);

//***************************************************
    //------------ quantities (at quadrature points) ---------------------
    vector<double>        sol_qp(n_unknowns);
    vector< vector<double> > sol_grad_qp(n_unknowns);

    std::fill(sol_qp.begin(), sol_qp.end(), 0.);
    for (unsigned  k = 0; k < n_unknowns; k++) {
        sol_grad_qp[k].resize(dim);
        std::fill(sol_grad_qp[k].begin(), sol_grad_qp[k].end(), 0.);
    }


//********************* DATA ************************
    const double u_des = DesiredTarget();
    const double alpha = ALPHA_CTRL_VOL;
    const double beta  = BETA_CTRL_VOL;
    const double penalty_strong_ctrl = 1.e30;
    const double penalty_strong_u = 1.e30;
    const double penalty_interface = 1.e10;         //penalty for u=q
//***************************************************

    RES->zero();
    if (assembleMatrix)  KK->zero();


    // element loop: each process loops only on the elements that owns
    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        int group_flag         = msh->GetElementGroup(iel);
        short unsigned ielGeom = msh->GetElementType(iel);    // element geometry type
//    std::cout << " ======= grp_flag === " << group_flag << " ================== " << std::endl;
//     int face_no         = msh->GetElementFaceNumber(iel);
//     std::cout << " ======= face# === " << face_no << " ================== " << std::endl;

//******************** GEOMETRY *********************
        unsigned nDofx = msh->GetElementDofNumber(iel, coords_fe_type);    // number of coordinate element dofs
        for (int i = 0; i < dim; i++)  coords_at_dofs[i].resize(nDofx);
        // local storage of coordinates
        for (unsigned i = 0; i < nDofx; i++) {
            unsigned xDof  = msh->GetSolutionDof(i, iel, coords_fe_type);  // global to global mapping between coordinates node and coordinate dof

            for (unsigned jdim = 0; jdim < dim; jdim++) {
                coords_at_dofs[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);      // global extraction and local storage for the element coordinates
            }
        }

        // elem average point
        vector < double > elem_center(dim);
        for (unsigned j = 0; j < dim; j++) {
            elem_center[j] = 0.;
        }
        for (unsigned j = 0; j < dim; j++) {
            for (unsigned i = 0; i < nDofx; i++) {
                elem_center[j] += coords_at_dofs[j][i];
            }
        }

        for (unsigned j = 0; j < dim; j++) {
            elem_center[j] = elem_center[j]/nDofx;
        }
//***************************************************

//****** set target domain flag *********************
        int target_flag = 0;
        target_flag = ElementTargetFlag(elem_center);
//***************************************************

        //all vars###################################################################
        for (unsigned  k = 0; k < n_unknowns; k++) {
            unsigned  ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
            Sol_n_el_dofs[k] = ndofs_unk;
            sol_eldofs[k].resize(ndofs_unk);
            L2G_dofmap[k].resize(ndofs_unk);
            for (unsigned i = 0; i < ndofs_unk; i++) {
                unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);                        // global to global mapping between solution node and solution dof // via local to global solution node
                sol_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);                            // global extraction and local storage for the solution
                L2G_dofmap[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
            }
        }
        //all vars###################################################################


//************** update active set flag for current nonlinear iteration ****************************
// 0: inactive; 1: active_a; 2: active_b
        assert(Sol_n_el_dofs[pos_mu] == Sol_n_el_dofs[pos_ctrl]);
        sol_actflag.resize(Sol_n_el_dofs[pos_mu]);
        ctrl_lower.resize(Sol_n_el_dofs[pos_mu]);
        ctrl_upper.resize(Sol_n_el_dofs[pos_mu]);
        std::fill(sol_actflag.begin(), sol_actflag.end(), 0);
        std::fill(ctrl_lower.begin(), ctrl_lower.end(), 0.);
        std::fill(ctrl_upper.begin(), ctrl_upper.end(), 0.);

        for (unsigned i = 0; i < sol_actflag.size(); i++) {
            std::vector<double> node_coords_i(dim,0.);
            for (unsigned d = 0; d < dim; d++) node_coords_i[d] = coords_at_dofs[d][i];
            ctrl_lower[i] = InequalityConstraint(node_coords_i,false);
            ctrl_upper[i] = InequalityConstraint(node_coords_i,true);

            if      ( (sol_eldofs[pos_mu][i] + c_compl * (sol_eldofs[pos_ctrl][i] - ctrl_lower[i] )) < 0 )  sol_actflag[i] = 1;
            else if ( (sol_eldofs[pos_mu][i] + c_compl * (sol_eldofs[pos_ctrl][i] - ctrl_upper[i] )) > 0 )  sol_actflag[i] = 2;
        }

//************** act flag ****************************
        unsigned nDof_act_flag  = msh->GetElementDofNumber(iel, solFEType_act_flag);    // number of solution element dofs

        for (unsigned i = 0; i < nDof_act_flag; i++) {
            unsigned solDof_mu = msh->GetSolutionDof(i, iel, solFEType_act_flag);
            (sol->_Sol[solIndex_act_flag])->set(solDof_mu,sol_actflag[i]);
        }



//******************** ALL VARS *********************
        //******************** ALL VARS *********************
        unsigned sum_Sol_n_el_dofs = 0;
        for (unsigned  k = 0; k < n_unknowns; k++) {
            sum_Sol_n_el_dofs += Sol_n_el_dofs[k];
        }
// TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable
        int nDof_max    =  0;
        for (unsigned  k = 0; k < n_unknowns; k++)     {
            if(Sol_n_el_dofs[k] > nDof_max)    nDof_max = Sol_n_el_dofs[k];
        }

        Res.resize(sum_Sol_n_el_dofs);
        std::fill(Res.begin(), Res.end(), 0.);
        Jac.resize(sum_Sol_n_el_dofs * sum_Sol_n_el_dofs);
        std::fill(Jac.begin(), Jac.end(), 0.);

        L2G_dofmap_AllVars.resize(0);
        for (unsigned  k = 0; k < n_unknowns; k++)     L2G_dofmap_AllVars.insert(L2G_dofmap_AllVars.end(),L2G_dofmap[k].begin(),L2G_dofmap[k].end());
//***************************************************


        std::vector <double> interface_node_flag(Sol_n_el_dofs[pos_state]); //flag for boundary interface
        std::fill(interface_node_flag.begin(), interface_node_flag.end(), 0.);

        //************ Boundary loops ***************************************
        vector<double> normal_qp(dim,0);

        for(unsigned jface=0; jface < msh->GetElementFaceNumber(iel); jface++) {

//========= compute coordinates of boundary nodes on each element ==========================================
            unsigned nDofx_bdry    = msh->GetElementFaceDofNumber(iel,jface,coords_fe_type);
            unsigned nDofu_bdry    = msh->GetElementFaceDofNumber(iel,jface,SolFEType[pos_state]);
            unsigned nDofctrl_bdry = msh->GetElementFaceDofNumber(iel,jface,SolFEType[pos_ctrl]);
            if (nDofu_bdry != nDofctrl_bdry) {
                std::cout << "State and control need to have the same FE space" << std::endl;
                abort();
            }
            for (unsigned idim = 0; idim < dim; idim++) {
                coords_at_dofs_bdry[idim].resize(nDofx_bdry);
            }
            const unsigned felt_bdry = msh->GetElementFaceType(iel, jface);
            for(unsigned i=0; i < nDofx_bdry; i++) {
                unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);
                unsigned iDof = msh->GetSolutionDof(i_vol, iel, coords_fe_type);
                for(unsigned idim=0; idim<dim; idim++) {
                    coords_at_dofs_bdry[idim][i]=(*msh->_topology->_Sol[idim])(iDof);
                }
            }
//===================================================

//============== face elem average point =====================================
            vector < double > elem_center_bdry(dim);
            for (unsigned j = 0; j < dim; j++) {
                elem_center_bdry[j] = 0.;
            }
            for (unsigned j = 0; j < dim; j++) {
                for (unsigned i = 0; i < nDofx_bdry; i++) {
                    elem_center_bdry[j] += coords_at_dofs_bdry[j][i];
                }
            }
            for (unsigned j = 0; j < dim; j++) {
                elem_center_bdry[j] = elem_center_bdry[j]/nDofx_bdry;
            }
//===================================================


//============ find interface boundary elements (now we do with coordinates, later we can do also with flag in gambit) =======================================
            double my_eps = 1.e-6;
            if (elem_center_bdry[0] > 1. - my_eps   && elem_center_bdry[0] < 1.   + my_eps  &&
                    elem_center_bdry[1] > 0.25 - my_eps && elem_center_bdry[1] < 0.75 + my_eps
               ) {
                std::cout << " bdry elem on interface with center " << "(" << elem_center_bdry[0] << "," << elem_center_bdry[1] << ")" << std::endl;
                for (int i_bdry = 0; i_bdry < nDofu_bdry; i_bdry++)  {
                    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
                    interface_node_flag[i_vol] = 1.;
                }


//========= initialize gauss quantities on the boundary ============================================

                for(unsigned ig_bdry=0; ig_bdry < msh->_finiteElement[felt_bdry][SolFEType[pos_ctrl]]->GetGaussPointNumber(); ig_bdry++) {

                    // *** get gauss point weight, test function and test function partial derivatives ***
                    for(int fe=0; fe < NFE_FAMS; fe++) {
                        msh->_finiteElement[felt_bdry][fe]->JacobianSur(coords_at_dofs_bdry, ig_bdry, weight_qp_bdry, phi_fe_qp_bdry[fe], phi_x_fe_qp_bdry[fe], normal_qp);
                    }
                    //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
                    msh->_finiteElement[felt_bdry][coords_fe_type]->JacobianSur(coords_at_dofs_bdry, ig_bdry, weight_qp_bdry, phi_fe_qp_bdry[coords_fe_type], phi_x_fe_qp_bdry[coords_fe_type], normal_qp);

                    if (ielGeom != QUAD) {
                        std::cout << "VolumeShapeAtBoundary not implemented" << std::endl;
                        abort();
                    }
                    msh->_finiteElement[ielGeom][SolFEType[pos_adj]]->VolumeShapeAtBoundary(coords_at_dofs, coords_at_dofs_bdry, jface, ig_bdry, phi_adj_vol_at_bdry, phi_adj_x_vol_at_bdry);
                    msh->_finiteElement[ielGeom][SolFEType[pos_adj]]->VolumeShapeAtBoundary(coords_at_dofs, coords_at_dofs_bdry, jface, ig_bdry, phi_adj_ext_vol_at_bdry, phi_adj_ext_x_vol_at_bdry);

                    
//=============== grad dot n for residual =========================================
//     compute gauss quantities on the boundary through VOLUME interpolation
                    std::fill(sol_adj_x_vol_at_bdry_gss.begin(), sol_adj_x_vol_at_bdry_gss.end(), 0.);
                    for (int iv = 0; iv < Sol_n_el_dofs[pos_adj]; iv++)  {

                        for (int d = 0; d < dim; d++) {
//    std::cout << " ivol " << iv << std::endl;
//    std::cout << " adj dofs " << sol_adj[iv] << std::endl;
                            sol_adj_x_vol_at_bdry_gss[d] += sol_eldofs[pos_adj][iv] * phi_adj_x_vol_at_bdry[iv * dim + d];//notice that the convention of the orders x y z is different from vol to bdry
                        }
                    }

                    double grad_adj_dot_n_res = 0.;
                    for (unsigned d = 0; d < dim; d++) {
                        grad_adj_dot_n_res += sol_adj_x_vol_at_bdry_gss[d] * normal_qp[d];
                    }
//=============== grad dot n  for residual =========================================

//=============== grad dot n for residual =========================================
//     compute gauss quantities on the boundary through VOLUME interpolation
                    std::fill(sol_adj_ext_x_vol_at_bdry_gss.begin(), sol_adj_ext_x_vol_at_bdry_gss.end(), 0.);
                    for (int iv = 0; iv < Sol_n_el_dofs[pos_adj_ext]; iv++)  {

                        for (int d = 0; d < dim; d++) {
                            sol_adj_ext_x_vol_at_bdry_gss[d] += sol_eldofs[pos_adj_ext][iv] * phi_adj_ext_x_vol_at_bdry[iv * dim + d];//notice that the convention of the orders x y z is different from vol to bdry
                        }
                    }

                    double grad_adj_ext_dot_n_res = 0.;
                    for (unsigned d = 0; d < dim; d++) {
                        grad_adj_ext_dot_n_res += sol_adj_ext_x_vol_at_bdry_gss[d] * normal_qp[d];
                    }
//=============== grad dot n  for residual =========================================


                    for (int i_bdry = 0; i_bdry < nDofu_bdry; i_bdry++)  {
                        unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);

//============ Bdry Residuals ==================

                        if ( group_flag == GROUP_INTERNAL ) Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_state,i_vol)  ]  +=  -  weight_qp_bdry *  ( grad_adj_dot_n_res * phi_fe_qp_bdry[SolFEType[pos_state]][i_bdry] );

                        
                        if ( group_flag == GROUP_EXTERNAL ) Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_ctrl,i_vol)  ]  +=  -  weight_qp_bdry *  ( grad_adj_ext_dot_n_res * phi_fe_qp_bdry[SolFEType[pos_ctrl]][i_bdry] );

                        if ( group_flag == GROUP_INTERNAL ) Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_adj,i_vol)   ]  +=  -  penalty_interface * ( sol_eldofs[pos_state][i_vol] - sol_eldofs[pos_ctrl][i_vol] ) ;    // u = q
                        
                        if ( group_flag == GROUP_EXTERNAL ) Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_adj_ext,i_vol)   ]  +=  - penalty_interface * weight_qp_bdry * phi_fe_qp_bdry[SolFEType[pos_adj_ext]][i_bdry] * ( grad_adj_dot_n_res - grad_adj_ext_dot_n_res ) ;
//============ Bdry Residuals ==================

//============ Bdry Jacobians ==================
                        for(unsigned j_bdry=0; j_bdry < nDofu_bdry; j_bdry ++) {
                            unsigned int j_vol = msh->GetLocalFaceVertexIndex(iel, jface, j_bdry);

//============ u = q =============================
                            if (i_vol == j_vol)  {
                                if ( group_flag == GROUP_INTERNAL ) Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_adj, pos_state, i_vol, j_vol) ]  += penalty_interface *  ( 1.);
                                if ( group_flag == GROUP_INTERNAL ) Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_adj, pos_ctrl, i_vol, j_vol) ]   += penalty_interface *  (-1.);
                            }
//============ u = q =============================

                        } //end j_bdry

//===================loop over j in the VOLUME (while i is in the boundary)
                        for(unsigned j=0; j < nDof_max; j ++) {

//=============== grad dot n  =========================================
                            double grad_adj_dot_n_mat = 0.;
                            for(unsigned d=0; d<dim; d++) {
                                grad_adj_dot_n_mat += phi_adj_x_vol_at_bdry[j * dim + d] * normal_qp[d];  //notice that the convention of the orders x y z is different from vol to bdry
                            }
//=============== grad dot n  =========================================

//=============== grad dot n  =========================================
                            double grad_adj_ext_dot_n_mat = 0.;
                            for(unsigned d=0; d<dim; d++) {
                                grad_adj_ext_dot_n_mat += phi_adj_ext_x_vol_at_bdry[j * dim + d] * normal_qp[d];  //notice that the convention of the orders x y z is different from vol to bdry
                            }
//=============== grad dot n  =========================================

                            if ( group_flag == GROUP_INTERNAL ) {
                                Jac[  assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_state, pos_adj, i_vol, j) ]  += weight_qp_bdry * grad_adj_dot_n_mat * phi_fe_qp_bdry[SolFEType[pos_state]][i_bdry];
                            }
                            
                            if ( group_flag == GROUP_EXTERNAL ) {
                                Jac[  assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_ctrl, pos_adj_ext, i_vol, j) ]  += weight_qp_bdry * grad_adj_ext_dot_n_mat * phi_fe_qp_bdry[SolFEType[pos_ctrl]][i_bdry];
                            }
                            
                            if (i_vol == j)  {
                                if ( group_flag == GROUP_EXTERNAL ) Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_adj_ext, pos_adj, i_vol, j) ]  += penalty_interface * weight_qp_bdry * phi_fe_qp_bdry[SolFEType[pos_adj_ext]][i_bdry] * ( 1.) * grad_adj_dot_n_mat;
                                if ( group_flag == GROUP_EXTERNAL ) Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_adj_ext, pos_adj_ext, i_vol, j) ]   += penalty_interface * weight_qp_bdry * phi_fe_qp_bdry[SolFEType[pos_adj_ext]][i_bdry] * ( -1.) * grad_adj_ext_dot_n_mat;
                            }

                        } //end j

                    } //end i_bdry

                }  //end ig_bdry loop


            }  //bdry interface elements

        }  //end boundary face loop

        //************ Boundary loops ***************************************


// *** Gauss point loop ***
        for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][coords_fe_type]->GetGaussPointNumber(); ig++) {

            // *** get gauss point weight, test function and test function partial derivatives ***
            for(int fe=0; fe < NFE_FAMS; fe++) {
                msh->_finiteElement[ielGeom][fe]->Jacobian(coords_at_dofs,ig,weight_qp,phi_fe_qp[fe],phi_x_fe_qp[fe],phi_xx_fe_qp[fe]);
            }
            //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
            msh->_finiteElement[ielGeom][coords_fe_type]->Jacobian(coords_at_dofs,ig,weight_qp,phi_fe_qp[coords_fe_type],phi_x_fe_qp[coords_fe_type],phi_xx_fe_qp[coords_fe_type]);

//========= fill gauss value xyz ==================
            std::fill(coord_at_qp.begin(), coord_at_qp.end(), 0.);
            for (unsigned  d = 0; d < dim; d++) {
                for (unsigned i = 0; i < coords_at_dofs[d].size(); i++) {
                    coord_at_qp[d] += coords_at_dofs[d][i] * phi_fe_qp[coords_fe_type][i];
                }
            }
            //========= fill gauss value xyz ==================

//========= fill gauss value quantities ==================
            std::fill(sol_qp.begin(), sol_qp.end(), 0.);
            for (unsigned  k = 0; k < n_unknowns; k++) {
                std::fill(sol_grad_qp[k].begin(), sol_grad_qp[k].end(), 0.);
            }

            for (unsigned  k = 0; k < n_unknowns; k++) {
                for (unsigned i = 0; i < Sol_n_el_dofs[k]; i++) {
                    sol_qp[k]    += sol_eldofs[k][i] *   phi_fe_qp[SolFEType[k]][i];
                    for (unsigned d = 0; d < dim; d++)   sol_grad_qp[k][d] += sol_eldofs[k][i] * phi_x_fe_qp[SolFEType[k]][i * dim + d];
                }
            }
//========= fill gauss value quantities ==================





//==========FILLING WITH THE EQUATIONS ===========
            // *** phi_i loop ***
            for (unsigned i = 0; i < nDof_max; i++) {

                double laplace_rhs_du_adj_i = 0.;
                double laplace_rhs_dadj_u_i = 0.;
                double laplace_rhs_dctrl_adj_ext_i = 0.;
                double laplace_rhs_dadj_ext_ctrl_i = 0.;
                double laplace_rhs_dctrl_ctrl_i = 0.;

                for (unsigned kdim = 0; kdim < dim; kdim++) {
                    if ( i < Sol_n_el_dofs[pos_state] )       laplace_rhs_du_adj_i          +=  (phi_x_fe_qp[SolFEType[pos_state]]   [i * dim + kdim] * sol_grad_qp[pos_adj][kdim]);
                    if ( i < Sol_n_el_dofs[pos_adj] )         laplace_rhs_dadj_u_i          +=  (phi_x_fe_qp[SolFEType[pos_adj]]     [i * dim + kdim] * sol_grad_qp[pos_state][kdim]);
                    if ( i < Sol_n_el_dofs[pos_ctrl] )        laplace_rhs_dctrl_adj_ext_i   +=  (phi_x_fe_qp[SolFEType[pos_ctrl]]    [i * dim + kdim] * sol_grad_qp[pos_adj_ext][kdim]);
                    if ( i < Sol_n_el_dofs[pos_adj_ext] )     laplace_rhs_dadj_ext_ctrl_i   +=  (phi_x_fe_qp[SolFEType[pos_adj_ext]] [i * dim + kdim] * sol_grad_qp[pos_ctrl][kdim]);
                    if ( i < Sol_n_el_dofs[pos_ctrl] )        laplace_rhs_dctrl_ctrl_i      +=  (phi_x_fe_qp[SolFEType[pos_ctrl]]    [i * dim + kdim] * sol_grad_qp[pos_ctrl][kdim]);
                }

//======================Volume Residuals=======================
                if ( group_flag == GROUP_INTERNAL )     {
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_state,i) ]  += - weight_qp * (target_flag * phi_fe_qp[SolFEType[pos_state]][i] * ( sol_qp[pos_state] - u_des) - laplace_rhs_du_adj_i);
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_ctrl,i)]    += - (1 - interface_node_flag[i]) * penalty_strong_ctrl * (sol_eldofs[pos_ctrl][i] - 0.);
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_adj,i)]     += - weight_qp *  ( - laplace_rhs_dadj_u_i    - 0.) ;
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_adj_ext,i)] += - (1 - interface_node_flag[i]) * penalty_strong_ctrl * ( sol_eldofs[pos_adj_ext][i] - 0.);
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_mu,i)]      += - (1 - interface_node_flag[i]) * penalty_strong_ctrl * ( sol_eldofs[pos_mu][i] - 0.);
                }

                else if ( group_flag == GROUP_EXTERNAL )  {
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_state,i) ]  += - (1 - interface_node_flag[i]) * penalty_strong_u * (sol_eldofs[pos_state][i] - 0.);
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_ctrl,i) ]   += - weight_qp * ( alpha * phi_fe_qp[SolFEType[pos_ctrl]][i] * sol_qp[pos_ctrl]
                            + beta * laplace_rhs_dctrl_ctrl_i
                            - laplace_rhs_dctrl_adj_ext_i /*+ 7000. * phi_fe_qp[SolFEType[pos_ctrl]][i]*/);
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_adj,i)]     += - (1 - interface_node_flag[i]) * penalty_strong_u * (  sol_eldofs[pos_adj][i] - 0.);
                    Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs,pos_adj_ext,i)] += - weight_qp *  ( - laplace_rhs_dadj_ext_ctrl_i - 0.) ;
                }
//======================Volume Residuals=======================

                if (assembleMatrix) {

                    // *** phi_j loop ***
                    for (unsigned j = 0; j < nDof_max; j++) {

                        double laplace_mat_du_adj = 0.;
                        double laplace_mat_dadj_u = 0.;
                        double laplace_mat_dctrl_adj_ext = 0.;
                        double laplace_mat_dadj_ext_ctrl = 0.;
                        double laplace_mat_dctrl_ctrl = 0.;

                        for (unsigned kdim = 0; kdim < dim; kdim++) {
                            if ( i < Sol_n_el_dofs[pos_state]   && j < Sol_n_el_dofs[pos_adj] )      laplace_mat_du_adj        += (phi_x_fe_qp[SolFEType[pos_state]]  [i * dim + kdim] * phi_x_fe_qp[SolFEType[pos_adj]]    [j * dim + kdim]);
                            if ( i < Sol_n_el_dofs[pos_adj]     && j < Sol_n_el_dofs[pos_state] )    laplace_mat_dadj_u        += (phi_x_fe_qp[SolFEType[pos_adj]]    [i * dim + kdim] * phi_x_fe_qp[SolFEType[pos_state]]  [j * dim + kdim]);  //equal to the previous
                            if ( i < Sol_n_el_dofs[pos_ctrl]    && j < Sol_n_el_dofs[pos_adj_ext] )  laplace_mat_dctrl_adj_ext += (phi_x_fe_qp[SolFEType[pos_ctrl]]   [i * dim + kdim] * phi_x_fe_qp[SolFEType[pos_adj_ext]][j * dim + kdim]);
                            if ( i < Sol_n_el_dofs[pos_adj_ext] && j < Sol_n_el_dofs[pos_ctrl] )     laplace_mat_dadj_ext_ctrl += (phi_x_fe_qp[SolFEType[pos_adj_ext]][i * dim + kdim] * phi_x_fe_qp[SolFEType[pos_ctrl]]   [j * dim + kdim]);  //equal to the previous
                            if ( i < Sol_n_el_dofs[pos_ctrl]    && j < Sol_n_el_dofs[pos_ctrl] )     laplace_mat_dctrl_ctrl    += (phi_x_fe_qp[SolFEType[pos_ctrl]]   [i * dim + kdim] * phi_x_fe_qp[SolFEType[pos_ctrl]]   [j * dim + kdim]);
                        }


                        if ( group_flag == GROUP_INTERNAL ) {

                            Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_state, pos_state, i, j) ] += weight_qp * target_flag * phi_fe_qp[SolFEType[pos_state]][j] * phi_fe_qp[SolFEType[pos_state]][i];
                            Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_state, pos_adj, i, j) ]   += weight_qp * (-1.) * laplace_mat_du_adj;

                            if ( i==j ) Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_ctrl, pos_ctrl, i, j) ]      += (1 - interface_node_flag[i]) * penalty_strong_ctrl;
                            
                                        Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_adj, pos_state, i, j) ]   += weight_qp * (-1.) * laplace_mat_dadj_u;

                            if ( i==j ) Jac[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_adj_ext, pos_adj_ext, i, j)]  += (1 - interface_node_flag[i]) * penalty_strong_ctrl;
                            if ( i==j ) Jac[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_mu, pos_mu, i, j) ]           += (1 - interface_node_flag[i]) * penalty_strong_ctrl;
                            

                        }

                        else if ( group_flag == GROUP_EXTERNAL ) {

                            if (  i==j ) Jac[  assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_state, pos_state, i, j) ]  += (1 - interface_node_flag[i]) * 1. * penalty_strong_u;
                            
                            Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_ctrl, pos_ctrl, i, j) ]     += weight_qp * ( 
                                 alpha * phi_fe_qp[SolFEType[pos_ctrl]][i] * phi_fe_qp[SolFEType[pos_ctrl]][j]  
                                + beta * laplace_mat_dctrl_ctrl
                                    );

                            Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_ctrl, pos_adj_ext, i, j) ]  += weight_qp * (-1.) * laplace_mat_dctrl_adj_ext;

                            if ( i==j ) Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_adj, pos_adj, i, j) ]       += (1 - interface_node_flag[i]) * 1. * penalty_strong_u;

                            Jac[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, pos_adj_ext, pos_ctrl, i, j) ]  += weight_qp * (-1.) * laplace_mat_dadj_ext_ctrl;

                        }

                    } // end phi_j loop

                } // endif assemble_matrix

            } // end phi_i loop

        } // end gauss point loop


//========== sum-based part ================================

        RES->add_vector_blocked(Res, L2G_dofmap_AllVars);
        if (assembleMatrix) KK->add_matrix_blocked(Jac, L2G_dofmap_AllVars, L2G_dofmap_AllVars);


//========== dof-based part, without summation

//============= delta_mu row ===============================
        std::vector<double> Res_mu (Sol_n_el_dofs[pos_mu]);
        std::fill(Res_mu.begin(),Res_mu.end(), 0.);
        for (unsigned i = 0; i < sol_actflag.size(); i++) {
            if (sol_actflag[i] == 0) {  //inactive
                Res_mu [i] = - ineq_flag * ( 1. * sol_eldofs[pos_mu][i] - 0. );
            }
            else if (sol_actflag[i] == 1) {  //active_a
                Res_mu [i] = - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl][i] - c_compl * ctrl_lower[i]);
            }
            else if (sol_actflag[i] == 2) {  //active_b
                Res_mu [i]  = - ineq_flag * ( c_compl * sol_eldofs[pos_ctrl][i] - c_compl * ctrl_upper[i]);
            }
        }


        RES->insert(Res_mu, L2G_dofmap[pos_mu]);

        //============= delta_ctrl-delta_mu row ===============================
        KK->matrix_set_off_diagonal_values_blocked( L2G_dofmap[pos_ctrl], L2G_dofmap[pos_mu], ineq_flag * 1.);

        //============= delta_mu-delta_ctrl row ===============================
        for (unsigned i = 0; i < sol_actflag.size(); i++) if (sol_actflag[i] != 0 ) sol_actflag[i] = ineq_flag * c_compl;

        KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_ctrl], sol_actflag);

        //============= delta_mu-delta_mu row ===============================
        for (unsigned i = 0; i < sol_actflag.size(); i++) sol_actflag[i] =  ineq_flag * (1 - sol_actflag[i]/c_compl)  + (1-ineq_flag) * 1.;

        KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_mu], sol_actflag );

    assemble_jacobian<double,double>::print_element_residual(iel, Res, Sol_n_el_dofs, 9, 5);
    assemble_jacobian<double,double>::print_element_jacobian(iel, Jac, Sol_n_el_dofs, 9, 5);
    
    } //end element loop for each process


    // ***************** END ASSEMBLY *******************
    //this part is to put a bunch of 1's in the residual, so it is not based on the element loop
    unsigned int ctrl_index = mlPdeSys->GetSolPdeIndex("control");

    unsigned int global_ctrl_size = pdeSys->KKoffset[ctrl_index+1][iproc] - pdeSys->KKoffset[ctrl_index][iproc];

    std::vector<double>  one_times_mu(global_ctrl_size, 0.);
    std::vector<int>    positions(global_ctrl_size);

    for (unsigned i = 0; i < positions.size(); i++) {
        positions[i] = pdeSys->KKoffset[ctrl_index][iproc] + i;
        one_times_mu[i] = ineq_flag * 1. * (*sol->_Sol[SolIndex[pos_mu]])(i/*position_mu_i*/) ;
    }
    RES->add_vector_blocked(one_times_mu, positions);


    RES->close();
    KK->close();

     //print JAC and RES to files
    const unsigned nonlin_iter = mlPdeSys->GetNonlinearIt();
    assemble_jacobian< double, double >::print_global_jacobian(assembleMatrix, ml_prob, KK, nonlin_iter);
    assemble_jacobian< double, double >::print_global_residual(ml_prob, RES, nonlin_iter);

    return;
}



void ComputeIntegral(const MultiLevelProblem& ml_prob)    {


    const NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys  = &ml_prob.get_system<NonLinearImplicitSystemWithPrimalDualActiveSetMethod> ("LiftRestr");   // pointer to the linear implicit system named "LiftRestr"
    const unsigned level         = mlPdeSys->GetLevelToAssemble();

    Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);            // pointer to the mesh (level) object

    MultiLevelSolution*    ml_sol = ml_prob._ml_sol;                             // pointer to the multilevel solution object
    Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

    const unsigned           dim = msh->GetDimension();                         // get the domain dimension of the problem
    const unsigned          dim2 = (3 * (dim - 1) + !(dim - 1));                // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
    const unsigned       max_size = static_cast< unsigned >(ceil(pow(3, dim)));  // conservative: based on line3, quad9, hex27

    const unsigned         iproc = msh->processor_id();                         // get the process_id (for parallel computation)

//***************************************************
    vector < vector < double > > coords_at_dofs(dim);   // local coordinates
    unsigned coords_fe_type  = BIQUADR_FE;  // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
    for (unsigned i = 0; i < dim; i++) {
        coords_at_dofs[i].reserve(max_size);
    }
//***************************************************

//***************************************************
    double weight_qp; // gauss point weight

//***************************************************
    double alpha = ALPHA_CTRL_VOL;
    double beta  = BETA_CTRL_VOL;

//******************** state ************************
//***************************************************
    vector <double> phi_u;    // local test function
    vector <double> phi_u_x;  // local test function first order partial derivatives
    vector <double> phi_u_xx; // local test function second order partial derivatives

    phi_u.reserve(max_size);
    phi_u_x.reserve(max_size * dim);
    phi_u_xx.reserve(max_size * dim2);


    unsigned solIndex_u;
    solIndex_u = ml_sol->GetIndex("state");                    // get the position of "state" in the ml_sol object
    unsigned solType_u = ml_sol->GetSolutionType(solIndex_u);  // get the finite element type for "state"

    vector < double >  sol_u; // local solution
    sol_u.reserve(max_size);

    double u_gss = 0.;
//***************************************************
//***************************************************

//******************** control **********************
//***************************************************
    vector <double> phi_ctrl;    // local test function
    vector <double> phi_ctrl_x;  // local test function first order partial derivatives
    vector <double> phi_ctrl_xx; // local test function second order partial derivatives

    phi_ctrl.reserve(max_size);
    phi_ctrl_x.reserve(max_size * dim);
    phi_ctrl_xx.reserve(max_size * dim2);

    unsigned solIndex_ctrl;
    solIndex_ctrl = ml_sol->GetIndex("control");
    unsigned solType_ctrl = ml_sol->GetSolutionType(solIndex_ctrl);

    vector < double >  sol_ctrl; // local solution
    sol_ctrl.reserve(max_size);

    double ctrl_gss = 0.;
    double ctrl_x_gss = 0.;
//***************************************************
//***************************************************


//********************* desired *********************
//***************************************************
    vector <double> phi_udes;    // local test function
    vector <double> phi_udes_x;  // local test function first order partial derivatives
    vector <double> phi_udes_xx; // local test function second order partial derivatives

    phi_udes.reserve(max_size);
    phi_udes_x.reserve(max_size * dim);
    phi_udes_xx.reserve(max_size * dim2);

    vector < double >  sol_udes; // local solution
    sol_udes.reserve(max_size);

    double udes_gss = 0.;
//***************************************************
//***************************************************

//********************* DATA ************************
    double u_des = DesiredTarget();
//***************************************************

    double integral_target = 0.;
    double integral_alpha  = 0.;
    double integral_beta   = 0.;


    // element loop: each process loops only on the elements that owns
    for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

        int group_flag         = msh->GetElementGroup(iel);      // element group flag (Exterior = GROUP_EXTERNAL, Interior = GROUP_INTERNAL)
        short unsigned ielGeom = msh->GetElementType(iel);       // element geometry type

//***************** GEOMETRY ************************
        unsigned nDofx = msh->GetElementDofNumber(iel, coords_fe_type); // number of coordinate element dofs
        for (int i = 0; i < dim; i++)  coords_at_dofs[i].resize(nDofx);
        // local storage of coordinates
        for (unsigned i = 0; i < nDofx; i++) {
            unsigned xDof  = msh->GetSolutionDof(i, iel, coords_fe_type); // global to global mapping between coordinates node and coordinate dof

            for (unsigned jdim = 0; jdim < dim; jdim++) {
                coords_at_dofs[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);  // global extraction and local storage for the element coordinates
            }
        }

        // elem average point
        vector < double > elem_center(dim);
        for (unsigned j = 0; j < dim; j++) {
            elem_center[j] = 0.;
        }
        for (unsigned j = 0; j < dim; j++) {
            for (unsigned i = 0; i < nDofx; i++) {
                elem_center[j] += coords_at_dofs[j][i];
            }
        }

        for (unsigned j = 0; j < dim; j++) {
            elem_center[j] = elem_center[j]/nDofx;
        }
//***************************************************

//************** set target domain flag *************
        int target_flag = 0;
        target_flag = ElementTargetFlag(elem_center);
//***************************************************


//**************** state ****************************
        unsigned nDof_u     = msh->GetElementDofNumber(iel, solType_u);     // number of solution element dofs
        sol_u    .resize(nDof_u);
        // local storage of global mapping and solution
        for (unsigned i = 0; i < sol_u.size(); i++) {
            unsigned solDof_u = msh->GetSolutionDof(i, iel, solType_u);        // global to global mapping between solution node and solution dof
            sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);              // global extraction and local storage for the solution
        }
//***************************************************


//************** control ****************************
        unsigned nDof_ctrl  = msh->GetElementDofNumber(iel, solType_ctrl);  // number of solution element dofs
        sol_ctrl    .resize(nDof_ctrl);
        for (unsigned i = 0; i < sol_ctrl.size(); i++) {
            unsigned solDof_ctrl = msh->GetSolutionDof(i, iel, solType_ctrl); // global to global mapping between solution node and solution dof
            sol_ctrl[i] = (*sol->_Sol[solIndex_ctrl])(solDof_ctrl);           // global extraction and local storage for the solution
        }
//***************************************************


//**************** u_des ****************************
        unsigned nDof_udes  = msh->GetElementDofNumber(iel, solType_u);    // number of solution element dofs
        sol_udes    .resize(nDof_udes);
        for (unsigned i = 0; i < sol_udes.size(); i++) {
            sol_udes[i] = u_des;  //dof value
        }
//***************************************************


//******************* ALL VARS **********************
        int nDof_max    =  nDof_u;   // TODO COMPUTE MAXIMUM maximum number of element dofs for one scalar variable

        if(nDof_udes > nDof_max)
        {
            nDof_max = nDof_udes;
        }

        if(nDof_ctrl > nDof_max)
        {
            nDof_max = nDof_ctrl;
        }

//***************************************************

        // *** Gauss point loop ***
        for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][coords_fe_type]->GetGaussPointNumber(); ig++) {

            // *** get gauss point weight, test function and test function partial derivatives ***
            //  ==== State
            msh->_finiteElement[ielGeom][solType_u]   ->Jacobian(coords_at_dofs, ig, weight_qp, phi_u, phi_u_x, phi_u_xx);
            //  ==== Adjoint
            msh->_finiteElement[ielGeom][solType_u/*solTypeTdes*/]->Jacobian(coords_at_dofs, ig, weight_qp, phi_udes, phi_udes_x, phi_udes_xx);
            //  ==== Control
            msh->_finiteElement[ielGeom][solType_ctrl]  ->Jacobian(coords_at_dofs, ig, weight_qp, phi_ctrl, phi_ctrl_x, phi_ctrl_xx);

            u_gss = 0.;
            for (unsigned i = 0; i < nDof_u; i++) u_gss += sol_u[i] * phi_u[i];
            ctrl_gss = 0.;
            for (unsigned i = 0; i < nDof_ctrl; i++) ctrl_gss += sol_ctrl[i] * phi_ctrl[i];
            udes_gss  = 0.;
            for (unsigned i = 0; i < nDof_udes; i++)  udes_gss  += sol_udes[i]  * phi_udes[i];
            ctrl_x_gss  = 0.;
            for (unsigned i = 0; i < nDof_ctrl; i++)
            {
                for (unsigned idim = 0; idim < dim; idim ++) ctrl_x_gss  += sol_ctrl[i] * phi_ctrl_x[i + idim * nDof_ctrl];
            }

            integral_target += target_flag * weight_qp * (u_gss - udes_gss) * (u_gss - udes_gss);
            integral_alpha  += (group_flag - GROUP_INTERNAL) * weight_qp * ctrl_gss * ctrl_gss;
            integral_beta   += (group_flag - GROUP_INTERNAL) * weight_qp * ctrl_x_gss * ctrl_x_gss;

        } // end gauss point loop
    } //end element loop

    std::cout << "The value of the integral_target is " << std::setw(11) << std::setprecision(10) << integral_target << std::endl;
    std::cout << "The value of the integral_alpha  is " << std::setw(11) << std::setprecision(10) << integral_alpha << std::endl;
    std::cout << "The value of the integral_beta   is " << std::setw(11) << std::setprecision(10) << integral_beta << std::endl;
    std::cout << "The value of the total integral  is " << std::setw(11) << std::setprecision(10) << 0.5 * integral_target + 0.5 * alpha * integral_alpha + 0.5 * beta * integral_beta << std::endl;

    return;

}




#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "LinearImplicitSystem.hpp"
#include "Parameter.hpp"
#include "paral.hpp"//to get iproc HAVE_MPI is inside here
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Files.hpp"
#include "FE_convergence.hpp"
#include "Assemble_jacobian.hpp"
#include "PetscMatrix.hpp"

#define FACE_FOR_CONTROL             3  //we do control on the right (=2) face
#define AXIS_DIRECTION_CONTROL_SIDE  0 //=0 if horizontal face, = 1 if vertical face (change this accordingly to the other variable above)


#include   "../solidopt_params.hpp"

#define CTRL_FACE_IDX  3

using namespace femus;


 double penalty_outside_control_boundary = 1.e50;       // penalty for zero control outside Gamma_c
 double penalty_ctrl = 1.e10;         //penalty for u=g
 double theta_value_outside_fake_element = 0.;


double Solution_set_initial_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char name[]) {
         
           double value = 0.;

             if(!strcmp(name,"DX")) {
                 value = 0.;
             }
             else if(!strcmp(name,"DY")) {
                 value = 0.;
             }
             else if(!strcmp(name,"DZ")) {
                 value = 0.;
             }
             else if(!strcmp(name,"P")) {
                 value = 0.;
             }

             else if(!strcmp(name,"DX_ADJ")) {
                 value = 0.;
             }
             else if(!strcmp(name,"DY_ADJ")) {
                 value = 0.;
             }
             else if(!strcmp(name,"DZ_ADJ")) {
                 value = 0.;
             }
             else if(!strcmp(name,"P_ADJ")) {
                 value = 0.;
             }
             
             else if(!strcmp(name,"GX")) {
                 value = 0.;
             }
             else if(!strcmp(name,"GY")) {
                 value = 0.;
             }
             else if(!strcmp(name,"GZ")) {
                 value = 0.;
             }
             else if(!strcmp(name,"THETA")) {
                 value = 0.;
             }
          
    
      return value;   
}  
  

bool Solution_set_boundary_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  //1: bottom (y=y_min) //2: right (x=x_max)  //3: top (y=y_max)  //4: left (x=x_min)  (in 2D)
  
  bool dirichlet  = true; 
  value = 0.;
  
 
                if (!strcmp(SolName, "GX"))       { if (facename == CTRL_FACE_IDX) dirichlet = false; }
                if (!strcmp(SolName, "GY"))       { if (facename == CTRL_FACE_IDX) dirichlet = false; }
                if (!strcmp(SolName, "GZ"))       { if (facename == CTRL_FACE_IDX) dirichlet = false; }
                if (!strcmp(SolName, "THETA"))    { dirichlet = false; }

               if (!strcmp(SolName, "DX"))       { if (facename == CTRL_FACE_IDX) dirichlet = false; }
                if (!strcmp(SolName, "DY"))       { if (facename == CTRL_FACE_IDX) dirichlet = false; }
                if (!strcmp(SolName, "DZ"))       { if (facename == CTRL_FACE_IDX) dirichlet = false; }
  
      
 return dirichlet;
}


template < class system_type, class real_num, class real_num_mov >
void AssembleSolidMech(MultiLevelProblem& ml_prob);


template < class system_type, class real_num, class real_num_mov >
void AssembleSolidMech(MultiLevelProblem& ml_prob,
                       const Solid & solid_in,
                       system_type * mlPdeSys,
                       const std::vector< Unknown > &  unknowns);



template < class system_type, class real_num, class real_num_mov >
void ComputeIntegral(const MultiLevelProblem& ml_prob);

template < class system_type, class real_num, class real_num_mov >
void ComputeIntegral(  const MultiLevelProblem& ml_prob,
                       const system_type * mlPdeSys,
                       const std::vector< Unknown > &  unknowns);


 //Unknown definition  ==================
 const std::vector< Unknown >  provide_list_of_unknowns(const unsigned int dimension) {
     
     
  std::vector< FEFamily > feFamily;
  std::vector< FEOrder >   feOrder;

                        feFamily.push_back(LAGRANGE);   //state
                        feFamily.push_back(LAGRANGE);
  if (dimension == 3)   feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/);
                        feFamily.push_back(LAGRANGE);   //adjoint
                        feFamily.push_back(LAGRANGE);
  if (dimension == 3)   feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/);
                        feFamily.push_back(LAGRANGE);   //control
                        feFamily.push_back(LAGRANGE);
  if (dimension == 3)   feFamily.push_back(LAGRANGE);
                        feFamily.push_back(DISCONTINUOUS_POLYNOMIAL);
 
                        feOrder.push_back(SECOND);
                        feOrder.push_back(SECOND);
  if (dimension == 3)   feOrder.push_back(SECOND);
                        feOrder.push_back(FIRST);
                        feOrder.push_back(SECOND);
                        feOrder.push_back(SECOND);
  if (dimension == 3)   feOrder.push_back(SECOND);
                        feOrder.push_back(FIRST);
                        feOrder.push_back(SECOND);
                        feOrder.push_back(SECOND);
  if (dimension == 3)   feOrder.push_back(SECOND);
                        feOrder.push_back(ZERO);
 

  assert( feFamily.size() == feOrder.size() );
 
 std::vector< Unknown >  unknowns(feFamily.size());
 
  const int adj_pos_begin   = dimension + 1;
  const int ctrl_pos_begin  = 2 * (dimension + 1);

                                        unknowns[0]._name      = "DX";
                                        unknowns[1]._name      = "DY";
  if (dimension == 3)                   unknowns[2]._name      = "DZ";
                                unknowns[dimension]._name      = "P";
                        unknowns[adj_pos_begin + 0]._name      = "DX_ADJ";
                        unknowns[adj_pos_begin + 1]._name      = "DY_ADJ";
  if (dimension == 3)   unknowns[adj_pos_begin + 2]._name      = "DZ_ADJ";
                unknowns[adj_pos_begin + dimension]._name      = "P_ADJ";
                       unknowns[ctrl_pos_begin + 0]._name      = "GX";
                       unknowns[ctrl_pos_begin + 1]._name      = "GY";
  if (dimension == 3)  unknowns[ctrl_pos_begin + 2]._name      = "GZ";
               unknowns[ctrl_pos_begin + dimension]._name      = "THETA";

     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
              unknowns[u]._fe_family  = feFamily[u];
              unknowns[u]._fe_order   = feOrder[u];
              unknowns[u]._time_order = 0;
              unknowns[u]._is_pde_unknown = true;
              
     }
 
 
   return unknowns;
     
}


template < class real_num > 
class My_main_single_level : public Main_single_level {
    
public:
    
const MultiLevelSolution  run_on_single_level(const Files & files,
                                              MultiLevelProblem & ml_prob,
                                              const std::vector< Unknown > & unknowns,
                                              const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition_in,
                                              const MultiLevelSolution::InitFuncMLProb SetInitialCondition_in,
                                              MultiLevelMesh & ml_mesh,
                                              const unsigned i) const;
  
};
 

void set_phi_gradient_surface(const unsigned int axis_direction_control_side, 
                      const Mesh* msh,
                      const unsigned int ig_bd, 
                      const std::vector<unsigned int> SolFEType, 
                      const unsigned int felt_bd, 
                      const unsigned int nve_bd,
                      const std::vector< std::vector< double > > & coords_hat_bd,
                      std::vector< std::vector< double > > & phi_hat_bd_x_dof_qp,
                      const unsigned int dim,  //dim of the (vector) unknown
                      const unsigned int ctrl_pos_begin //unknown
                     ) {
    
//========== temporary soln for surface gradient on a face parallel to the axes ===================
//           const unsigned int axis_direction_control_side = AXIS_DIRECTION_CONTROL_SIDE;
		  double dx_dcurv_abscissa = 0.;
		 const elem_type_1D * myeltype = static_cast<const elem_type_1D*>(msh->_finiteElement[felt_bd][SolFEType[ctrl_pos_begin]]);
		 const double * myptr = myptr = myeltype->GetDPhiDXi(ig_bd);
		      for (int inode = 0; inode < nve_bd; inode++) dx_dcurv_abscissa += myptr[inode] * coords_hat_bd[axis_direction_control_side][inode];
  
		    for (unsigned  kdim = 0; kdim < dim; kdim++) {
		      for (int inode = 0; inode < nve_bd; inode++) {
                            for (int d = 0; d < coords_hat_bd.size(); d++) {
                              if (d == axis_direction_control_side ) phi_hat_bd_x_dof_qp[kdim + ctrl_pos_begin][inode + d*nve_bd] = myptr[inode]* (1./ dx_dcurv_abscissa);
                              else  phi_hat_bd_x_dof_qp[kdim + ctrl_pos_begin][inode + d*nve_bd] = 0.;
                         }
                     }
            }
            
}


int main(int argc, char** args) {

  // ======= Init ==========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Files ==================
  Files files;
        files.CheckIODirectories();
        files.RedirectCout();

  // ======= Quad Rule ========================
  std::string quad_rule_order("seventh");

  // ======= Problem ========================
    MultiLevelProblem ml_prob;
  ml_prob.SetQuadratureRuleAllGeomElems(quad_rule_order);
  ml_prob.SetFilesHandler(&files);
    
  // ======= Mesh ==================
  const ElemType geom_elem_type = QUAD9;
  const std::vector< unsigned int > nsub = {NSUB, NSUB, 0};
  const std::vector< double >      xyz_min = {0., 0., 0.};
  const std::vector< double >      xyz_max = {1., 1., 0.};

  MultiLevelMesh ml_mesh;
  ml_mesh.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,quad_rule_order.c_str());

  // ======= Unknowns ========================
  const unsigned int dimension = ml_mesh.GetDimension();  
  std::vector< Unknown > unknowns = provide_list_of_unknowns(dimension);
  
  // ======= Normal run ========================   //if you don't want the convergence study
  My_main_single_level< adept::adouble > my_main;
  const unsigned int n_levels = 1;
  my_main.run_on_single_level(files, ml_prob, unknowns, Solution_set_boundary_conditions, Solution_set_initial_conditions, ml_mesh, n_levels); 
 
  
  
//   // ======= Convergence study ========================
//     
//    //set coarse storage mesh (should write the copy constructor or "=" operator to copy the previous mesh) ==================
//    MultiLevelMesh ml_mesh_all_levels;
//    ml_mesh_all_levels.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,quad_rule_order.c_str());
//    //   ml_mesh_all_levels.ReadCoarseMesh(infile.c_str(),quad_rule_order.c_str(),1.);
//  
//    // convergence choices ================  
//    unsigned int max_number_of_meshes;               // set total number of levels ================  
// 
//    if (nsub[2] == 0)   max_number_of_meshes = 3;
//    else                max_number_of_meshes = 4;
//   
// //    My_exact_solution<> exact_sol;                //provide exact solution, if available ==============
//    const unsigned conv_order_flag = 0;              //Choose how to compute the convergence order ========= //0: incremental 1: absolute (with analytical sol)  2: absolute (with projection of finest sol)...
//    const unsigned norm_flag = 1;                    //Choose what norms to compute (//0 = only L2: //1 = L2 + H1) ==============
// 
//    
//    // object ================  
//     FE_convergence<>  fe_convergence;
//     
//     fe_convergence.convergence_study(files, ml_prob, unknowns, Solution_set_boundary_conditions, Solution_set_initial_conditions, ml_mesh, ml_mesh_all_levels, max_number_of_meshes, norm_flag, conv_order_flag, my_main);
  
    
  return 0;
}



template < class system_type, class real_num, class real_num_mov = double >
void AssembleSolidMech(MultiLevelProblem& ml_prob) {
    

  AssembleSolidMech< system_type, real_num, real_num_mov > (  ml_prob, 
                                                              ml_prob.parameters.get<Solid>("Solid"), 
                                                            & ml_prob.get_system< system_type >(0), 
                                                              ml_prob.get_system< system_type >(0).get_unknown_list_for_assembly());

}



template < class system_type, class real_num, class real_num_mov = double >
void AssembleSolidMech(MultiLevelProblem& ml_prob,
                       const Solid & solid_in,
                       system_type * mlPdeSys,
                       const std::vector< Unknown > &  unknowns) {
    
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled


  const unsigned 		level 		    = mlPdeSys->GetLevelToAssemble();
  bool 			assembleMatrix 		    = mlPdeSys->GetAssembleMatrix(); 
  const char* 			pdename         = mlPdeSys->name().c_str();

  Mesh*          		msh    		= ml_prob._ml_msh->GetLevel(level);
  elem*          		el     		= msh->el;

  MultiLevelSolution*  		ml_sol  = ml_prob._ml_sol;
  Solution*    			sol      	= ml_prob._ml_sol->GetSolutionLevel(level); 


  LinearEquationSolver* 	pdeSys  = mlPdeSys->_LinSolver[level];
  SparseMatrix*    		JAC    		= pdeSys->_KK;
  NumericVector*   		RES         = pdeSys->_RES;

  MatSetOption(static_cast< PetscMatrix* >(JAC)->mat(),MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

  const unsigned 	dim  = msh->GetDimension();
  const unsigned   dim2  = (3 * (dim - 1) + !(dim - 1));
  const unsigned    nel  = msh->GetNumberOfElements();
  const unsigned  igrid  = msh->GetLevel();
  const unsigned  iproc  = msh->processor_id();
  // reserve memory for the local standard vectors
  const unsigned max_size_elem_dofs = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

 
  
  //----------- unknowns ------------------------------
   const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();
//      enum Sol_pos {pos_dx = 0, pos_dy, pos_p};  //these are known at compile-time
//      enum Sol_pos {pos_dx = 0, pos_dy, pos_dz, pos_p};
//   constexpr unsigned int pos_dx = 0;  
//   constexpr unsigned int pos_dy = 1;  
//   constexpr unsigned int pos_dp = 2;  
//   if (dim == 3) {std::cout << "Uncomment the 3d enum and recompile"; abort(); }
   
  constexpr int sol_index_displ      = 0;     //known at compile time
  const     int sol_index_press      = sol_index_displ + dim;      //known at run time
  constexpr int state_pos_begin      = sol_index_displ;   //known at compile time

  const     int adj_pos_begin        = dim + 1;
  const     int sol_index_adj        = adj_pos_begin;      
  const     int sol_index_press_adj  = sol_index_adj + dim;     

  const     int ctrl_pos_begin       = 2 * (dim + 1);
  const     int sol_index_ctrl       = ctrl_pos_begin;     
  const     int theta_index          = sol_index_ctrl + dim;      

  const int  n_state_vars = dim + 1;
  
  vector < std::string > Solname(n_unknowns);     
  vector < unsigned int > SolIndex(n_unknowns);  
  vector < unsigned int > SolPdeIndex(n_unknowns);
  vector < unsigned int > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
       Solname[ivar]    = unknowns[ivar]._name; 
      SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(Solname[ivar].c_str());
//       assert(ivar == SolPdeIndex[ivar]);  //I would like this ivar order to coincide either with SolIndex or with SolPdeIndex, otherwise too many orders... it has to be  with SolPdeIndex, because SolIndex is related to the Solution object which can have more fields... This has to match the order of the unknowns[] argument
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

   constexpr int solFEType_max = BIQUADR_FE;  //biquadratic, this is the highest-order FE 

   vector < unsigned int > Sol_n_el_dofs(n_unknowns);
  
  //----------- of dofs and at quadrature points ------------------------------
  vector < vector < double > > phi_dof_qp(n_unknowns);
  vector < vector < real_num_mov > > phi_x_dof_qp(n_unknowns);   //the derivatives depend on the Jacobian, which depends on the unknown, so we have to be adept here  //some of these should be real_num, some real_num_mov...
  vector < vector < real_num_mov > > phi_xx_dof_qp(n_unknowns);  //the derivatives depend on the Jacobian, which depends on the unknown, so we have to be adept here

  vector < vector < double > > phi_hat_dof_qp(n_unknowns);
  vector < vector < double > > phi_x_hat_dof_qp(n_unknowns);
  vector < vector < double > > phi_xx_hat_dof_qp(n_unknowns);

  for(int unk=0; unk < n_unknowns; unk++) {  
        phi_dof_qp[unk].reserve(max_size_elem_dofs);
    phi_hat_dof_qp[unk].reserve(max_size_elem_dofs);
      phi_x_dof_qp[unk].reserve(max_size_elem_dofs * dim);
  phi_x_hat_dof_qp[unk].reserve(max_size_elem_dofs * dim);
     phi_xx_dof_qp[unk].reserve(max_size_elem_dofs * dim2);
 phi_xx_hat_dof_qp[unk].reserve(max_size_elem_dofs * dim2);
   }
   
//----------boundary adjoint & ctrl shape functions----------------------  
  vector < vector < double > > phi_hat_bd_dof_qp(n_unknowns);
  vector < vector < double > > phi_hat_bd_x_dof_qp(n_unknowns);

  //bdry vol adj  evaluated at bdry points
   vector < vector < double > > phi_hat_vol_at_bdry_dof(n_unknowns);
   vector < vector < double > > phi_hat_x_vol_at_bdry_dof(n_unknowns);

    for(int unk=0; unk < n_unknowns ; unk++) {  
        phi_hat_bd_dof_qp[unk].reserve(max_size_elem_dofs);
      phi_hat_bd_x_dof_qp[unk].reserve(max_size_elem_dofs * dim);
  //bdry vol adj  evaluated at bdry points
         phi_hat_vol_at_bdry_dof[unk].reserve(max_size_elem_dofs);
       phi_hat_x_vol_at_bdry_dof[unk].reserve(max_size_elem_dofs * dim);     
    }

  vector < vector < real_num > >  sol_hat_adj_x_vol_at_bdry_qp(dim);
  vector < real_num > grad_dot_n_adj_res;
  vector < double > grad_adj_dot_n;
  for (int ldim =0; ldim < dim; ldim++) sol_hat_adj_x_vol_at_bdry_qp[ldim].reserve(max_size_elem_dofs);
  grad_dot_n_adj_res.reserve(max_size_elem_dofs);
  grad_adj_dot_n.reserve(max_size_elem_dofs);



    //----------- at dofs ------------------------------
  vector < double >   Jac;   Jac.reserve( n_unknowns * max_size_elem_dofs * n_unknowns * max_size_elem_dofs);
  vector < real_num > Res;   Res.reserve( n_unknowns * max_size_elem_dofs);
           vector < int >       L2G_dofmap_AllVars;   L2G_dofmap_AllVars.reserve( n_unknowns *max_size_elem_dofs);
  vector < vector < int > >         L2G_dofmap(n_unknowns);  for(int i = 0; i < n_unknowns; i++) {    L2G_dofmap[i].reserve(max_size_elem_dofs); }
  vector < vector < real_num > > SolVAR_eldofs(n_unknowns);  for(int k = 0; k < n_unknowns; k++) { SolVAR_eldofs[k].reserve(max_size_elem_dofs); }
 //=================state_only - for Jac retrieval via AD
  vector < double >   Jac_state_true;   Jac_state_true.reserve( n_state_vars * max_size_elem_dofs * n_state_vars * max_size_elem_dofs);
  vector < double >   Jac_state_false;   Jac_state_false.reserve( n_state_vars * max_size_elem_dofs * n_state_vars * max_size_elem_dofs);
  vector < real_num > Res_state;   Res_state.reserve( n_state_vars * max_size_elem_dofs);
 
  //----------------------------------------------------------------------------------------------------- 
  vector < real_num > Res_bd;   Res_bd.reserve( n_unknowns * max_size_elem_dofs);
  vector < double >   Jac_g;   Jac_g.reserve( n_unknowns * max_size_elem_dofs * n_unknowns * max_size_elem_dofs);
  vector < real_num > Res_g;   Res_g.reserve( n_unknowns * max_size_elem_dofs);
  vector < vector < double > > Jac_outer(dim);
  vector < vector < real_num > > Res_dctrl_dot_n(dim);
  vector < real_num > Res_outer(1);

         for(int i = 0; i < dim; i++) {  
             Jac_outer[i].reserve(max_size_elem_dofs); 
             Res_dctrl_dot_n[i].reserve(max_size_elem_dofs);
        }

  // geometry (at dofs) --------------------------------
  vector< vector < real_num_mov > >   coords(dim);
  vector  < vector  < double > >  coords_hat(dim);
  vector  < vector  < double > >  coords_hat_bd(dim);
  unsigned coordsType = BIQUADR_FE; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  for(int i = 0; i < dim; i++) {   
           coords[i].reserve(max_size_elem_dofs); 
       coords_hat[i].reserve(max_size_elem_dofs); 
    coords_hat_bd[i].reserve(max_size_elem_dofs); 
 }
 
 //************** geometry phi **************************
  vector < vector < double > > phi_dof_qp_domain(dim);
  vector < vector < real_num_mov > > phi_x_dof_qp_domain(dim);
  vector < vector < real_num_mov > > phi_xx_dof_qp_domain(dim);
  vector < vector < double > > phi_dof_qp_domain_hat(dim);
  vector < vector < double > > phi_x_dof_qp_domain_hat(dim);
  vector < vector < double > > phi_xx_dof_qp_domain_hat(dim);
 
  for(int dom_coord = 0; dom_coord < dim; dom_coord++) {
        phi_dof_qp_domain[dom_coord].reserve(max_size_elem_dofs);
      phi_x_dof_qp_domain[dom_coord].reserve(max_size_elem_dofs * dim);
     phi_xx_dof_qp_domain[dom_coord].reserve(max_size_elem_dofs * dim2);
        phi_dof_qp_domain_hat[dom_coord].reserve(max_size_elem_dofs);
      phi_x_dof_qp_domain_hat[dom_coord].reserve(max_size_elem_dofs * dim);
     phi_xx_dof_qp_domain_hat[dom_coord].reserve(max_size_elem_dofs * dim2);
   }
 
 
 
  // geometry ------------------------------------------
  

  //------------ at quadrature points ---------------------
  real_num_mov weight_qp = 0.;
    double weight_hat_qp = 0.;
    double weight_hat_bd_qp = 0.;
    vector < real_num > SolVAR_qp(n_unknowns);
    vector < real_num > SolVAR_hat_qp(n_unknowns);
    vector < vector < real_num > > gradSolVAR_qp(n_unknowns);
    vector < vector < real_num > > gradSolVAR_hat_qp(n_unknowns);
    for(int k = 0; k < n_unknowns; k++) { 
          gradSolVAR_qp[k].resize(dim);
      gradSolVAR_hat_qp[k].resize(dim);
    }
      

   // ------------------------------------------------------------------------
    // Physical parameters
    const int    solid_model	= solid_in.get_physical_model();
    double rhos         	 	= solid_in.get_density();
    const double mu_lame 	    = solid_in.get_lame_shear_modulus();
    const double lambda_lame 	= solid_in.get_lame_lambda();

    const bool incompressible   = (0.5  ==  solid_in.get_poisson_coeff()) ? 1 : 0;
    const bool penalty = solid_in.get_if_penalty();

    const double young_modulus = solid_in.get_young_module();
    const double Lref = solid_in._parameter->Get_reference_length();
    const double nondim_number = Lref/young_modulus;
    
    const double mus   =  /*nondim_number **/ mu_lame/rhos;
    // -----------------------------------------------------------------
 
    RES->zero();
  if (assembleMatrix) JAC->zero();
  
  adept::Stack & stack = FemusInit::_adeptStack;  // call the adept stack object for potential use of AD
  const assemble_jacobian< real_num, double/*real_num_mov*/ > assemble_jac;
 //=================state_only - for Jac retrieval via AD
  adept::Stack & stack_state = FemusInit::_adeptStack;  // call the adept stack object for potential use of AD
  const assemble_jacobian< real_num, double/*real_num_mov*/ > assemble_jac_state;


  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

   //all vars ###################################################################  
    unsigned int fake_iel_flag = 0;
    unsigned int global_row_index_bdry_constr = pdeSys->KKoffset[SolPdeIndex[theta_index]][iproc];
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
       Sol_n_el_dofs[k] = ndofs_unk;
       SolVAR_eldofs[k].resize(Sol_n_el_dofs[k]);
          L2G_dofmap[k].resize(Sol_n_el_dofs[k]); 
    for (unsigned i = 0; i < SolVAR_eldofs[k].size(); i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);
       SolVAR_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);
          L2G_dofmap[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
     if (k == SolPdeIndex[theta_index] && L2G_dofmap[k][i] == global_row_index_bdry_constr) {       fake_iel_flag = iel;  }
      }
    }
    
    L2G_dofmap_AllVars.resize(0);
      for (unsigned  k = 0; k < n_unknowns; k++)     L2G_dofmap_AllVars.insert(L2G_dofmap_AllVars.end(),L2G_dofmap[k].begin(),L2G_dofmap[k].end());

    unsigned sum_Sol_n_el_dofs = 0;
    for (unsigned  k = 0; k < n_unknowns; k++) { sum_Sol_n_el_dofs += Sol_n_el_dofs[k]; }
    Jac.resize(sum_Sol_n_el_dofs * sum_Sol_n_el_dofs);  std::fill(Jac.begin(), Jac.end(), 0.);
    Res.resize(sum_Sol_n_el_dofs);                      std::fill(Res.begin(), Res.end(), 0.);
   Res_bd.resize(sum_Sol_n_el_dofs);                      std::fill(Res_bd.begin(), Res_bd.end(), 0.);
    Jac_g.resize(sum_Sol_n_el_dofs * sum_Sol_n_el_dofs);  std::fill(Jac_g.begin(), Jac_g.end(), 0.);
    Res_g.resize(sum_Sol_n_el_dofs);                      std::fill(Res_g.begin(), Res_g.end(), 0.);
 
 //=================state_only - for Jac retrieval via AD
    unsigned state_sum_Sol_n_el_dofs = 0;
    for (unsigned  k = 0; k < n_state_vars; k++) { state_sum_Sol_n_el_dofs += Sol_n_el_dofs[k]; }
    Jac_state_true.resize(state_sum_Sol_n_el_dofs * state_sum_Sol_n_el_dofs);    std::fill(Jac_state_true.begin(), Jac_state_true.end(), 0.);
    Jac_state_false.resize(state_sum_Sol_n_el_dofs * state_sum_Sol_n_el_dofs);    std::fill(Jac_state_false.begin(), Jac_state_false.end(), 0.);
    Res_state.resize(state_sum_Sol_n_el_dofs);                              std::fill(Res_state.begin(), Res_state.end(), 0.);
 
     for(int ivar = 0; ivar < dim; ivar++)    {
        Jac_outer[ivar].resize(Sol_n_el_dofs[ivar + ctrl_pos_begin]);
        Res_dctrl_dot_n[ivar].resize(Sol_n_el_dofs[ivar + ctrl_pos_begin]);
    }
    Res_outer[0] = 0.;
    
    //all vars ###################################################################
  
  
//indices for restr --------------------------------------
           vector < int > ind_dtheta(theta_index - sol_index_press_adj);
            for (unsigned k = 0; k < ind_dtheta.size(); k++){
                ind_dtheta[k] = assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, SolPdeIndex[theta_index] , SolPdeIndex[ctrl_pos_begin + k], 0 , 0);
            }
    
           vector < int > ind_res_dctrl_dtheta(theta_index - sol_index_press_adj);
            for (unsigned k = 0; k < ind_res_dctrl_dtheta.size(); k++){
                ind_res_dctrl_dtheta[k] = assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[ctrl_pos_begin + k], 0);
            }
//indices for restr --------------------------------------
    
  // geometry *****************************
    short unsigned ielGeom = msh->GetElementType(iel);

      unsigned nDofsX = msh->GetElementDofNumber(iel, coordsType);    // number of coordinate element dofs
    
    for(int ivar=0; ivar<dim; ivar++) {
      coords[ivar].resize(nDofsX);
      coords_hat[ivar].resize(nDofsX);
    }
    
   for( unsigned i = 0; i < nDofsX; i++) {
      unsigned coordsDof  = msh->GetSolutionDof(i, iel, coordsType);    // global to global mapping between coordinates node and coordinate dof // via local to global solution node
      for(unsigned ivar = 0; ivar < dim; ivar++) {
          //Fixed coordinates (Reference frame)
	coords_hat[ivar][i] = (*msh->_topology->_Sol[ivar])(coordsDof);
      }
    }

  //***************************************  
    assert(nDofsX == Sol_n_el_dofs[sol_index_displ]);
    //Moving coordinates (Moving frame)
      for (unsigned idim = 0; idim < dim; idim++) {
        for (int j = 0; j < nDofsX; j++) {
          coords[idim][j] = coords_hat[idim][j] + SolVAR_eldofs[SolIndex[idim]][j];
        }
      }
  
     // elem average point 
    vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
  for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofsX; i++) {
         elem_center[j] += coords_hat[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofsX; }
  //*************************************** 
  
  //***** set target domain flag ********************************** 
   int target_flag = 0;
   target_flag = ElementTargetFlag < double > (elem_center);

  // geometry end *****************************

      
//************ set control flag *********************
    int control_el_flag = 0;
        control_el_flag = ControlDomainFlag(elem_center);
    std::vector< std::vector<int> > control_node_flag(dim);
	    for(unsigned idim=0; idim < dim; idim++) {
	          control_node_flag[idim].resize(Sol_n_el_dofs[ctrl_pos_begin]);
    std::fill(control_node_flag[idim].begin(), control_node_flag[idim].end(), 0);
	    }
 //*************************************************** 
      
   //************ set fake theta flag: this flag tells me what degrees of freedom of the current element are fake for the theta variable *********************
    std::vector<int>  bdry_int_constr_pos_vec(1,global_row_index_bdry_constr); /*KKoffset[SolPdeIndex[PADJ]][iproc]*/
    std::vector<int> fake_theta_flag(Sol_n_el_dofs[theta_index],0);
    for (unsigned i = 0; i < Sol_n_el_dofs[theta_index]; i++) {
      if ( L2G_dofmap[ SolPdeIndex[theta_index] ] [i] == bdry_int_constr_pos_vec[0]) { 	fake_theta_flag[i] = 1;       }
    }
 //************ end set fake theta flag *********************

     assemble_jac.prepare_before_integration_loop(stack);
     assemble_jac_state.prepare_before_integration_loop(stack_state); //=================state_only - for Jac retrieval via AD
     
  
//========BoundaryLoop=====================================================================

  // Perform face loop over elements that contain some control face
  if (control_el_flag == 1) {
	  
      double tau = 0.;
      vector < double > normal(dim,0);
	  
      for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
	  std::vector < double > xyz_bdc(3,0.);  //not being used, because the boundaries are identified by the face numbers
	  // look for boundary faces
	  if(el->GetFaceElementIndex(iel,jface) < 0) {
	      unsigned int face = -( msh->el->GetFaceElementIndex(iel,jface)+1);
	      if(  face == CTRL_FACE_IDX) { //control face
//=================================================== 
		   //we use the dirichlet flag to say: if dirichlet == true, we set 1 on the diagonal. if dirichlet == false, we put the boundary equation
		  std::vector < bool > dir_bool(dim);
		  for(unsigned idim=0; idim < dim; idim++) {
		      dir_bool[idim] = false; //mlSol->GetBdcFunction()(xyz_bdc,ctrl_name[idim].c_str(),tau,face,0.);
		  }
	  
//=================================================== 
		unsigned nve_bd = msh->GetElementFaceDofNumber(iel,jface, SolFEType[ctrl_pos_begin] ); //AAAAAAAAAAAAAAAAA
		const unsigned felt_bd = msh->GetElementFaceType(iel, jface);    
		for(unsigned i=0; i < nve_bd; i++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);
		    unsigned      iDof = msh->GetSolutionDof(i_vol, iel, coordsType);
		    for(unsigned idim=0; idim<dim; idim++) { coords_hat_bd[idim][i]=(*msh->_topology->_Sol[idim])(iDof); }
		}
		
//========= initialize gauss quantities on the boundary ============================================
		vector < real_num >                      SolVAR_bd_qp(n_unknowns);
		vector < vector < real_num > >       gradSolVAR_bd_qp(n_unknowns);
		for(int k=0; k<n_unknowns; k++) {  gradSolVAR_bd_qp[k].resize(dim);  }

//========= gauss_loop boundary===============================================================
		  for(unsigned ig_bd=0; ig_bd < ml_prob.GetQuadratureRule(felt_bd).GetGaussPointsNumber(); ig_bd++) {

              ml_prob._ml_msh->_finiteElement[felt_bd][SolFEType[theta_index]]->JacobianSur(coords_hat_bd,ig_bd,weight_hat_bd_qp,phi_hat_bd_dof_qp[theta_index],phi_hat_bd_x_dof_qp[theta_index],normal);
		    for (unsigned  kdim = 0; kdim < dim; kdim++) {
		      ml_prob._ml_msh->_finiteElement[felt_bd][SolFEType[kdim + ctrl_pos_begin]]->JacobianSur(coords_hat_bd,ig_bd,weight_hat_bd_qp,phi_hat_bd_dof_qp[kdim + ctrl_pos_begin],phi_hat_bd_x_dof_qp[kdim + ctrl_pos_begin],normal);
		      ml_prob._ml_msh->_finiteElement[ielGeom][SolFEType[kdim + adj_pos_begin]]->fill_volume_shape_funcs_at_boundary_quadrature_points_on_current_elem(coords_hat,coords_hat_bd,jface,ig_bd,phi_hat_vol_at_bdry_dof[kdim + adj_pos_begin],phi_hat_x_vol_at_bdry_dof[kdim + adj_pos_begin]);
            }
            
            
 set_phi_gradient_surface(AXIS_DIRECTION_CONTROL_SIDE, msh, ig_bd, SolFEType, felt_bd, nve_bd, coords_hat_bd, phi_hat_bd_x_dof_qp, dim, ctrl_pos_begin);
 

//========== compute gauss quantities on the boundary ===============================================
		    for (unsigned  kdim = 0; kdim < dim; kdim++) {
			    unsigned int ctrl_index = kdim + ctrl_pos_begin;
										  SolVAR_bd_qp[ SolPdeIndex[ctrl_index] ] = 0.;
			  for(unsigned ivar2=0; ivar2<dim; ivar2++) {  gradSolVAR_bd_qp[ SolPdeIndex[ctrl_index] ][ivar2] = 0.; }
            }
	  
		    for (unsigned  kdim = 0; kdim < dim; kdim++) {
			    unsigned int ctrl_index = kdim + ctrl_pos_begin;
			  for(int i_bd = 0; i_bd < nve_bd; i_bd++) {
		                  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
                                                                    SolVAR_bd_qp[SolPdeIndex[ctrl_index]]           += phi_hat_bd_dof_qp  [ ctrl_index ][i_bd]                  * SolVAR_eldofs[ SolPdeIndex[ ctrl_index ] ][i_vol];
			      for(unsigned ivar2=0; ivar2<dim; ivar2++) {  gradSolVAR_bd_qp[SolPdeIndex[ctrl_index]][ivar2] 	+= phi_hat_bd_x_dof_qp[ ctrl_index ][i_bd + ivar2 * nve_bd] * SolVAR_eldofs[ SolPdeIndex[ ctrl_index ] ][i_vol];      }
			  }
		    }//kdim
 //end unknowns eval at gauss points ********************************

 
//=============== grad dot n for residual ========================================= 
//     compute gauss quantities on the boundary through VOLUME interpolation
		for(unsigned ldim=0; ldim<dim; ldim++) {   sol_hat_adj_x_vol_at_bdry_qp[ldim].resize(dim); }
		grad_dot_n_adj_res.resize(dim);
		grad_adj_dot_n.resize(dim);
		for(unsigned ldim=0; ldim<dim; ldim++) {   std::fill(sol_hat_adj_x_vol_at_bdry_qp[ldim].begin(), sol_hat_adj_x_vol_at_bdry_qp[ldim].end(), 0.);  }
		
		for(unsigned ldim=0; ldim<dim; ldim++) {  
		  for (int iv = 0; iv < Sol_n_el_dofs[adj_pos_begin]; iv++)  {
                     for (int d = 0; d < dim; d++) {
			   sol_hat_adj_x_vol_at_bdry_qp[ldim][d] += SolVAR_eldofs[SolPdeIndex[ldim + adj_pos_begin]][iv] * phi_hat_x_vol_at_bdry_dof[ldim + adj_pos_begin][iv * dim + d];//notice that the convention of the orders x y z is different from vol to bdry
		     }
		  }  
		      
		  grad_dot_n_adj_res[ldim] = 0.;
		  for(unsigned d=0; d<dim; d++) {
		      grad_dot_n_adj_res[ldim] += sol_hat_adj_x_vol_at_bdry_qp[ldim][d]*normal[d];  
		  }
		}
//=============== grad dot n  for residual =========================================       
		  
//========== compute gauss quantities on the boundary ================================================

		
  // *** phi_i loop ***
		for(unsigned i_bdry=0; i_bdry < nve_bd; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
                 
//=============== construct control node flag field on the go  =========================================    
	      /* (control_node_flag[i])       picks nodes on \Gamma_c
	         (1 - control_node_flag[i])   picks nodes on \Omega \setminus \Gamma_c
	       */
		    for(unsigned idim=0; idim<dim; idim++) {
			if (dir_bool[idim] == false) { 
// 		std::cout << " found boundary control nodes ==== " << std::endl;
			    for(unsigned k=0; k < control_node_flag[idim].size(); k++) {
				  control_node_flag[idim][i_vol] = 1;
			    }
			}
		    }
//=============== construct control node flag field on the go  =========================================    
		
		  
//Boundary Residuals  and Jacobians ==================	

		  
//============ Boundary Residuals============================================================================================
		  
		      for (unsigned  kdim = 0; kdim < dim; kdim++) {
			
			    real_num lap_res_dctrl_ctrl_bd = 0.;
			      for (unsigned jdim = 0; jdim < dim; jdim++) {
				  lap_res_dctrl_ctrl_bd += gradSolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]][jdim]*phi_hat_bd_x_dof_qp[kdim + ctrl_pos_begin][i_bdry + jdim*nve_bd];
			      }//jdim

			
/*delta_state row */	    if(i_vol<Sol_n_el_dofs[state_pos_begin])   Res_bd[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[kdim], i_vol) ]                += control_node_flag[kdim][i_vol] * penalty_ctrl * (SolVAR_eldofs[SolPdeIndex[kdim + state_pos_begin]][i_vol] - SolVAR_eldofs[SolPdeIndex[kdim + ctrl_pos_begin]][i_vol]);	    //u-g
/*delta_adjoint row */     if(i_vol<Sol_n_el_dofs[adj_pos_begin])   Res_bd[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[kdim + adj_pos_begin], i_vol) ] += 0.;	   
/*delta_control row */     if(i_vol<Sol_n_el_dofs[ctrl_pos_begin])  Res_bd[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[kdim + ctrl_pos_begin], i_vol) ]  += control_node_flag[kdim][i_vol] * weight_hat_bd_qp * (
                                                                                          beta_val* SolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * phi_hat_bd_dof_qp[kdim +  ctrl_pos_begin][i_bdry]
                                                                                        + gamma_val* lap_res_dctrl_ctrl_bd
                                                                                        - nondim_number * mus * grad_dot_n_adj_res[kdim]  * phi_hat_bd_dof_qp[kdim +  ctrl_pos_begin][i_bdry]
                                                                                        );	    
		      }//kdim  

//============ Boundary Residuals  ==================================================================================================


                        }//end i_bdry loop

                    }  //end ig_bdry loop
	  
                }    //end if control face
            }  //end if boundary faces
        }  // loop over element faces //jface   
  } //end if control element flag

//End Boundary Residuals  and Jacobians ==================	


//======================= Loop without Integration =====================================================    
        //============ delta_theta - theta row ==================================================================================================
  for (unsigned i = 0; i < Sol_n_el_dofs[theta_index]; i++) {
             /* if ( fake_theta_flag[i] != 1 ) */       Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[theta_index], i) ]    = (1 - fake_theta_flag[i]) * ( theta_value_outside_fake_element - SolVAR_eldofs[SolPdeIndex[theta_index]][i]);  // Res_outer for the exact row (i.e. when fakeflag=1 , res =0(use Res_outer) and if not 1 this loop) and this is to take care of fake placement for the rest of dofs of theta values as 8
        }//i_theta loop
   
 //============ delta_theta row ==================================================================================================
 //======================= Loop without Integration =====================================================    

 
 
 
//======================= VolumeLoop with Integration (and fake boundary) =====================================================    
// *** Gauss point loop ***
    for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {

	// *** get Jacobian and test function and test function derivatives ***
      for(int unk=0; unk < n_unknowns; unk++) {
	msh->_finiteElement[ielGeom][SolFEType[unk]]->Jacobian(coords,    ig,weight_qp,    phi_dof_qp[unk],    phi_x_dof_qp[unk],    phi_xx_dof_qp[unk]);
	msh->_finiteElement[ielGeom][SolFEType[unk]]->Jacobian(coords_hat,ig,weight_hat_qp,phi_hat_dof_qp[unk],phi_x_hat_dof_qp[unk],phi_xx_hat_dof_qp[unk]);
      }
      
   //=============== Integration at qp ========================================
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN...
    for (unsigned int d = 0; d < dim; d++) {
         msh->_finiteElement[ielGeom][coordsType]->Jacobian(coords,    ig, weight_qp,    phi_dof_qp_domain[d],     phi_x_dof_qp_domain[d],     phi_xx_dof_qp_domain[d]);
         msh->_finiteElement[ielGeom][coordsType]->Jacobian(coords_hat,ig, weight_hat_qp,phi_dof_qp_domain_hat[d], phi_x_dof_qp_domain_hat[d], phi_xx_dof_qp_domain_hat[d]);
      }
     
  //begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk <  n_unknowns; unk++) {
	  SolVAR_qp[SolPdeIndex[unk]] = 0.;
	  SolVAR_hat_qp[SolPdeIndex[unk] ]= 0.;

	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR_qp[SolPdeIndex[unk]][ivar2] = 0.; 
	    gradSolVAR_hat_qp[SolPdeIndex[unk]][ivar2] = 0.; 
	  }
    }
	  
	for(unsigned unk = 0; unk <  n_unknowns; unk++) {
	  for(unsigned i = 0; i < Sol_n_el_dofs[unk]; i++) {
	    SolVAR_qp[SolPdeIndex[unk]]    += phi_dof_qp[ unk ][i] * SolVAR_eldofs[SolPdeIndex[unk]][i];
	    SolVAR_hat_qp[SolPdeIndex[unk]]    += phi_hat_dof_qp[ unk ][i] * SolVAR_eldofs[SolPdeIndex[unk]][i];
	    for(unsigned ivar2 = 0; ivar2 < dim; ivar2++) {
	      gradSolVAR_qp[SolPdeIndex[unk]][ivar2]     += phi_x_dof_qp[ unk ][i*dim+ivar2] * SolVAR_eldofs[SolPdeIndex[unk]][i]; 
	      gradSolVAR_hat_qp[SolPdeIndex[unk]][ivar2] += phi_x_hat_dof_qp[ unk ][i*dim+ivar2] * SolVAR_eldofs[SolPdeIndex[unk]][i]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************


//  assemble_jacobian< real_num, real_num_mov >::mass_residual (Res, Sol_n_el_dofs, sum_Sol_n_el_dofs, SolPdeIndex, phi_dof_qp, SolVAR_qp, weight_hat_qp); //identity


 //*******************************************************************************************************
   vector < vector < real_num_mov > > Cauchy(3); for (int i = 0; i < Cauchy.size(); i++) Cauchy[i].resize(3);
   real_num_mov J_hat;
   real_num_mov trace_e_hat;

 
 Cauchy = Solid::get_Cauchy_stress_tensor< real_num_mov >(solid_model, mus, lambda_lame, incompressible, dim, sol_index_displ, sol_index_press, gradSolVAR_hat_qp, SolVAR_qp, SolPdeIndex, J_hat, trace_e_hat);

    
  //STATE=====================================
//----------from displ_state and p_state
              //BEGIN residual Solid Momentum in moving domain
          for (unsigned i = 0; i < Sol_n_el_dofs[sol_index_displ]; i++) {

              real_num_mov Cauchy_direction[3] = {0., 0., 0.};

              for (int idim = 0.; idim < dim; idim++) {
                for (int jdim = 0.; jdim < dim; jdim++) {
                  Cauchy_direction[idim] += nondim_number * phi_x_dof_qp[ idim ][i * dim + jdim] * Cauchy[idim][jdim];
                }
              }

              for (int idim = 0; idim < dim; idim++) {
                Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[idim], i) ] += ( Cauchy_direction[idim] -  phi_dof_qp[ idim ][i] * _gravity[idim] ) * weight_qp;
              }

            }
              //END residual Solid Momentum in moving domain
//----------from displ_state and p_state
              

//----------from displ_state of divergence
              //BEGIN residual solid mass balance in reference domain
            for (unsigned i = 0; i < Sol_n_el_dofs[sol_index_press]; i++) {
                
              Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[sol_index_press], i) ] += 
                - weight_hat_qp * phi_hat_dof_qp[ sol_index_press ][i] * Solid::get_mass_balance_reference_domain< real_num_mov >(solid_model, penalty, incompressible, lambda_lame, trace_e_hat, J_hat, SolVAR_qp, SolPdeIndex, sol_index_press);
//               - weight_hat_qp * phi_dof_qp[ sol_index_press ][i] * Solid::get_mass_balance_moving_domain< real_num_mov >(gradSolVAR_qp, SolPdeIndex);
              
            }
              //END residual solid mass balance in reference domain
//----------from displ_state of div
  //STATE=====================================

 
 //ADJOINT=====================================
//----------from displ_state and displ_ctrl
         for (unsigned i = 0; i < Sol_n_el_dofs[sol_index_adj]; i++) {
               for (int idim = 0; idim < dim; idim++) {
   Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[idim + adj_pos_begin], i)] += 
    - alpha_val * target_flag * (SolVAR_hat_qp[SolPdeIndex[idim]] - TargetDisp[idim]) 
     * phi_dof_qp[ idim + adj_pos_begin ][i] * weight_qp;
               }
          }
//----------from displ_state and displ_ctrl
 //ADJOINT=====================================
 
               
 //CONTROL=====================================
 //============ delta_control row ==================================================================================================
// delta_control
    for (unsigned kdim = 0; kdim < dim; kdim++) { 
         for (unsigned i = 0; i < Sol_n_el_dofs[ctrl_pos_begin]; i++) {
       Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[kdim + ctrl_pos_begin], i) ] += - penalty_outside_control_boundary * ( (1 - control_node_flag[kdim][i]) * (  SolVAR_eldofs[SolPdeIndex[kdim + ctrl_pos_begin]][i] - 0.)  );              //enforce control zero outside the control boundary
        }//i_ctrl loop
    }  //kdim

 //============ delta_control row ==================================================================================================
 //CONTROL=====================================

        
        
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------

    
 for (int i = 0; i < Res_state.size(); i++) {
        Res_state[i] = Res[i];
  }
  
    std::vector < double > Res_state_double(Res_state.size());

    for (int i = 0; i < Res_state_double.size(); i++) {
      Res_state_double[i] = - Res_state[i].value();
    }


    stack_state.dependent(  & Res_state[0], Res_state_double.size());      // define the dependent variables
    for (unsigned  k = 0; k < n_state_vars; k++) {   
    stack_state.independent(&SolVAR_eldofs[k][0],       SolVAR_eldofs[k].size());    // define the independent variables
    }
    
    stack_state.jacobian(&Jac_state_true[0], true);    // get the jacobian matrix (ordered by row major for Ctrl block in delta_state) 
    stack_state.jacobian(&Jac_state_false[0], false);    // get the jacobian matrix (ordered by column major for Adjoint) // Adj = transpose of state block

    stack_state.clear_independents();
    stack_state.clear_dependents();    

 
  //ADJOINT=====================================
//----------from adj displ
      for(unsigned i_block = adj_pos_begin; i_block < adj_pos_begin + dim; i_block++) {
        for(unsigned i_dof = 0; i_dof < Sol_n_el_dofs[i_block]; i_dof++) {
            for(unsigned j_block = 0; j_block < dim; j_block++) {
               for(unsigned j_dof = 0; j_dof < Sol_n_el_dofs[j_block]; j_dof++) {
    Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, i_block, i_dof)] += 
        Jac_state_false[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, state_sum_Sol_n_el_dofs,i_block - adj_pos_begin , j_block, i_dof, j_dof )] 
            * SolVAR_eldofs[SolPdeIndex[j_block + adj_pos_begin]][j_dof];
                }
           }
        }
   }
//----------from adj displ 

//----------from p_adj only
   for(unsigned i_block = adj_pos_begin; i_block < adj_pos_begin + dim ; i_block++) {
        for(unsigned i_dof = 0; i_dof < Sol_n_el_dofs[i_block]; i_dof++) {
            for(unsigned j_block = dim; j_block < n_state_vars; j_block++) {
               for(unsigned j_dof = 0; j_dof < Sol_n_el_dofs[j_block]; j_dof++) {
    Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, i_block, i_dof)] += 
        Jac_state_true[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, state_sum_Sol_n_el_dofs,i_block - adj_pos_begin , j_block, i_dof, j_dof )] 
            * SolVAR_eldofs[SolPdeIndex[j_block + adj_pos_begin]][j_dof];
                }
           }
        }
   }
//----------from p_adj only

//----------from div of displ_adj
   for(unsigned i_block = sol_index_press_adj; i_block < ctrl_pos_begin ; i_block++) {
        for(unsigned i_dof = 0; i_dof < Sol_n_el_dofs[i_block]; i_dof++) {
            for(unsigned j_block = 0; j_block < dim; j_block++) {
               for(unsigned j_dof = 0; j_dof < Sol_n_el_dofs[j_block]; j_dof++) {
    Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, i_block, i_dof)] += 
        Jac_state_true[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, state_sum_Sol_n_el_dofs,i_block - adj_pos_begin , j_block, i_dof, j_dof )] 
            * SolVAR_eldofs[SolPdeIndex[j_block + adj_pos_begin]][j_dof];
                }
           }
        }
   }
// ----------from div of displ_adj
  //ADJOINT=====================================

 
    for(unsigned i_block = 0; i_block < n_unknowns; i_block++) {
            for(unsigned i_dof = 0; i_dof < Sol_n_el_dofs[i_block]; i_dof++) {
                Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, i_block, i_dof)] += 
                Res_bd[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, i_block, i_dof)];
            }
    }
 
  
      
     if (assembleMatrix) assemble_jac.compute_jacobian_outside_integration_loop(stack, SolVAR_eldofs, Res, Jac, L2G_dofmap_AllVars, RES, JAC);


//       //***************************************************************************************************************
     
  if (control_el_flag == 1) {
      
      vector < double > normal(dim,0);

        for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
                
           	  if(el->GetFaceElementIndex(iel,jface) < 0) {  //I am on the boundary

                    unsigned int face = -( msh->el->GetFaceElementIndex(iel,jface)+1);
                        if(  face == CTRL_FACE_IDX) { //I am a control face
              
		unsigned nve_bd = msh->GetElementFaceDofNumber(iel,jface, SolFEType[ctrl_pos_begin] ); //AAAAAAAAAAAAAAAAA
		const unsigned felt_bd = msh->GetElementFaceType(iel, jface);    

		vector < real_num >                      SolVAR_bd_qp(n_unknowns);

        for(unsigned ig_bd=0; ig_bd < ml_prob.GetQuadratureRule(felt_bd).GetGaussPointsNumber(); ig_bd++) {
              
              		    for (unsigned  kdim = 0; kdim < dim; kdim++) {
		      ml_prob._ml_msh->_finiteElement[felt_bd][SolFEType[kdim + ctrl_pos_begin]]->JacobianSur(coords_hat_bd,ig_bd,weight_hat_bd_qp,phi_hat_bd_dof_qp[kdim + ctrl_pos_begin],phi_hat_bd_x_dof_qp[kdim + ctrl_pos_begin],normal);
                        }
              
              
              		    for (unsigned  kdim = 0; kdim < dim; kdim++) {
			    unsigned int ctrl_index = kdim + ctrl_pos_begin;
										  SolVAR_bd_qp[ SolPdeIndex[ctrl_index] ] = 0.;
                                          
                        for(int i_bd = 0; i_bd < nve_bd; i_bd++) {
		                  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
                                SolVAR_bd_qp[SolPdeIndex[ctrl_index]] += phi_hat_bd_dof_qp  [ ctrl_index ][i_bd]  * SolVAR_eldofs[ SolPdeIndex[ ctrl_index ] ][i_vol];
			                  }                                         
                        }
              
              
  //============ Res _ Boundary Integral Constraint ============================================================================================
              //residual  -  delta_theta  * g dot n =========================================
    for (unsigned  kdim = 0; kdim < dim; kdim++) {
// 		for(unsigned i=0; i < Sol_n_el_dofs[theta_index]; i ++) { //avoid because it is an element dof
/*delta_theta row */ 	 /*Res_bd[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[theta_index], i) ] */Res_outer[0] += - /*fake_theta_flag[i] **/ weight_hat_bd_qp * SolVAR_bd_qp[SolPdeIndex[kdim + ctrl_pos_begin]] * normal[kdim] ;
// 		}  
	  }
              //residual  -  delta_theta * g dot n =========================================
              
              		for(unsigned i_bdry=0; i_bdry < nve_bd; i_bdry++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bdry);
          
              //residual  -  theta * delta_g dot n =========================================
              		      for (unsigned  kdim = 0; kdim < dim; kdim++) {
			
			
/*delta_control row */     if( i_vol < Sol_n_el_dofs[ctrl_pos_begin] )  Res_g[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[kdim + ctrl_pos_begin], i_vol) ]  += -
                                                                                       control_node_flag[kdim][i_vol] * weight_hat_bd_qp * (
                                                                                         SolVAR_eldofs[SolPdeIndex[theta_index]][0] /*(*sol->_Sol[SolIndex[theta_index]])(0)*/ * phi_hat_bd_dof_qp[kdim +  ctrl_pos_begin][i_bdry] * normal[kdim]  
                                                                                        //*sol->_Sol[SolIndex[theta_index]])(0) finds the global value from KKDof pos(63, 169,etc), SolVAReldof_theta gets the value in the boundary point which will be zero. Theta is just a const
                                                                                        );	    
		                        }//kdim  
              //residual  -  theta * delta_g dot n =========================================
              //residual - end =========================================
//============ End of Res _ Boundary Integral Constraint ============================================================================================
              
              //jacobian =========================================
              
  //============ Jac _ Boundary Integral Constraint ============================================================================================
		    for (unsigned  kdim = 0; kdim < dim; kdim++) { 
			  for(unsigned i =0; i < Sol_n_el_dofs[theta_index]; i ++) {
			    if(i_vol < Sol_n_el_dofs[ctrl_pos_begin]) {
				double temp = weight_hat_bd_qp * ( phi_hat_bd_dof_qp[kdim + ctrl_pos_begin][i_bdry] * normal[kdim]);
//ROW_BLOCK delta_theta - control -- loop over i in the VOLUME (while j(/i_vol) is in the boundary) -------------------------------------------------------------------------------------------------------------
			      Jac_g[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, theta_index , ctrl_pos_begin + kdim, i, i_vol) ]  +=  temp;
//COLUMN_BLOCK delta_control - theta ---- loop over j in the VOLUME (while i(/i_vol) is in the boundary) ---------------------------------------------------------------------------------------------------
			      Jac_g[ assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, sum_Sol_n_el_dofs, ctrl_pos_begin + kdim, theta_index , i_vol, i) ] +=  control_node_flag[kdim][i_vol] * temp;
			    }//endif
			  }// i 
		    }//kdim
//============ End of Jac _ Boundary Integral Constraint ============================================================================================
            
              
              
              //jacobian - end =========================================
              
                    } //i_bdry
                    
          }  //end gauss boundary loop

              
     for(unsigned kdim = 0; kdim < dim; kdim++){
          vector < double > Jac_dtheta_g(Jac_g.begin() + ind_dtheta[kdim], Jac_g.begin() + ind_dtheta[kdim + 1]);
          vector < real_num >  Res_dctrl_dot_n_theta(Res_g.begin() + ind_res_dctrl_dtheta[kdim], Res_g.begin() + ind_res_dctrl_dtheta[kdim + 1]);
          Jac_outer[kdim] = Jac_dtheta_g;
          Res_dctrl_dot_n[kdim] = Res_dctrl_dot_n_theta;
     }
             
                        }
                    }
                }      
            }
     

     vector < double > Res_outer_double(1,0.);
    for (int i=0; i < Res_outer_double.size(); i++){   Res_outer_double[i] = Res_outer[i].value();   }

     vector < vector < double > > Res_dctrl_dot_n_double(dim);
     for (unsigned kdim =0; kdim < dim; kdim++) Res_dctrl_dot_n_double[kdim].resize(Sol_n_el_dofs[ctrl_pos_begin + kdim]);
    for (unsigned kdim =0; kdim < dim; kdim++) {
        for (int i=0; i < Res_dctrl_dot_n_double[kdim].size(); i++){
            Res_dctrl_dot_n_double[kdim][i] = Res_dctrl_dot_n[kdim][i].value();
        }
     }

   //--------------------------------------------------------------------------------------------------------  
    // THE BLOCKS WITH THETA ROW OR COLUMN 
// 	/*delta_theta-theta*/    JAC->add_matrix_blocked( Jac[ SolPdeIndex[n_unknowns-1] ][ SolPdeIndex[n_unknowns-1] ], L2G_dofmap[n_unknowns-1], L2G_dofmap[n_unknowns-1]);
	    
     if (control_el_flag == 1) {
	      for (unsigned kdim = 0; kdim < dim; kdim++) {
                          /*delta_control*/       RES->add_vector_blocked(Res_dctrl_dot_n_double[kdim],L2G_dofmap[n_unknowns-2-kdim]); 
		if(assembleMatrix) {
                          /*delta_theta-control*/ JAC->add_matrix_blocked( Jac_outer[kdim], bdry_int_constr_pos_vec, L2G_dofmap[n_unknowns-2-kdim]);
                          /*delta_control-theta*/ JAC->add_matrix_blocked( Jac_outer[kdim], L2G_dofmap[n_unknowns-2-kdim], bdry_int_constr_pos_vec); 
		}
	      }  //kdim
     }  //add control boundary element contributions
     
     
          if (control_el_flag == 1) {
          /*delta_theta(bdry constr)*/         RES->add_vector_blocked(Res_outer_double,bdry_int_constr_pos_vec);
	  }
	  
//      /* if (L2G_dofmap[n_unknowns-1][0] != bdry_int_constr_pos_vec[0]) */ /*delta_theta(fake)*/          RES->add_vector_blocked( Res[ SolPdeIndex[n_unknowns-1]],       L2G_dofmap[n_unknowns-1]);
	  
   //--------------------------------------------------------------------------------------------------------  

       
//     assemble_jacobian<real_num,real_num_mov>::print_element_residual(iel, Res, Sol_n_el_dofs, 9, 3);
//     assemble_jacobian<real_num,real_num_mov>::print_element_jacobian(iel, Jac, Sol_n_el_dofs, 4, 1);
    
  } //end element loop for each process


 if (assembleMatrix) JAC->close();
                     RES->close();
  
  
  
    //print JAC and RES to files
//     const unsigned nonlin_iter = mlPdeSys->GetNonlinearIt();
//     assemble_jacobian< real_num_mov,double >::print_global_jacobian(assembleMatrix, ml_prob, JAC, nonlin_iter);
//     assemble_jacobian< real_num_mov,double >::print_global_residual(ml_prob, RES, nonlin_iter);

  
}




template < class real_num > 
const MultiLevelSolution  My_main_single_level< real_num >::run_on_single_level(const Files & files,
                                                                                MultiLevelProblem & ml_prob,
                                                                                const std::vector< Unknown > &  unknowns,
                                                                                const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition_in,
                                                                                const MultiLevelSolution::InitFuncMLProb SetInitialCondition_in,
                                                                                MultiLevelMesh & ml_mesh,
                                                                                const unsigned lev) const {
                                                                                    
                                                                                    
  // ======= Mesh  ==================
            unsigned numberOfUniformLevels = lev + 1;  //this has to have a + 1 for the convergence study
            unsigned numberOfSelectiveLevels = 0;
            ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
            ml_mesh.EraseCoarseLevels(numberOfUniformLevels - 1);

            ml_mesh.PrintInfo();
                                                                                    
  // ======= Solution  ==================
  MultiLevelSolution ml_sol(&ml_mesh);

  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);
  
  // ======= Problem ========================
  ml_prob.SetMultiLevelMeshAndSolution(& ml_mesh,& ml_sol);
  
  ml_prob.get_systems_map().clear();

  //material  ==================
              //Adimensional quantity (Lref,Uref)
            double Lref = 1.;
            double Uref = 1.;
           // *** apparently needed by non-AD assemble only **********************
            // add fluid material
            Parameter par(Lref,Uref);
            
            // Generate Solid Object
            double E = YOUNGS_MODULUS;
            double ni = NI;
            double rhos = SOLID_DENSITY;
            Solid solid;
            solid = Solid(par,E,ni,rhos,MODEL);

            std::cout << "Solid properties: " << std::endl;
            std::cout << solid << std::endl;
            
            
            ml_prob.parameters.set<Solid>("Solid") = solid;
  //end material  ==================

  // ======= Solution, II ==================

  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryCondition_in);
  for (unsigned int u = 0; u < unknowns.size(); u++)  { 
      ml_sol.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order, unknowns[u]._time_order, unknowns[u]._is_pde_unknown);
      ml_sol.Initialize(unknowns[u]._name.c_str(), SetInitialCondition_in, &ml_prob);
      ml_sol.GenerateBdc(unknowns[u]._name.c_str(), "Steady", & ml_prob);
  }
  
  // ======= System ========================
  NonLinearImplicitSystem& system = ml_prob.add_system < NonLinearImplicitSystem > ("SolidMech");

  for (unsigned int u = 0; u < unknowns.size(); u++)   system.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
 
  system.set_unknown_list_for_assembly(unknowns); 
            
  system.SetAssembleFunction( AssembleSolidMech< NonLinearImplicitSystem, adept::adouble, adept::adouble >);

  // initialize and solve the system
  system.init();
  
  // Solver and preconditioner
  system.SetOuterSolver(PREONLY);
  system.SetSolverFineGrids(PREONLY);
  system.SetPreconditionerFineGrids(LU_PRECOND);

  //for Vanka   
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");

  system.SetDebugNonlinear(true);
  system.SetDebugFunction( ComputeIntegral< NonLinearImplicitSystem, adept::adouble, adept::adouble >);

//   system.SetMaxNumberOfNonLinearIterations(3);
//   system.SetNonLinearConvergenceTolerance(1.e-2);
//   system.SetDebugLinear(true);
//   system.SetMaxNumberOfLinearIterations(4);
//   system.SetAbsoluteLinearConvergenceTolerance(1.e-10);
 
  system.MGsolve();

  // ======= Print ========================
  std::vector < std::string > displ;
  displ.push_back("DX");
  displ.push_back("DY");
  if ( ml_mesh.GetDimension() == 3 ) displ.push_back("DZ");
  
  // set DX as moving variable
  ml_sol.GetWriter()->SetMovingMesh(displ);
  
  // print solutions
  std::vector < std::string > variablesToBePrinted;  variablesToBePrinted.push_back("All");
  ml_sol.GetWriter()->Write(files.GetOutputPath(),"biquadratic", variablesToBePrinted, lev);
 
 return ml_sol;

}

template < class system_type, class real_num, class real_num_mov = double >
void ComputeIntegral(const MultiLevelProblem& ml_prob) {
    ComputeIntegral< system_type, real_num, real_num_mov > (  ml_prob, 
                                                            & ml_prob.get_system< system_type >(0), 
                                                              ml_prob.get_system< system_type >(0).get_unknown_list_for_assembly());
}

template < class system_type, class real_num, class real_num_mov = double >
void ComputeIntegral(const MultiLevelProblem& ml_prob,
                     const system_type * mlPdeSys,
                     const std::vector< Unknown > & unknowns) {

  const unsigned 		level 		= mlPdeSys->GetLevelToAssemble();

  Mesh*          		msh    		= ml_prob._ml_msh->GetLevel(level);
  elem*          		el     		= msh->el;

  MultiLevelSolution*  		ml_sol  = ml_prob._ml_sol;
  Solution*    			sol      	= ml_prob._ml_sol->GetSolutionLevel(level); 


  const unsigned 	dim  = msh->GetDimension();
  const unsigned   dim2  = (3 * (dim - 1) + !(dim - 1));
  const unsigned  iproc  = msh->processor_id();
  // reserve memory for the local standard vectors
  const unsigned max_size_elem_dofs = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27 

  
  //----------- unknowns ------------------------------
   const unsigned int n_unknowns = mlPdeSys->GetSolPdeIndex().size();
   
  constexpr int sol_index_displ      = 0;     //known at compile time
  const     int sol_index_press      = sol_index_displ + dim;      //known at run time
  constexpr int state_pos_begin      = sol_index_displ;   //known at compile time
  const     int adj_pos_begin        = dim + 1;
  const     int sol_index_adj        = adj_pos_begin;      
  const     int sol_index_press_adj  = sol_index_adj + dim;     
  const     int ctrl_pos_begin       = 2 * (dim + 1);
  const     int sol_index_ctrl       = ctrl_pos_begin;     
  const     int theta_index = sol_index_ctrl + dim;      

  
  vector < std::string > Solname(n_unknowns);     
  vector < unsigned int > SolIndex(n_unknowns);  
  vector < unsigned int > SolPdeIndex(n_unknowns);
  vector < unsigned int > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
       Solname[ivar]    = unknowns[ivar]._name; 
      SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
//     SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(Solname[ivar].c_str());
//       assert(ivar == SolPdeIndex[ivar]);  //I would like this ivar order to coincide either with SolIndex or with SolPdeIndex, otherwise too many orders... it has to be  with SolPdeIndex, because SolIndex is related to the Solution object which can have more fields... This has to match the order of the unknowns[] argument
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);//
  }

   constexpr int solFEType_max = BIQUADR_FE;  //biquadratic, this is the highest-order FE 

   vector < unsigned int > Sol_n_el_dofs(n_unknowns);
  
  //----------- of dofs and at quadrature points ------------------------------
  vector < vector < double > > phi_dof_qp(n_unknowns);
  vector < vector < real_num_mov > > phi_x_dof_qp(n_unknowns);   //the derivatives depend on the Jacobian, which depends on the unknown, so we have to be adept here  //some of these should be real_num, some real_num_mov...
  vector < vector < real_num_mov > > phi_xx_dof_qp(n_unknowns);  //the derivatives depend on the Jacobian, which depends on the unknown, so we have to be adept here

  vector < vector < double > > phi_hat_dof_qp(n_unknowns);
  vector < vector < double > > phi_x_hat_dof_qp(n_unknowns);
  vector < vector < double > > phi_xx_hat_dof_qp(n_unknowns);

  for(int unk=0; unk < n_unknowns; unk++) {  
        phi_dof_qp[unk].reserve(max_size_elem_dofs);
    phi_hat_dof_qp[unk].reserve(max_size_elem_dofs);
      phi_x_dof_qp[unk].reserve(max_size_elem_dofs * dim);
  phi_x_hat_dof_qp[unk].reserve(max_size_elem_dofs * dim);
     phi_xx_dof_qp[unk].reserve(max_size_elem_dofs * dim2);
 phi_xx_hat_dof_qp[unk].reserve(max_size_elem_dofs * dim2);
   }
   
//----------boundary shape functions----------------------  
  vector < vector < double > > phi_hat_bd_dof_qp(n_unknowns);
  vector < vector < double > > phi_hat_bd_x_dof_qp(n_unknowns);


    for(int unk=0; unk < n_unknowns ; unk++) {  
        phi_hat_bd_dof_qp[unk].reserve(max_size_elem_dofs);
      phi_hat_bd_x_dof_qp[unk].reserve(max_size_elem_dofs * dim);
    }

  //----------- at dofs ------------------------------
  vector < vector < real_num > > SolVAR_eldofs(n_unknowns);  for(int k = 0; k < n_unknowns; k++) { SolVAR_eldofs[k].reserve(max_size_elem_dofs); }
 

  // geometry (at dofs) --------------------------------//
  vector< vector < real_num_mov > >   coords(dim);
  vector  < vector  < double > >  coords_hat(dim);
  vector  < vector  < double > >  coords_hat_bd(dim);
  unsigned coordsType = BIQUADR_FE; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  for(int i = 0; i < dim; i++) {   
           coords[i].reserve(max_size_elem_dofs); 
       coords_hat[i].reserve(max_size_elem_dofs); 
    coords_hat_bd[i].reserve(max_size_elem_dofs); 
 }
 
 //************** geometry phi **************************
  vector < vector < double > > phi_dof_qp_domain(dim);
  vector < vector < real_num_mov > > phi_x_dof_qp_domain(dim);
  vector < vector < real_num_mov > > phi_xx_dof_qp_domain(dim);
  vector < vector < double > > phi_dof_qp_domain_hat(dim);
  vector < vector < double > > phi_x_dof_qp_domain_hat(dim);
  vector < vector < double > > phi_xx_dof_qp_domain_hat(dim);
 
  for(int dom_coord = 0; dom_coord < dim; dom_coord++) {
        phi_dof_qp_domain[dom_coord].reserve(max_size_elem_dofs);
      phi_x_dof_qp_domain[dom_coord].reserve(max_size_elem_dofs * dim);
     phi_xx_dof_qp_domain[dom_coord].reserve(max_size_elem_dofs * dim2);
        phi_dof_qp_domain_hat[dom_coord].reserve(max_size_elem_dofs);
      phi_x_dof_qp_domain_hat[dom_coord].reserve(max_size_elem_dofs * dim);
     phi_xx_dof_qp_domain_hat[dom_coord].reserve(max_size_elem_dofs * dim2);
   }
 
   // geometry ------------------------------------------
  

  //------------ at quadrature points ---------------------
  real_num_mov weight_qp = 0.;//
    double weight_hat_qp = 0.;
    double weight_hat_bd_qp = 0.;
    vector < real_num > SolVAR_qp(n_unknowns);
    vector < real_num > SolVAR_hat_qp(n_unknowns);
    vector < vector < real_num > > gradSolVAR_qp(n_unknowns);
    vector < vector < real_num > > gradSolVAR_hat_qp(n_unknowns);
    for(int k = 0; k < n_unknowns; k++) { 
          gradSolVAR_qp[k].resize(dim);
      gradSolVAR_hat_qp[k].resize(dim);
    }
      


// Displ_desired##################################################################
  vector < double > phi_des_dof_qp;
  vector < real_num_mov > phi_des_x_dof_qp;   //the derivatives depend on the Jacobian, which depends on the unknown, so we have to be adept here  //some of these should be real_num, some real_num_mov...
  vector < real_num_mov > phi_des_xx_dof_qp;  

  vector < double > phi_des_hat_dof_qp;
  vector < double > phi_des_x_hat_dof_qp;
  vector < double > phi_des_xx_hat_dof_qp;

        phi_des_dof_qp.reserve(max_size_elem_dofs);
      phi_des_x_dof_qp.reserve(max_size_elem_dofs * dim);
     phi_des_xx_dof_qp.reserve(max_size_elem_dofs * dim2);
    phi_des_hat_dof_qp.reserve(max_size_elem_dofs);
  phi_des_x_hat_dof_qp.reserve(max_size_elem_dofs * dim);
 phi_des_xx_hat_dof_qp.reserve(max_size_elem_dofs * dim2);
   
   
  vector < real_num >  SolDisp_eldofs(dim,0.);
  vector < real_num > SolDisp_des_qp(dim, 0.);  
  
  vector < real_num > SolDisp_des_hat_qp(dim, 0.);  
// Displ_desired##################################################################




real_num  integral_target_alpha = 0.;
real_num	integral_beta   = 0.;
real_num	integral_gamma  = 0.;
real_num   integral_g_dot_n = 0.;

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {


   //all vars ###################################################################  
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
       Sol_n_el_dofs[k] = ndofs_unk;
       SolVAR_eldofs[k].resize(Sol_n_el_dofs[k]);
    for (unsigned i = 0; i < SolVAR_eldofs[k].size(); i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);
       SolVAR_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);
      }
    }
   //all vars ###################################################################
  
    
  // geometry *****************************
    short unsigned ielGeom = msh->GetElementType(iel);

      unsigned nDofsX = msh->GetElementDofNumber(iel, coordsType);    // number of coordinate element dofs
    
    for(int ivar=0; ivar<dim; ivar++) {
      coords[ivar].resize(nDofsX);
      coords_hat[ivar].resize(nDofsX);
    }
    
   for( unsigned i = 0; i < nDofsX; i++) {
      unsigned coordsDof  = msh->GetSolutionDof(i, iel, coordsType);    // global to global mapping between coordinates node and coordinate dof // via local to global solution node
      for(unsigned ivar = 0; ivar < dim; ivar++) {
          //Fixed coordinates (Reference frame)
	coords_hat[ivar][i] = (*msh->_topology->_Sol[ivar])(coordsDof);
      }
    }

  //***************************************  
    assert(nDofsX == Sol_n_el_dofs[sol_index_displ]);
    //Moving coordinates (Moving frame)
      for (unsigned idim = 0; idim < dim; idim++) {
        for (int j = 0; j < nDofsX; j++) {
          coords[idim][j] = coords_hat[idim][j] + SolVAR_eldofs[SolIndex[idim]][j];
        }
      }
  
     // elem average point 
    vector < double > elem_center(dim);   
    for (unsigned j = 0; j < dim; j++) {  elem_center[j] = 0.;  }
  for (unsigned j = 0; j < dim; j++) {  
      for (unsigned i = 0; i < nDofsX; i++) {
         elem_center[j] += coords_hat[j][i];
       }
    }
    
   for (unsigned j = 0; j < dim; j++) { elem_center[j] = elem_center[j]/nDofsX; }
  //*************************************** 
  
  //***** set target domain flag ********************************** 
   int target_flag = 0;
   target_flag = ElementTargetFlag < double >(elem_center); //domain is considered wrt moving

  // geometry end *****************************

 //************ set control flag *********************
    int control_el_flag = 0;
        control_el_flag = ControlDomainFlag(elem_center);
    std::vector< std::vector<int> > control_node_flag(dim);
	    for(unsigned idim=0; idim < dim; idim++) {
	          control_node_flag[idim].resize(Sol_n_el_dofs[ctrl_pos_begin]);
    std::fill(control_node_flag[idim].begin(), control_node_flag[idim].end(), 0);
	    }
 //*************************************************** 
  

  //DESIRED DISPL###################################################################  
    // displacement ************
//     for (unsigned i = 0; i < Sol_n_el_dofs[sol_index_displ]; i++) {
//       unsigned solDdesDof = msh->GetSolutionDof(i, iel, solDType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < SolDisp_eldofs.size() /*dim*/; k++) {
        SolDisp_eldofs[k]/*[i]*/ = TargetDisp[k] /*(*sol->_Sol[solDIndex[k]])(solDdesDof)*/;      // global extraction and local storage for the solution
      }
//     }
 //DESIRED DISPL###################################################################

 //========BoundaryLoop=====================================================================

  // Perform face loop over elements that contain some control face
  if (control_el_flag == 1) {
	  
      double tau = 0.;
      vector < double > normal(dim,0);
	  
      for(unsigned jface = 0; jface < msh->GetElementFaceNumber(iel); jface++) {
	  std::vector < double > xyz_bdc(3,0.);  //not being used, because the boundaries are identified by the face numbers
	  // look for boundary faces
	  if(el->GetFaceElementIndex(iel,jface) < 0) {
	      unsigned int face = -( msh->el->GetFaceElementIndex(iel,jface)+1);
	      if(  face == CTRL_FACE_IDX) { //control face
//=================================================== 
		unsigned nve_bd = msh->GetElementFaceDofNumber(iel,jface, SolFEType[ctrl_pos_begin] ); //AAAAAAAAAAAAAAAAA
		const unsigned felt_bd = msh->GetElementFaceType(iel, jface);    
		for(unsigned i=0; i < nve_bd; i++) {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i);
		    unsigned      iDof = msh->GetSolutionDof(i_vol, iel, coordsType);
		    for(unsigned idim=0; idim<dim; idim++) { coords_hat_bd[idim][i]=(*msh->_topology->_Sol[idim])(iDof); }
		}
		
//========= initialize gauss quantities on the boundary ============================================
		vector < real_num >                      SolVAR_bd_qp(n_unknowns);
		vector < vector < real_num > >       gradSolVAR_bd_qp(n_unknowns);
		for(int k=0; k<n_unknowns; k++) {    gradSolVAR_bd_qp[k].resize(dim);  }

//========= gauss_loop boundary===============================================================
		  for(unsigned ig_bd=0; ig_bd < ml_prob.GetQuadratureRule(felt_bd).GetGaussPointsNumber(); ig_bd++) {
		    for (unsigned  kdim = 0; kdim < dim; kdim++) {
		      ml_prob._ml_msh->_finiteElement[felt_bd][SolFEType[ctrl_pos_begin]]->JacobianSur(coords_hat_bd,ig_bd,weight_hat_bd_qp,phi_hat_bd_dof_qp[kdim + ctrl_pos_begin],phi_hat_bd_x_dof_qp[kdim + ctrl_pos_begin],normal);
            }
            
 set_phi_gradient_surface(AXIS_DIRECTION_CONTROL_SIDE, msh, ig_bd, SolFEType, felt_bd, nve_bd, coords_hat_bd, phi_hat_bd_x_dof_qp, dim, ctrl_pos_begin);
            
		  
//========== compute gauss quantities on the boundary ===============================================
		    for (unsigned  kdim = 0; kdim < dim; kdim++) {
			    unsigned int ctrl_index = kdim + ctrl_pos_begin;
										  SolVAR_bd_qp[ SolPdeIndex[ctrl_index] ] = 0.;
			  for(unsigned ivar2=0; ivar2<dim; ivar2++) {  gradSolVAR_bd_qp[ SolPdeIndex[ctrl_index] ][ivar2] = 0.; }
		    }//kdim
	  
		    for (unsigned  kdim = 0; kdim < dim; kdim++) {
			    unsigned int ctrl_index = kdim + ctrl_pos_begin;
			  for(int i_bd = 0; i_bd < nve_bd; i_bd++) {
		                  unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, jface, i_bd);
                        SolVAR_bd_qp[SolPdeIndex[ctrl_index]]   += phi_hat_bd_dof_qp  [ ctrl_index ][i_bd] * SolVAR_eldofs[ SolPdeIndex[ ctrl_index ] ][i_vol];
		      for(unsigned ivar2=0; ivar2<dim; ivar2++) {  
                      gradSolVAR_bd_qp[SolPdeIndex[ctrl_index]][ivar2] 	+= phi_hat_bd_x_dof_qp[ ctrl_index ][i_bd + ivar2 * nve_bd] * SolVAR_eldofs[ SolPdeIndex[ ctrl_index ] ][i_vol];      
                }
			  }
		    }//kdim
 //end unknowns eval at gauss points ********************************

 //========== compute gauss quantities on the boundary ================================================
      for (unsigned  k = 0; k < dim; k++) {
        integral_beta	 += SolVAR_bd_qp[SolPdeIndex[k + ctrl_pos_begin]] * SolVAR_bd_qp[SolPdeIndex[k + ctrl_pos_begin]] * weight_hat_bd_qp;
        integral_g_dot_n += SolVAR_bd_qp[SolPdeIndex[k + ctrl_pos_begin]] * normal[k] * weight_hat_bd_qp; 
      }
      for (unsigned  k = 0; k < dim; k++) {
            for (unsigned  j = 0; j < dim; j++) {	
		integral_gamma	  += gradSolVAR_bd_qp[SolPdeIndex[k + ctrl_pos_begin]][j] * gradSolVAR_bd_qp[SolPdeIndex[k + ctrl_pos_begin]][j] * weight_hat_bd_qp;
            }
      }
		  

                }  //end ig_bdry loop
	  
             }    //end if control face
	 }  //end if boundary faces
      }  // loop over element faces //jface   
  } //end if control element flag


     // *** Gauss point loop ***
      for (unsigned ig = 0; ig < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); ig++) {

	// *** get Jacobian and test function and test function derivatives ***
      for(int unk=0; unk < n_unknowns; unk++) {
	msh->_finiteElement[ielGeom][SolFEType[unk]]->Jacobian(coords,    ig,weight_qp,    phi_dof_qp[unk],    phi_x_dof_qp[unk],    phi_xx_dof_qp[unk]);
	msh->_finiteElement[ielGeom][SolFEType[unk]]->Jacobian(coords_hat,ig,weight_hat_qp,phi_hat_dof_qp[unk],phi_x_hat_dof_qp[unk],phi_xx_hat_dof_qp[unk]);
      }
      
   //=============== Integration at qp ========================================
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN...
    for (unsigned int d = 0; d < dim; d++) {
         msh->_finiteElement[ielGeom][coordsType]->Jacobian(coords,    ig, weight_qp,    phi_dof_qp_domain[d],     phi_x_dof_qp_domain[d],     phi_xx_dof_qp_domain[d]);
         msh->_finiteElement[ielGeom][coordsType]->Jacobian(coords_hat,ig, weight_hat_qp,phi_dof_qp_domain_hat[d], phi_x_dof_qp_domain_hat[d], phi_xx_dof_qp_domain_hat[d]);
      }
     
  //begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk <  n_unknowns; unk++) {
	  SolVAR_qp[SolPdeIndex[unk]] = 0.;
	  SolVAR_hat_qp[SolPdeIndex[unk]] = 0.;

	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR_qp[SolPdeIndex[unk]][ivar2] = 0.; 
	    gradSolVAR_hat_qp[SolPdeIndex[unk]][ivar2] = 0.; 
	  }
    }
	  
	for(unsigned unk = 0; unk <  n_unknowns; unk++) {
	  for(unsigned i = 0; i < Sol_n_el_dofs[unk]; i++) {
	    SolVAR_qp[SolPdeIndex[unk]]    += phi_dof_qp[ unk ][i] * SolVAR_eldofs[SolPdeIndex[unk]][i];
	    SolVAR_hat_qp[SolPdeIndex[unk]]    += phi_hat_dof_qp[ unk ][i] * SolVAR_eldofs[SolPdeIndex[unk]][i];
	    for(unsigned ivar2 = 0; ivar2 < dim; ivar2++) {
	      gradSolVAR_qp[SolPdeIndex[unk]][ivar2]     += phi_x_dof_qp[ unk ][i*dim+ivar2] * SolVAR_eldofs[SolPdeIndex[unk]][i]; 
	      gradSolVAR_hat_qp[SolPdeIndex[unk]][ivar2] += phi_x_hat_dof_qp[ unk ][i*dim+ivar2] * SolVAR_eldofs[SolPdeIndex[unk]][i]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************

	msh->_finiteElement[ielGeom][SolFEType[sol_index_displ]  /*SolDisp_eldofs*/]->Jacobian(coords, ig, weight_qp, phi_des_dof_qp, phi_des_x_dof_qp, phi_des_xx_dof_qp);
	msh->_finiteElement[ielGeom][SolFEType[sol_index_displ]  /*SolDisp_eldofs*/]->Jacobian(coords_hat, ig, weight_hat_qp, phi_des_hat_dof_qp, phi_des_x_hat_dof_qp, phi_des_xx_hat_dof_qp);

    for (unsigned  k = 0; k < dim; k++) {      SolDisp_des_qp[k]    = 0.;     }
    
      for (unsigned  k = 0; k < dim; k++) {
        for (unsigned i = 0; i < Sol_n_el_dofs[sol_index_displ]; i++) {
                SolDisp_des_qp[k] += SolDisp_eldofs[k]/*[i]*/ * phi_des_dof_qp[i];
            SolDisp_des_hat_qp[k] += SolDisp_eldofs[k]/*[i]*/ * phi_des_hat_dof_qp[i];
		}
      }
	
      
          
      for (unsigned  k = 0; k < dim; k++) {	 
          integral_target_alpha += target_flag * (SolVAR_qp[SolPdeIndex[k]]  - SolDisp_des_qp[k]) * 
                                                 (SolVAR_qp[SolPdeIndex[k]]  - SolDisp_des_qp[k]) * weight_qp; 
                }
   
  
      }// end gauss point loop
    } //end element loop  

       std::ostringstream filename_out; filename_out << ml_prob.GetFilesHandler()->GetOutputPath() << "/" << "Integral_computation"  << ".txt";

       std::ofstream intgr_fstream;
  if (paral::get_rank() == 0 ) {
      intgr_fstream.open(filename_out.str().c_str(),std::ios_base::app);
      intgr_fstream << " ***************************** Non Linear Iteration "<< mlPdeSys->GetNonlinearIt() << " *********************************** " <<  std::endl << std::endl;
      intgr_fstream << "The value of the target functional for " << "alpha " <<   std::setprecision(0) << std::scientific << alpha_val << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_target_alpha << std::endl;
      intgr_fstream << "The value of the L2 control for        " << "beta  " <<   std::setprecision(0) << std::scientific << beta_val  << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_beta         << std::endl;
      intgr_fstream << "The value of the H1 control for        " << "gamma " <<   std::setprecision(0) << std::scientific << gamma_val << " is " <<  std::setw(11) << std::setprecision(10) <<  integral_gamma        << std::endl;
      intgr_fstream << "The value of the integral of g.n "<<    integral_g_dot_n << std::endl;
      intgr_fstream << "The value of the theta is                             " <<    std::setw(11) << std::setprecision(10) <<  SolVAR_eldofs[SolPdeIndex[theta_index]][0] << std::endl;
      intgr_fstream << "The value of the total integral is " << std::setw(11) << std::setprecision(10) <<  integral_target_alpha * alpha_val*0.5  + integral_beta *beta_val*0.5 + integral_gamma *gamma_val*0.5 << std::endl;
      intgr_fstream <<  std::endl;
      intgr_fstream.close();  //you have to close to disassociate the file from the stream
}  

    return; 
	  
  
}

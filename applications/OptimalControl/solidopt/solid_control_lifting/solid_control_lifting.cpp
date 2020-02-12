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

#include   "../solidopt_params.hpp"

using namespace femus;



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
             
             else if(!strcmp(name,"DX_CTRL")) {
                 value = 0.;
             }
             else if(!strcmp(name,"DY_CTRL")) {
                 value = 0.;
             }
             else if(!strcmp(name,"DZ_CTRL")) {
                 value = 0.;
             }
             else if(!strcmp(name,"P_CTRL")) {
                 value = 0.;
             }
          
    
      return value;   
}  
  

bool Solution_set_boundary_conditions(const MultiLevelProblem * ml_prob, const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  //1: bottom (y=y_min) //2: right (x=x_max)  //3: top (y=y_max)  //4: left (x=x_min)  (in 2D)
  //1: bottom (y=y_min) //2: right (x=x_max)  //3: top (y=y_max)  //4: left (x=x_min) //5: (z=z_min)  //6:  (z=z_max) (in 3D, I guess...)
  
  bool dirichlet  = true; 
  value = 0.;
  
        if (facename == 2) {
       if (!strcmp(SolName, "DX_CTRL"))   { dirichlet = false; }
  else if (!strcmp(SolName, "DY_CTRL"))   { dirichlet = false; } 
  else if (!strcmp(SolName, "DZ_CTRL"))   { dirichlet = false; } 
  	
      }

  
// //       if (facename == 1) {
// //        if (!strcmp(SolName, "DX"))        { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY"))        { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ"))        { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DX_ADJ"))    { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY_ADJ"))    { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ_ADJ"))    { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DX_CTRL"))   { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY_CTRL"))   { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ_CTRL"))   { dirichlet = true; value = 0.; } 
// //   	
// //       }
// // 
// //       if (facename == 2) {
// //        if (!strcmp(SolName, "DX"))        { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY"))        { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ"))        { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DX_ADJ"))    { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY_ADJ"))    { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ_ADJ"))    { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DX_CTRL"))   { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY_CTRL"))   { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ_CTRL"))   { dirichlet = true; value = 0.; } 
// //   	
// //       }
// // 
// //       if (facename == 3) {
// //        if (!strcmp(SolName, "DX"))        { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY"))        { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ"))        { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DX_ADJ"))    { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY_ADJ"))    { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ_ADJ"))    { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DX_CTRL"))   { dirichlet = false; value = 0.; }
// //   else if (!strcmp(SolName, "DY_CTRL"))   { dirichlet = false; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ_CTRL"))   { dirichlet = false; value = 0.; } 
// //   	
// //       }
// // 
// //       if (facename == 4) {
// //        if (!strcmp(SolName, "DX"))        { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY"))        { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ"))        { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DX_ADJ"))    { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY_ADJ"))    { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ_ADJ"))    { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DX_CTRL"))   { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY_CTRL"))   { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ_CTRL"))   { dirichlet = true; value = 0.; } 
// //   	
// //       }
// //       
// //       if (facename == 5) {
// //        if (!strcmp(SolName, "DX"))        { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY"))        { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ"))        { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DX_ADJ"))    { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DY_ADJ"))    { dirichlet = true; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ_ADJ"))    { dirichlet = true; value = 0.; }
// //   else if (!strcmp(SolName, "DX_CTRL"))   { dirichlet = false; value = 0.; }
// //   else if (!strcmp(SolName, "DY_CTRL"))   { dirichlet = false; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ_CTRL"))   { dirichlet = true; value = 0.; } 
// //      }
// //       
// //       if (facename == 6) {
// //        if (!strcmp(SolName, "DX"))        { dirichlet = true; value = 0. ; }
// //   else if (!strcmp(SolName, "DY"))        { dirichlet = true; value = 0. ; } 
// //   else if (!strcmp(SolName, "DZ"))        { dirichlet = true; value = 0. ; }
// //   else if (!strcmp(SolName, "DX_ADJ"))    { dirichlet = true; value = 0. ; }
// //   else if (!strcmp(SolName, "DY_ADJ"))    { dirichlet = true; value = 0. ; } 
// //   else if (!strcmp(SolName, "DZ_ADJ"))    { dirichlet = true; value = 0. ; }
// //   else if (!strcmp(SolName, "DX_CTRL"))   { dirichlet = false; value = 0.; }
// //   else if (!strcmp(SolName, "DY_CTRL"))   { dirichlet = false; value = 0.; } 
// //   else if (!strcmp(SolName, "DZ_CTRL"))   { dirichlet = true; value = 0.; } 
// //       }
      
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
                        feFamily.push_back(LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/);
 
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
                        feOrder.push_back(FIRST);
 

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
                       unknowns[ctrl_pos_begin + 0]._name      = "DX_CTRL";
                       unknowns[ctrl_pos_begin + 1]._name      = "DY_CTRL";
  if (dimension == 3)  unknowns[ctrl_pos_begin + 2]._name      = "DZ_CTRL";
               unknowns[ctrl_pos_begin + dimension]._name      = "P_CTRL";

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
 


int main(int argc, char** args) {

  // ======= Init ==========================
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);

  // ======= Files ==================
  Files files;
        files.CheckIODirectories();
        files.RedirectCout();

    std::string mesh_folder_file = "input/";
    std::string input_file = "cyl.med";
  std::ostringstream mystream; mystream << "./" << /*DEFAULT_INPUTDIR*/ mesh_folder_file << input_file;
  const std::string infile = mystream.str();

  // ======= Quad Rule ========================
  std::string quad_rule_order("fifth");

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
//   ml_mesh.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,quad_rule_order.c_str());
   ml_mesh.ReadCoarseMesh(infile.c_str(),quad_rule_order.c_str(), 1./*Lref*/);

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
  const     int sol_index_press_ctrl = sol_index_ctrl + dim;      

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
 

  // geometry (at dofs) --------------------------------
  vector< vector < real_num_mov > >   coords(dim);
  vector  < vector  < double > >  coords_hat(dim);
  unsigned coordsType = BIQUADR_FE; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  for(int i = 0; i < dim; i++) {   
           coords[i].reserve(max_size_elem_dofs); 
       coords_hat[i].reserve(max_size_elem_dofs); 
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
    const double mus   = mu_lame/rhos;

    const bool incompressible   = (0.5  ==  solid_in.get_poisson_coeff()) ? 1 : 0;
    const bool penalty = solid_in.get_if_penalty();

    const double young_modulus = solid_in.get_young_module();
    const double Lref = solid_in._parameter->Get_reference_length();
    const double nondim_number = Lref/young_modulus;
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
  for (unsigned  k = 0; k < n_unknowns; k++) {
    unsigned ndofs_unk = msh->GetElementDofNumber(iel, SolFEType[k]);
       Sol_n_el_dofs[k] = ndofs_unk;
       SolVAR_eldofs[k].resize(Sol_n_el_dofs[k]);
          L2G_dofmap[k].resize(Sol_n_el_dofs[k]); 
    for (unsigned i = 0; i < SolVAR_eldofs[k].size(); i++) {
       unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);
       SolVAR_eldofs[k][i] = (*sol->_Sol[SolIndex[k]])(solDof);
          L2G_dofmap[k][i] = pdeSys->GetSystemDof(SolIndex[k], SolPdeIndex[k], i, iel);    // global to global mapping between solution node and pdeSys dof
      }
    }
    
    L2G_dofmap_AllVars.resize(0);
      for (unsigned  k = 0; k < n_unknowns; k++)     L2G_dofmap_AllVars.insert(L2G_dofmap_AllVars.end(),L2G_dofmap[k].begin(),L2G_dofmap[k].end());

    unsigned sum_Sol_n_el_dofs = 0;
    for (unsigned  k = 0; k < n_unknowns; k++) { sum_Sol_n_el_dofs += Sol_n_el_dofs[k]; }
    Jac.resize(sum_Sol_n_el_dofs * sum_Sol_n_el_dofs);  std::fill(Jac.begin(), Jac.end(), 0.);
    Res.resize(sum_Sol_n_el_dofs);                      std::fill(Res.begin(), Res.end(), 0.);
 
 //=================state_only - for Jac retrieval via AD
    unsigned state_sum_Sol_n_el_dofs = 0;
    for (unsigned  k = 0; k < n_state_vars; k++) { state_sum_Sol_n_el_dofs += Sol_n_el_dofs[k]; }
    Jac_state_true.resize(state_sum_Sol_n_el_dofs * state_sum_Sol_n_el_dofs);    std::fill(Jac_state_true.begin(), Jac_state_true.end(), 0.);
    Jac_state_false.resize(state_sum_Sol_n_el_dofs * state_sum_Sol_n_el_dofs);    std::fill(Jac_state_false.begin(), Jac_state_false.end(), 0.);
    Res_state.resize(state_sum_Sol_n_el_dofs);                              std::fill(Res_state.begin(), Res_state.end(), 0.);
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
   target_flag = ElementTargetFlag < double > (elem_center);

  // geometry end *****************************

      
     assemble_jac.prepare_before_integration_loop(stack);
     assemble_jac_state.prepare_before_integration_loop(stack_state); //=================state_only - for Jac retrieval via AD
     
  
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
	  SolVAR_qp[unk] = 0.;
	  SolVAR_hat_qp[unk] = 0.;

	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR_qp[unk][ivar2] = 0.; 
	    gradSolVAR_hat_qp[unk][ivar2] = 0.; 
	  }
	  
	  for(unsigned i = 0; i < Sol_n_el_dofs[unk]; i++) {
	    SolVAR_qp[unk]    += phi_dof_qp[ unk ][i] * SolVAR_eldofs[unk][i];
	    SolVAR_hat_qp[unk]    += phi_hat_dof_qp[ unk ][i] * SolVAR_eldofs[unk][i];
	    for(unsigned ivar2 = 0; ivar2 < dim; ivar2++) {
	      gradSolVAR_qp[unk][ivar2]     += phi_x_dof_qp[ unk ][i*dim+ivar2] * SolVAR_eldofs[unk][i]; 
	      gradSolVAR_hat_qp[unk][ivar2] += phi_x_hat_dof_qp[ unk ][i*dim+ivar2] * SolVAR_eldofs[unk][i]; 
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
    - alpha_val * target_flag * (SolVAR_hat_qp[SolPdeIndex[idim]] + SolVAR_hat_qp[SolPdeIndex[idim + ctrl_pos_begin]] - TargetDisp[idim]) 
     * phi_dof_qp[ idim + adj_pos_begin ][i] * weight_qp;
               }
          }
//----------from displ_state and displ_ctrl
 //ADJOINT=====================================
 
        
  //CONTROL=====================================
 //----------from displ_state and displ_ctrl
        for (unsigned i = 0; i < Sol_n_el_dofs[sol_index_ctrl]; i++) {
               for (int idim = 0; idim < dim; idim++) {
                    for (unsigned jdim = 0; jdim < dim; jdim++) {
   Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[idim + ctrl_pos_begin], i)] += 
                    gamma_val * ( gradSolVAR_hat_qp[idim + ctrl_pos_begin][jdim] + gradSolVAR_hat_qp[jdim + ctrl_pos_begin][idim] )
                    * phi_x_dof_qp[idim + ctrl_pos_begin][i * dim + jdim] * weight_qp;                   
                    }
   Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[idim + ctrl_pos_begin], i)] += 
            (   ( alpha_val * target_flag * (SolVAR_hat_qp[SolPdeIndex[idim]] + SolVAR_hat_qp[SolPdeIndex[idim + ctrl_pos_begin]] - TargetDisp[idim])  
                + beta_val * SolVAR_hat_qp[SolPdeIndex[idim + ctrl_pos_begin]] )* phi_dof_qp[ idim + ctrl_pos_begin ][i] 
            ) * weight_qp;
               }
          }
 //----------from displ_state and displ_ctrl
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

 
  //STATE=====================================
//----------from displ_ctrl only

      for(unsigned i_block = state_pos_begin; i_block < state_pos_begin + dim; i_block++) {
        for(unsigned i_dof = 0; i_dof < Sol_n_el_dofs[i_block]; i_dof++) {
            for(unsigned j_block = state_pos_begin; j_block < state_pos_begin + dim; j_block++) {
               for(unsigned j_dof = 0; j_dof < Sol_n_el_dofs[j_block]; j_dof++) {
    Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, i_block, i_dof)] += 
        Jac_state_true[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, state_sum_Sol_n_el_dofs,i_block , j_block, i_dof, j_dof )] 
            * SolVAR_eldofs[j_block + ctrl_pos_begin][j_dof];
                }
           }
        }
   }
//----------from displ_ctrl only
  //STATE=====================================


    
  //ADJOINT=====================================
//----------from adj displ
      for(unsigned i_block = adj_pos_begin; i_block < adj_pos_begin + dim; i_block++) {
        for(unsigned i_dof = 0; i_dof < Sol_n_el_dofs[i_block]; i_dof++) {
            for(unsigned j_block = 0; j_block < dim; j_block++) {
               for(unsigned j_dof = 0; j_dof < Sol_n_el_dofs[j_block]; j_dof++) {
    Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, i_block, i_dof)] += 
        Jac_state_false[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, state_sum_Sol_n_el_dofs,i_block - adj_pos_begin , j_block, i_dof, j_dof )] 
            * SolVAR_eldofs[j_block + adj_pos_begin][j_dof];
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
            * SolVAR_eldofs[j_block + adj_pos_begin][j_dof];
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
            * SolVAR_eldofs[j_block + adj_pos_begin][j_dof];
                }
           }
        }
   }
// ----------from div of displ_adj
  //ADJOINT=====================================

   //CONTROL=====================================

//----------from displ__adj only
      for(unsigned i_block = ctrl_pos_begin; i_block < ctrl_pos_begin + dim ; i_block++) {
        for(unsigned i_dof = 0; i_dof < Sol_n_el_dofs[i_block]; i_dof++) {
            for(unsigned j_block = 0; j_block < dim; j_block++) {
               for(unsigned j_dof = 0; j_dof < Sol_n_el_dofs[j_block]; j_dof++) {
    Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, i_block, i_dof)] += 
        - Jac_state_false[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, state_sum_Sol_n_el_dofs,i_block - ctrl_pos_begin , j_block, i_dof, j_dof )] 
          * SolVAR_eldofs[j_block + adj_pos_begin][j_dof];
                }
           }
        }
   }
//----------from displ__adj only

//----------from p_ctrl only
   for(unsigned i_block = ctrl_pos_begin; i_block < ctrl_pos_begin + dim ; i_block++) {
        for(unsigned i_dof = 0; i_dof < Sol_n_el_dofs[i_block]; i_dof++) {
            for(unsigned j_block = dim; j_block < n_state_vars; j_block++) {
               for(unsigned j_dof = 0; j_dof < Sol_n_el_dofs[j_block]; j_dof++) {
    Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, i_block, i_dof)] += 
        Jac_state_true[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, state_sum_Sol_n_el_dofs,i_block - ctrl_pos_begin , j_block, i_dof, j_dof )] 
            * SolVAR_eldofs[j_block + ctrl_pos_begin][j_dof];
                }
           }
        }
   }
//----------from p_ctrl only

//----------from div of displ_ctrl
   for(unsigned i_block = sol_index_press_ctrl; i_block < n_unknowns ; i_block++) {
        for(unsigned i_dof = 0; i_dof < Sol_n_el_dofs[i_block]; i_dof++) {
            for(unsigned j_block = 0; j_block < dim; j_block++) {
               for(unsigned j_dof = 0; j_dof < Sol_n_el_dofs[j_block]; j_dof++) {
    Res[assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, i_block, i_dof)] += 
        Jac_state_true[assemble_jacobian<double,double>::jac_row_col_index(Sol_n_el_dofs, state_sum_Sol_n_el_dofs,i_block - ctrl_pos_begin , j_block, i_dof, j_dof )] 
            * SolVAR_eldofs[j_block + ctrl_pos_begin][j_dof];
                }
           }
        }
   }
//----------from div of displ_ctrl
  //CONTROL=====================================
  
  
      
     if (assembleMatrix) assemble_jac.compute_jacobian_outside_integration_loop(stack, SolVAR_eldofs, Res, Jac, L2G_dofmap_AllVars, RES, JAC);



       
//     assemble_jacobian<real_num,real_num_mov>::print_element_residual(iel, Res, Sol_n_el_dofs, 9, 3);
//     assemble_jacobian<real_num,real_num_mov>::print_element_jacobian(iel, Jac, Sol_n_el_dofs, 9, 3);
    
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

//   system.SetMaxNumberOfNonLinearIterations(2);
//   system.SetNonLinearConvergenceTolerance(1.e-2);
//   system.SetDebugLinear(true);
//   system.SetMaxNumberOfLinearIterations(4);
//   system.SetAbsoluteLinearConvergenceTolerance(1.e-10);
 
  system.MGsolve();

  // ======= Print ========================
  std::vector < std::string > displ_tot;
  displ_tot.push_back("DX_TOT");
  displ_tot.push_back("DY_TOT");
  if ( ml_mesh.GetDimension() == 3 ) displ_tot.push_back("DZ_TOT");
  
  std::vector < std::string > displ_hom;
  displ_hom.push_back("DX");
  displ_hom.push_back("DY");
  if ( ml_mesh.GetDimension() == 3 ) displ_hom.push_back("DZ");
  
  std::vector < std::string > displ_ctrl;
  displ_ctrl.push_back("DX_CTRL");
  displ_ctrl.push_back("DY_CTRL");
  if ( ml_mesh.GetDimension() == 3 ) displ_ctrl.push_back("DZ_CTRL");
  
  const unsigned dim = ml_mesh.GetDimension();
  
  for (unsigned d = 0; d < dim; d++)   assert(  ml_sol.GetSolutionType( displ_hom[d].c_str() ) == ml_sol.GetSolutionType( displ_ctrl[d].c_str() ) );
     
  for (unsigned d = 0; d < dim; d++)   {
     ml_sol.AddSolution(displ_tot[d].c_str(),  ml_sol.GetSolutionFamily( displ_hom[d] ), ml_sol.GetSolutionOrder( displ_hom[d] ), ml_sol.GetSolutionTimeOrder( displ_hom[d] ), false);
     ml_sol.Initialize(displ_tot[d].c_str()); //minimum that is needed to add and print a Solution without errors
  }
  
  // DX_TOT = DX + DX_CTRL
  for (unsigned d = 0; d < dim; d++)   {
      /**(*/ ml_sol.GetSolutionLevel(lev-1)->GetSolutionName(displ_tot[d].c_str()) /*)*/  = /**(*/ ml_sol.GetSolutionLevel(lev-1)->GetSolutionName(displ_hom[d].c_str()) /*)*/; 
      /**(*/ ml_sol.GetSolutionLevel(lev-1)->GetSolutionName(displ_tot[d].c_str()) /*)*/ += /**(*/ ml_sol.GetSolutionLevel(lev-1)->GetSolutionName(displ_ctrl[d].c_str()) /*)*/; 
  }
  

  // set DX_TOT as moving variable
  ml_sol.GetWriter()->SetMovingMesh(displ_tot);
  
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
  const     int sol_index_press_ctrl = sol_index_ctrl + dim;      

  
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
   
  //----------- at dofs ------------------------------
  vector < vector < real_num > > SolVAR_eldofs(n_unknowns);  for(int k = 0; k < n_unknowns; k++) { SolVAR_eldofs[k].reserve(max_size_elem_dofs); }
 

  // geometry (at dofs) --------------------------------//
  vector< vector < real_num_mov > >   coords(dim);
  vector  < vector  < double > >  coords_hat(dim);
  unsigned coordsType = BIQUADR_FE; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  for(int i = 0; i < dim; i++) {   
           coords[i].reserve(max_size_elem_dofs); 
       coords_hat[i].reserve(max_size_elem_dofs); 
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


  //DESIRED DISPL###################################################################  
    // displacement ************
//     for (unsigned i = 0; i < Sol_n_el_dofs[sol_index_displ]; i++) {
//       unsigned solDdesDof = msh->GetSolutionDof(i, iel, solDType);    // global to global mapping between solution node and solution dof

      for (unsigned  k = 0; k < SolDisp_eldofs.size() /*dim*/; k++) {
        SolDisp_eldofs[k]/*[i]*/ = TargetDisp[k] /*(*sol->_Sol[solDIndex[k]])(solDdesDof)*/;      // global extraction and local storage for the solution
      }
//     }
 //DESIRED DISPL###################################################################

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
	  SolVAR_qp[unk] = 0.;
	  SolVAR_hat_qp[unk] = 0.;

	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR_qp[unk][ivar2] = 0.; 
	    gradSolVAR_hat_qp[unk][ivar2] = 0.; 
	  }
	  
	  for(unsigned i = 0; i < Sol_n_el_dofs[unk]; i++) {
	    SolVAR_qp[unk]    += phi_dof_qp[ unk ][i] * SolVAR_eldofs[unk][i];
	    SolVAR_hat_qp[unk]    += phi_hat_dof_qp[ unk ][i] * SolVAR_eldofs[unk][i];
	    for(unsigned ivar2 = 0; ivar2 < dim; ivar2++) {
	      gradSolVAR_qp[unk][ivar2]     += phi_x_dof_qp[ unk ][i*dim+ivar2] * SolVAR_eldofs[unk][i]; 
	      gradSolVAR_hat_qp[unk][ivar2] += phi_x_hat_dof_qp[ unk ][i*dim+ivar2] * SolVAR_eldofs[unk][i]; 
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
          integral_target_alpha += target_flag * (SolVAR_qp[SolIndex[k]] + SolVAR_qp[SolIndex[k + ctrl_pos_begin]] - SolDisp_des_qp[k]) * 
                                                 (SolVAR_qp[SolIndex[k]] + SolVAR_qp[SolIndex[k + ctrl_pos_begin]] - SolDisp_des_qp[k]) * weight_qp; 
          integral_beta	 += SolVAR_qp[SolIndex[k + ctrl_pos_begin]] *SolVAR_qp[SolIndex[k + ctrl_pos_begin]] * weight_qp;

                }
      for (unsigned  k = 0; k < dim; k++) {
        for (unsigned  j = 0; j < dim; j++) {	
            integral_gamma	  +=  gradSolVAR_qp[SolIndex[k + ctrl_pos_begin]][j] * gradSolVAR_qp[SolIndex[k + ctrl_pos_begin]][j] * weight_qp;
        }
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
      intgr_fstream << "The value of the total integral is " << std::setw(11) << std::setprecision(10) <<  integral_target_alpha * alpha_val*0.5  + integral_beta *beta_val*0.5 + integral_gamma *gamma_val*0.5 << std::endl;
      intgr_fstream <<  std::endl;
      intgr_fstream.close();  //you have to close to disassociate the file from the stream
}  

    return; 
	  
  
}

#include "FemusInit.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
#include "LinearImplicitSystem.hpp"
#include "Parameter.hpp"
#include "Fluid.hpp"
#include "Solid.hpp"
#include "Files.hpp"
#include "FE_convergence.hpp"
#include "Assemble_jacobian.hpp"


#define MODEL "Linear_elastic"
// #define MODEL "Mooney-Rivlin" 
// #define MODEL "Neo-Hookean"

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
           
    
      return value;   
}  
  
  
  

bool Solution_set_boundary_conditions(const std::vector < double >& x, const char SolName[], double& value, const int facename, const double time) {
  //1: bottom (y=y_min) //2: right (x=x_max)  //3: top (y=y_max)  //4: left (x=)  (in 2D)
  //1: bottom (y=y_min) //2: right (x=x_max)  //3: top (y=y_max)  //4: left (x=) //5: (z=z_min)  //6:  (z=z_max) (in 3D, I guess...)
  
  bool dirichlet; 
  
      if (facename == 1) {
       if (!strcmp(SolName, "DX"))    { dirichlet = true; value = 0.; }
  else if (!strcmp(SolName, "DY"))    { dirichlet = true; value = 0.; } 
  else if (!strcmp(SolName, "DZ"))    { dirichlet = true; value = 0.; } 
  	
      }

      if (facename == 2) {
       if (!strcmp(SolName, "DX"))    { dirichlet = false/*true*/; value = 0.; }
  else if (!strcmp(SolName, "DY"))    { dirichlet = false/*true*/; value = 0.; } 
  else if (!strcmp(SolName, "DZ"))    { dirichlet = true; value = 0.; } 
  	
      }

      if (facename == 3) {
       if (!strcmp(SolName, "DX"))    { dirichlet = true; value = 0.; }
  else if (!strcmp(SolName, "DY"))    { dirichlet = true; value = 0.; } 
  else if (!strcmp(SolName, "DZ"))    { dirichlet = true; value = 0.; } 
  	
      }

      if (facename == 4) {
       if (!strcmp(SolName, "DX"))    { dirichlet = false/*true*/; value = 0.; }
  else if (!strcmp(SolName, "DY"))    { dirichlet = false/*true*/; value = 0.; } 
  else if (!strcmp(SolName, "DZ"))    { dirichlet = true; value = 0.; } 
  	
      }
      
      if (facename == 5) {
       if (!strcmp(SolName, "DX"))    { dirichlet = false; value = 0.; }
  else if (!strcmp(SolName, "DY"))    { dirichlet = false; value = 0.; } 
  else if (!strcmp(SolName, "DZ"))    { dirichlet = true; value = 0.; } 
      }
      
      if (facename == 6) {
       if (!strcmp(SolName, "DX"))    { dirichlet = false; value = 0.; }
  else if (!strcmp(SolName, "DY"))    { dirichlet = false; value = 0.; } 
  else if (!strcmp(SolName, "DZ"))    { dirichlet = true; value = 0.; } 
      }
      
  return dirichlet;
}



template < class system_type, class real_num, class real_num_mov >
void AssembleSolidMech(MultiLevelProblem& ml_prob);


template < class system_type, class real_num, class real_num_mov >
void AssembleSolidMech(MultiLevelProblem& ml_prob,
                       system_type * mlPdeSys,
                       const std::vector< Math::Unknown > &  unknowns);





 //Unknown definition  ==================
 const std::vector< Math::Unknown >  provide_list_of_unknowns(const unsigned int dimension) {
     
     
  std::vector< FEFamily > feFamily;
  std::vector< FEOrder >   feOrder;

                        feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE);
  if (dimension == 3)   feFamily.push_back(LAGRANGE);
                        feFamily.push_back(LAGRANGE/*DISCONTINOUS_POLYNOMIAL*/);
  
                        feOrder.push_back(SECOND);
                        feOrder.push_back(SECOND);
  if (dimension == 3)   feOrder.push_back(SECOND);
                        feOrder.push_back(FIRST);
  
  assert( feFamily.size() == feOrder.size() );
 
 std::vector< Math::Unknown >  unknowns(feFamily.size());

                        unknowns[0]._name      = "DX";
                        unknowns[1]._name      = "DY";
  if (dimension == 3)   unknowns[2]._name      = "DZ";
                unknowns[dimension]._name      = "P";
 
     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
              unknowns[u]._fe_family = feFamily[u];
              unknowns[u]._fe_order  = feOrder[u];
              
     }
 
 
   return unknowns;
     
}




template < class real_num > 
class My_main_single_level : public Main_single_level {
    
public:
    
const MultiLevelSolution  run_on_single_level(const Files & files, 
                                                   const std::vector< Math::Unknown > & unknowns,  
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

  // ======= Quad Rule ========================
  std::string fe_quad_rule("seventh");

  // ======= Mesh ==================
  const ElemType geom_elem_type = QUAD9;
  const std::vector< unsigned int > nsub = {4, 4, 0};
  const std::vector< double >      xyz_min = {0., 0., 0.};
  const std::vector< double >      xyz_max = {1., 1., 0.};

  MultiLevelMesh ml_mesh;
  ml_mesh.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,fe_quad_rule.c_str());

  // ======= Unknowns ========================
  const unsigned int dimension = ml_mesh.GetDimension();  
  std::vector< Math::Unknown > unknowns = provide_list_of_unknowns(dimension);
  
  // ======= Normal run ========================   //if you don't want the convergence study
  My_main_single_level< /*adept::a*/double > my_main;
//   const unsigned int n_levels = 1;
//   my_main.run_on_single_level(files, unknowns, ml_mesh, n_levels); 
 
  
  
  // ======= Convergence study ========================
    
   //set coarse storage mesh (should write the copy constructor or "=" operator to copy the previous mesh) ==================
   MultiLevelMesh ml_mesh_all_levels;
   ml_mesh_all_levels.GenerateCoarseBoxMesh(nsub[0],nsub[1],nsub[2],xyz_min[0],xyz_max[0],xyz_min[1],xyz_max[1],xyz_min[2],xyz_max[2],geom_elem_type,fe_quad_rule.c_str());
   //   ml_mesh_all_levels.ReadCoarseMesh(infile.c_str(),fe_quad_rule.c_str(),1.);
 
   // convergence choices ================  
   unsigned int max_number_of_meshes;               // set total number of levels ================  

   if (nsub[2] == 0)   max_number_of_meshes = 6;
   else                max_number_of_meshes = 4;
  
//    My_exact_solution<> exact_sol;                //provide exact solution, if available ==============
   const unsigned conv_order_flag = 0;              //Choose how to compute the convergence order ========= //0: incremental 1: absolute (with analytical sol)  2: absolute (with projection of finest sol)...
   const unsigned norm_flag = 1;                    //Choose what norms to compute (//0 = only L2: //1 = L2 + H1) ==============


   // object ================  
    FE_convergence<>  fe_convergence;
    
    fe_convergence.convergence_study(files, unknowns, Solution_set_boundary_conditions, ml_mesh, ml_mesh_all_levels, max_number_of_meshes, norm_flag, conv_order_flag, my_main);
  
  return 0;
}



template < class system_type, class real_num, class real_num_mov = double >
void AssembleSolidMech(MultiLevelProblem& ml_prob) {
    
  AssembleSolidMech< system_type, real_num, real_num_mov > (ml_prob, & ml_prob.get_system< system_type >(0), ml_prob.get_unknown_list_for_assembly());

}



template < class system_type, class real_num, class real_num_mov = double >
void AssembleSolidMech(MultiLevelProblem& ml_prob,
                       system_type * mlPdeSys,
                       const std::vector< Math::Unknown > &  unknowns) {
    
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


  const unsigned 	 dim = msh->GetDimension();
  const unsigned    dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned   nel   = msh->GetNumberOfElements();
  const unsigned igrid   = msh->GetLevel();
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
   
  constexpr int sol_index_displ = 0;     //known at compile time
  const     int sol_index_press = dim;      //known at run time
  constexpr int state_pos_begin = sol_index_displ;   //known at compile time


  vector < std::string > Solname(n_unknowns);     for(unsigned ivar=0; ivar < n_unknowns; ivar++) { Solname[ivar] = unknowns[ivar]._name; }
//   Solname              [state_pos_begin + 0] =                "DX";
//   Solname              [state_pos_begin + 1] =                "DY";
//   if (dim == 3) Solname[state_pos_begin + 2] =                "DZ";
//   Solname              [state_pos_begin + sol_index_press] = "P";
  
  vector < unsigned int > SolIndex(n_unknowns);  
  vector < unsigned int > SolPdeIndex(n_unknowns);
  vector < unsigned int > SolFEType(n_unknowns);  


  for(unsigned ivar=0; ivar < n_unknowns; ivar++) {
    SolIndex[ivar]	= ml_sol->GetIndex        (Solname[ivar].c_str());
    SolPdeIndex[ivar]	= mlPdeSys->GetSolPdeIndex(Solname[ivar].c_str());
    SolFEType[ivar]	= ml_sol->GetSolutionType(SolIndex[ivar]);
  }

   vector < unsigned int > Sol_n_el_dofs(n_unknowns);
  
  //----------- of dofs and at quadrature points ------------------------------
  vector < vector < double > > phi_dof_qp(NFE_FAMS);
  vector < vector < double > > phi_hat_dof_qp(NFE_FAMS);
  vector < vector < real_num_mov > > phi_x_dof_qp(NFE_FAMS);   //the derivatives depend on the Jacobian, which depends on the unknown, so we have to be adept here  //some of these should be real_num, some real_num_mov...
  vector < vector < real_num_mov > > phi_xx_dof_qp(NFE_FAMS);  //the derivatives depend on the Jacobian, which depends on the unknown, so we have to be adept here
  vector < vector < double > > phi_x_hat_dof_qp(NFE_FAMS);
  vector < vector < double > > phi_xx_hat_dof_qp(NFE_FAMS);

  for(int fe=0; fe < NFE_FAMS; fe++) {  
        phi_dof_qp[fe].reserve(max_size_elem_dofs);
    phi_hat_dof_qp[fe].reserve(max_size_elem_dofs);
      phi_x_dof_qp[fe].reserve(max_size_elem_dofs * dim);
  phi_x_hat_dof_qp[fe].reserve(max_size_elem_dofs * dim);
     phi_xx_dof_qp[fe].reserve(max_size_elem_dofs * dim2);
 phi_xx_hat_dof_qp[fe].reserve(max_size_elem_dofs * dim2);
   }
   
  //----------- at dofs ------------------------------
  vector < double >   Jac;   Jac.reserve( n_unknowns * max_size_elem_dofs * n_unknowns * max_size_elem_dofs);
  vector < real_num > Res;   Res.reserve( n_unknowns * max_size_elem_dofs);
           vector < int >       L2G_dofmap_AllVars;   L2G_dofmap_AllVars.reserve( n_unknowns *max_size_elem_dofs);
  vector < vector < int > >         L2G_dofmap(n_unknowns);  for(int i = 0; i < n_unknowns; i++) {    L2G_dofmap[i].reserve(max_size_elem_dofs); }
  vector < vector < real_num > > SolVAR_eldofs(n_unknowns);  for(int k = 0; k < n_unknowns; k++) { SolVAR_eldofs[k].reserve(max_size_elem_dofs); }
  

  // geometry (at dofs) --------------------------------
  vector< vector < real_num_mov > >   coords(dim);
  vector  < vector  < double > >  coords_hat(dim);
  unsigned coordsType = BIQUADR_FE; // get the finite element type for "x", it is always 2 (LAGRANGE TENSOR-PRODUCT-QUADRATIC)
  for(int i = 0; i < dim; i++) {   
           coords[i].reserve(max_size_elem_dofs); 
       coords_hat[i].reserve(max_size_elem_dofs); 
 }
  // geometry ------------------------------------------
  

  //------------ at quadrature points ---------------------
  real_num_mov weight_qp = 0.;
    double weight_hat_qp = 0.;
    vector < real_num > SolVAR_qp(n_unknowns);
    vector < vector < real_num > > gradSolVAR_qp(n_unknowns);
    vector < vector < real_num > > gradSolVAR_hat_qp(n_unknowns);
    for(int k = 0; k < n_unknowns; k++) { 
          gradSolVAR_qp[k].resize(dim);
      gradSolVAR_hat_qp[k].resize(dim);
    }
      

   // ------------------------------------------------------------------------
    // Physical parameters
    const int    solid_model  = ml_prob.parameters.get < Solid>("Solid").get_physical_model();
    const double mu_lame 	  = ml_prob.parameters.get < Solid>("Solid").get_lame_shear_modulus();
    const double lambda_lame  = ml_prob.parameters.get < Solid>("Solid").get_lame_lambda();

    const bool incompressible = (0.5  ==  ml_prob.parameters.get < Solid>("Solid").get_poisson_coeff()) ? 1 : 0;
    const bool penalty = ml_prob.parameters.get < Solid>("Solid").get_if_penalty();

    // gravity
    double _gravity[3] = {1., 0., 0.};
    // -----------------------------------------------------------------
 
    RES->zero();
  if (assembleMatrix) JAC->zero();
  
  adept::Stack & stack = FemusInit::_adeptStack;  // call the adept stack object for potential use of AD
  const assemble_jacobian< real_num, double/*real_num_mov*/ > assemble_jac;


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
  
  // geometry end *****************************

      
     assemble_jac.prepare_before_integration_loop(stack);
     
  
    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][SolFEType[sol_index_displ]]->GetGaussPointNumber(); ig++) {

	// *** get Jacobian and test function and test function derivatives ***
      for(int fe=0; fe < NFE_FAMS; fe++) {
	ml_prob._ml_msh->_finiteElement[ielGeom][fe]->Jacobian(coords,    ig,weight_qp,    phi_dof_qp[fe],    phi_x_dof_qp[fe],    phi_xx_dof_qp[fe]);
	ml_prob._ml_msh->_finiteElement[ielGeom][fe]->Jacobian(coords_hat,ig,weight_hat_qp,phi_hat_dof_qp[fe],phi_x_hat_dof_qp[fe],phi_xx_hat_dof_qp[fe]);
      }
         //HAVE TO RECALL IT TO HAVE BIQUADRATIC JACOBIAN
  	ml_prob._ml_msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coords,ig,weight_qp,phi_dof_qp[BIQUADR_FE],phi_x_dof_qp[BIQUADR_FE],phi_xx_dof_qp[BIQUADR_FE]);
  	ml_prob._ml_msh->_finiteElement[ielGeom][BIQUADR_FE]->Jacobian(coords_hat,ig,weight_hat_qp,phi_hat_dof_qp[BIQUADR_FE],phi_x_hat_dof_qp[BIQUADR_FE],phi_xx_hat_dof_qp[BIQUADR_FE]);


  //begin unknowns eval at gauss points ********************************
	for(unsigned unk = 0; unk <  n_unknowns; unk++) {
	  SolVAR_qp[unk] = 0.;

	  for(unsigned ivar2=0; ivar2<dim; ivar2++){ 
	    gradSolVAR_qp[unk][ivar2] = 0.; 
	    gradSolVAR_hat_qp[unk][ivar2] = 0.; 
	  }
	  
	  for(unsigned i = 0; i < Sol_n_el_dofs[unk]; i++) {
	    SolVAR_qp[unk]    += phi_dof_qp[ SolFEType[unk] ][i] * SolVAR_eldofs[unk][i];
	    for(unsigned ivar2=0; ivar2<dim; ivar2++) {
	      gradSolVAR_qp[unk][ivar2] += phi_x_dof_qp[ SolFEType[unk] ][i*dim+ivar2] * SolVAR_eldofs[unk][i]; 
	      gradSolVAR_hat_qp[unk][ivar2] += phi_x_hat_dof_qp[ SolFEType[unk] ][i*dim+ivar2] * SolVAR_eldofs[unk][i]; 
	    }
	  }
	  
	}  
 //end unknowns eval at gauss points ********************************


//    assemble_jacobian< real_num, real_num_mov >::mass_residual (Res, Sol_n_el_dofs, sum_Sol_n_el_dofs, SolPdeIndex, SolFEType, phi_dof_qp, SolVAR_qp, weight_hat_qp);
 

 //*******************************************************************************************************
   vector < vector < real_num_mov > > Cauchy(3); for (int i = 0; i < Cauchy.size(); i++) Cauchy[i].resize(3);
   real_num_mov J_hat;
   real_num_mov trace_e_hat;

    Cauchy = Solid::get_Cauchy_stress_tensor< real_num_mov >(solid_model, mu_lame, lambda_lame, dim, sol_index_press, gradSolVAR_hat_qp, SolVAR_qp, SolPdeIndex, J_hat, trace_e_hat);

    

              //BEGIN residual Solid Momentum in moving domain
          for (unsigned i = 0; i < Sol_n_el_dofs[sol_index_displ]; i++) {

              real_num_mov Cauchy_direction[3] = {0., 0., 0.};

              for (int idim = 0.; idim < dim; idim++) {
                for (int jdim = 0.; jdim < dim; jdim++) {
                  Cauchy_direction[idim] += phi_x_dof_qp[SolFEType[idim]][i * dim + jdim] * Cauchy[idim][jdim];
                }
              }

              for (int idim = 0; idim < dim; idim++) {
                Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[idim], i) ] += ( Cauchy_direction[idim] -  phi_dof_qp[SolFEType[idim]][i] * _gravity[idim] ) * weight_qp;
              }

            }
              //END residual Solid Momentum in moving domain
              

              //BEGIN residual solid mass balance in reference domain
            for (unsigned i = 0; i < Sol_n_el_dofs[sol_index_press]; i++) {
                
              Res[ assemble_jacobian<double,double>::res_row_index(Sol_n_el_dofs, SolPdeIndex[sol_index_press], i) ] += 
              weight_hat_qp * phi_dof_qp[SolFEType[sol_index_press]][i] * Solid::get_mass_balance_reference_domain< real_num_mov >(solid_model, penalty, incompressible, lambda_lame, trace_e_hat, J_hat, SolVAR_qp, SolPdeIndex, sol_index_press);
//               weight_qp * phi_dof_qp[SolFEType[sol_index_press]][i] * Solid::get_mass_balance_moving_domain< real_num_mov >(gradSolVAR_qp, SolPdeIndex);
                
            }
              //END residual solid mass balance in reference domain


    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------


    if (assembleMatrix) assemble_jac.compute_jacobian_outside_integration_loop(stack, SolVAR_eldofs, Res, Jac, L2G_dofmap_AllVars, RES, JAC);

    assemble_jacobian<real_num,real_num_mov>::print_element_residual(iel, Res, Sol_n_el_dofs, 9, 5);
    assemble_jacobian<real_num,real_num_mov>::print_element_jacobian(iel, Jac, Sol_n_el_dofs, 9, 5);
    
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
                                                                                const std::vector< Math::Unknown > &  unknowns,  
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

  // ======= Problem ========================
  MultiLevelProblem ml_prob(&ml_sol);

  //material  ==================
              //Nondimensional quantity (Lref,Uref)
            double Lref = 1.;
            double Uref = 1.;
            Parameter par(Lref,Uref);
            
            // Generate Solid Object
            const double E = 1500000;
            const double ni = 0.5;
            const double rhos = 1000;
            const Solid solid(par, E, ni, rhos, MODEL);

            std::cout << "Solid properties: " << std::endl;
            std::cout << solid << std::endl;
            
            ml_prob.parameters.set<Solid>("Solid") = solid;
  //end material  ==================

  ml_prob.SetFilesHandler(&files);

  // ======= Solution, II ==================

  for (unsigned int u = 0; u < unknowns.size(); u++)  ml_sol.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order);

  //initial conditions
  ml_sol.Initialize("All");
  for (unsigned int u = 0; u < unknowns.size(); u++)  ml_sol.Initialize(unknowns[u]._name.c_str(), Solution_set_initial_conditions, &ml_prob);
  
  //boundary conditions
  ml_sol.AttachSetBoundaryConditionFunction(Solution_set_boundary_conditions);
  ml_sol.GenerateBdc("All");


  // ======= System ========================
  NonLinearImplicitSystem& system = ml_prob.add_system < NonLinearImplicitSystem > ("SolidMech");

  for (unsigned int u = 0; u < unknowns.size(); u++)   system.AddSolutionToSystemPDE(unknowns[u]._name.c_str());
 
  ml_prob.set_unknown_list_for_assembly(unknowns); 
            
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

  ml_sol.SetWriter(VTK);
  ml_sol.GetWriter()->SetDebugOutput(true);
  system.SetDebugNonlinear(true);

//   system.SetMaxNumberOfNonLinearIterations(2);
//   system.SetNonLinearConvergenceTolerance(1.e-30);
//   system.SetDebugLinear(true);
//   system.SetMaxNumberOfLinearIterations(4);
//   system.SetAbsoluteLinearConvergenceTolerance(1.e-10);
 
  system.MGsolve();

  // ======= Print ========================
  // set moving variables
  std::vector < std::string > mov_vars;
  mov_vars.push_back("DX");
  mov_vars.push_back("DY");
  if ( ml_mesh.GetDimension() == 3 ) mov_vars.push_back("DZ");
  ml_sol.GetWriter()->SetMovingMesh(mov_vars);
  
  // print solutions
  std::vector < std::string > variablesToBePrinted;  variablesToBePrinted.push_back("All");
  ml_sol.GetWriter()->Write(files.GetOutputPath(),"biquadratic", variablesToBePrinted, lev);
 
 return ml_sol;

}

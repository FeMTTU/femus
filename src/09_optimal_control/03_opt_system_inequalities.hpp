#ifndef __03_OPT_SYSTEM_INEQUALITIES_HPP__
#define __03_OPT_SYSTEM_INEQUALITIES_HPP__
 

#include "NonLinearImplicitSystemWithPrimalDualActiveSetMethod.hpp"
#include "LinearEquationSolver.hpp"

//*********************** Domain and Mesh Independent - BEGIN *****************************************
  
//*********************** Inequality, in equations - BEGIN *******************************************************
#define  INEQ_FLAG 1
#define  C_COMPL 1.
//*********************** Inequality, in equations - END *******************************************************


namespace femus {




   
namespace ctrl {
    


template < class INEQUALITY_CONSTRAINT_CLASS >
class mixed_state_or_ctrl_inequality {

public:
    
 static void update_active_set_flag_for_current_nonlinear_iteration(const femus::Mesh* msh,
                                                             const femus::Solution* sol,
                                                             const unsigned int iel,
                                                             const std::vector < std::vector < double > > & coords_at_dofs,
                                                             const std::vector < std::vector < double > > sol_eldofs,  ///@todo why not a reference?
                                                             const std::vector < unsigned int > & Sol_n_el_dofs,
                                                             const double c_compl,
                                                             const std::vector<unsigned int>  pos_mu,
                                                             const std::vector<unsigned int>  pos_ctrl,
                                                             const std::vector<unsigned int>  solIndex_act_flag,
                                                                 std::vector <   std::vector < double > > & ctrl_lower_dofs,
                                                                 std::vector <   std::vector < double > > & ctrl_upper_dofs,
                                                                 std::vector <   std::vector < double > > & sol_actflag_dofs) {

        const unsigned int dim = coords_at_dofs.size();

    const unsigned int   n_components_ctrl = pos_ctrl.size();


              for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {


// 0: inactive; 1: active_a; 2: active_b
        assert(Sol_n_el_dofs[pos_mu[kdim] ] == Sol_n_el_dofs[ pos_ctrl[kdim] ]);
        sol_actflag_dofs[kdim].resize(Sol_n_el_dofs[pos_mu[kdim] ]);
        ctrl_lower_dofs[kdim].resize(Sol_n_el_dofs[pos_mu[kdim] ]);
        ctrl_upper_dofs[kdim].resize(Sol_n_el_dofs[pos_mu[kdim] ]);
        std::fill(sol_actflag_dofs[kdim].begin(), sol_actflag_dofs[kdim].end(), 0);
        std::fill(ctrl_lower_dofs[kdim].begin(),   ctrl_lower_dofs[kdim].end(), 0.);
        std::fill(ctrl_upper_dofs[kdim].begin(),   ctrl_upper_dofs[kdim].end(), 0.);

// // //         std::cout << " mu dofs " << std::endl;
// // //                 for (unsigned i = 0; i < sol_actflag_dofs.size(); i++) {
// // //                 std::cout << sol_eldofs[pos_mu][i] << " ";
// // //         }
// // //
// // //         std::cout << std::endl;

// // //         std::cout << " ctrl dofs " << std::endl;
// // //         for (unsigned i = 0; i < sol_actflag_dofs.size(); i++) {
// // //                 std::cout << sol_eldofs[pos_ctrl][i] << " ";
// // //         }
// // //         std::cout << std::endl;

        for (unsigned i = 0; i < sol_actflag_dofs[kdim].size(); i++) {
            std::vector<double> node_coords_i(dim, 0.);
            for (unsigned d = 0; d < dim; d++) node_coords_i[d] = coords_at_dofs[d][i];

            ctrl_lower_dofs[kdim][i] = /*femus:: ctrl::*/ INEQUALITY_CONSTRAINT_CLASS ::InequalityConstraint(n_components_ctrl, node_coords_i, false)[kdim];
            ctrl_upper_dofs[kdim][i] = /*femus:: ctrl::*/ INEQUALITY_CONSTRAINT_CLASS ::InequalityConstraint(n_components_ctrl, node_coords_i, true)[kdim];

            const double lower_test_value = sol_eldofs[ pos_mu[kdim] ][i] + c_compl * ( sol_eldofs[ pos_ctrl[kdim] ][i] - ctrl_lower_dofs[kdim][i] );
            const double upper_test_value = sol_eldofs[ pos_mu[kdim] ][i] + c_compl * ( sol_eldofs[ pos_ctrl[kdim] ][i] - ctrl_upper_dofs[kdim][i] );

            if      ( lower_test_value < 0. )  {
//                 std::cout << "Found active node below" << std::endl;
//                 std::cout << "The current value of mu is " <<  sol_eldofs[ pos_mu[kdim] ][i] << std::endl;
                   sol_actflag_dofs[kdim][i] = 1;
            }
            else if ( upper_test_value > 0. )  {
//                 std::cout << "Found active node above" << std::endl;
//                 std::cout << "The current value of mu is " <<  sol_eldofs[ pos_mu[kdim] ][i] << std::endl;
                sol_actflag_dofs[kdim][i] = 2;
            }
        }

//************** local to global act flag ***************************
     const unsigned solFEType_act_flag = sol->GetSolutionType(solIndex_act_flag[kdim]);

        for (unsigned i = 0; i < sol_actflag_dofs[kdim].size(); i++) {
            unsigned solDof_mu = msh->GetSolutionDof(i, iel, solFEType_act_flag);
            (sol->_Sol[ solIndex_act_flag[kdim] ])->set(solDof_mu, sol_actflag_dofs[kdim][i]);
        }


      }


   }





 static void update_active_set_flag_for_current_nonlinear_iteration_bdry(const femus::Mesh* msh,
                                                             const femus::Solution* sol,
                                                             const unsigned int iel,
                                                             const unsigned int iface,
                                                             const std::vector < std::vector < double > > & coords_at_dofs,
                                                             const std::vector < std::vector < double > > sol_eldofs,  ///@todo why not a reference?
                                                             const std::vector < unsigned int > & Sol_n_el_dofs,
                                                             const double c_compl,
                                                             const std::vector<unsigned int> pos_mu,
                                                             const std::vector<unsigned int> pos_ctrl,
                                                             const std::vector<unsigned int>  solIndex_act_flag,
                                                              std::vector < std::vector < double > > & ctrl_lower_dofs,
                                                              std::vector < std::vector < double > > & ctrl_upper_dofs,
                                                              std::vector < std::vector < double > > & sol_actflag_dofs) {

         const unsigned dim = coords_at_dofs.size();

      const unsigned int   n_components_ctrl = pos_ctrl.size();

       for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {



     const unsigned solFEType_act_flag = sol->GetSolutionType(solIndex_act_flag[kdim]);

     const unsigned ndofs_bdry = msh->GetElementFaceDofNumber(iel, iface, solFEType_act_flag);


         // 0: inactive; 1: active_a; 2: active_b
        assert(Sol_n_el_dofs[ pos_mu[kdim] ] == Sol_n_el_dofs[ pos_ctrl[kdim] ]);///@todo More appropriately,
            sol_actflag_dofs[kdim].resize(ndofs_bdry/*nDof_mu*/);
            ctrl_lower_dofs[kdim].resize(ndofs_bdry/*nDof_mu*/);
            ctrl_upper_dofs[kdim].resize(ndofs_bdry/*nDof_mu*/);
           std::fill(sol_actflag_dofs[kdim].begin(), sol_actflag_dofs[kdim].end(), 0);
           std::fill(ctrl_lower_dofs[kdim].begin(), ctrl_lower_dofs[kdim].end(), 0.);
           std::fill(ctrl_upper_dofs[kdim].begin(), ctrl_upper_dofs[kdim].end(), 0.);


      for (unsigned int i_bdry = 0; i_bdry < sol_actflag_dofs[kdim].size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
        std::vector<double> node_coords_i(dim, 0.);
        for (unsigned d = 0; d < dim; d++) node_coords_i[d] = coords_at_dofs[d][i_bdry];

        ctrl_lower_dofs[kdim][i_bdry] = /*femus:: ctrl::*/ INEQUALITY_CONSTRAINT_CLASS ::InequalityConstraint(n_components_ctrl, node_coords_i, false)[kdim];
        ctrl_upper_dofs[kdim][i_bdry] = /*femus:: ctrl::*/ INEQUALITY_CONSTRAINT_CLASS ::InequalityConstraint(n_components_ctrl, node_coords_i, true)[kdim];

        const double lower_test_value = sol_eldofs[ pos_mu[kdim] ][i_vol] + c_compl * ( sol_eldofs[ pos_ctrl[kdim] ][i_vol] - ctrl_lower_dofs[kdim][i_bdry] );
        const double upper_test_value = sol_eldofs[ pos_mu[kdim] ][i_vol] + c_compl * ( sol_eldofs[ pos_ctrl[kdim] ][i_vol] - ctrl_upper_dofs[kdim][i_bdry] );

        if      ( lower_test_value < 0. )  sol_actflag_dofs[kdim][i_bdry] = 1;
        else if ( upper_test_value > 0. )  sol_actflag_dofs[kdim][i_bdry] = 2;
            }

//************** local to global act flag ***************************
      for (unsigned int i_bdry = 0; i_bdry < sol_actflag_dofs[kdim].size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
      unsigned solDof_actflag = msh->GetSolutionDof(i_vol, iel, solFEType_act_flag);
      (sol->_Sol[solIndex_act_flag[kdim]])->set(solDof_actflag, sol_actflag_dofs[kdim][i_bdry]);
           }

       }


    }




 static void node_insertion(const unsigned int iel,
                         const   Mesh* msh,
                         const     std::vector < std::vector < int > > & L2G_dofmap_Mat,
                                                             const std::vector<unsigned int> pos_mat_mu,
                                                             const std::vector<unsigned int> pos_mat_ctrl,
                         const std::vector < std::vector < double > > & sol_eldofs_Mat,
                         const std::vector < unsigned int > & Sol_n_el_dofs,
                         std::vector < std::vector < double > > & sol_actflag_dofs,
                         const std::vector < std::vector < double > > & ctrl_lower_dofs,
                         const std::vector < std::vector < double > > & ctrl_upper_dofs,
                         const double ineq_flag,
                         const double c_compl,
                         SparseMatrix*             KK,
                         NumericVector* RES,
                          const bool assembleMatrix
                        ) {


   const unsigned int   n_components_ctrl = pos_mat_ctrl.size();

       for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {



      //***************************************************
    std::vector < int > l2GMap_mu(sol_actflag_dofs[kdim].size());
    std::vector < int > l2GMap_ctrl(sol_actflag_dofs[kdim].size());
    for (unsigned i = 0; i < sol_actflag_dofs[kdim].size(); i++) {
      l2GMap_mu[i]   = L2G_dofmap_Mat[ pos_mat_mu[kdim]  ][i];   //pdeSys->GetSystemDof(solIndex_mu, solPdeIndex_mu, i, iel);
      l2GMap_ctrl[i] = L2G_dofmap_Mat[ pos_mat_ctrl[kdim] ][i]; //pdeSys->GetSystemDof(solIndex_ctrl, solPdeIndex_ctrl, i, iel);
    }
 //***************************************************

 //============= delta_mu row - BEGIN ===============================
      std::vector<double> Res_mu (sol_actflag_dofs[kdim].size()); std::fill(Res_mu.begin(),Res_mu.end(), 0.);

    for (unsigned i = 0; i < sol_actflag_dofs[kdim].size(); i++) {
      if (sol_actflag_dofs[kdim][i] == 0){  //inactive
         Res_mu [i] = - ineq_flag * ( 1. * sol_eldofs_Mat[pos_mat_mu[kdim]][i] - 0. );
// 	 Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i];
      }
      else if (sol_actflag_dofs[kdim][i] == 1){  //active_a
	 Res_mu [i] = - ineq_flag * ( c_compl *  sol_eldofs_Mat[pos_mat_ctrl[kdim]][i] - c_compl * ctrl_lower_dofs[kdim][i]);
      }
      else if (sol_actflag_dofs[kdim][i] == 2){  //active_b
	Res_mu [i]  =  - ineq_flag * ( c_compl *  sol_eldofs_Mat[pos_mat_ctrl[kdim]][i] - c_compl * ctrl_upper_dofs[kdim][i]);
      }
    }
//          Res[nDof_u + nDof_ctrl + nDof_adj + i]  = c_compl * (  (2 - sol_actflag_dofs[i]) * (ctrl_lower_dofs[i] - sol_ctrl[i]) + ( sol_actflag_dofs[i] - 1 ) * (ctrl_upper_dofs[i] - sol_ctrl[i])  ) ;
//          Res_mu [i] = Res[nDof_u + nDof_ctrl + nDof_adj + i] ;


    RES->insert(Res_mu, l2GMap_mu);
//     RES->insert(Res_ctrl, l2GMap_ctrl);
//     RES->insert(Res_u, l2GMap_u);
//     RES->insert(Res_adj, l2GMap_adj);
 //============= delta_mu row - END ===============================

//  //============= delta_state-delta_state row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_u, l2GMap_u, 1.);

//  //============= delta_ctrl-delta_ctrl row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_ctrl, l2GMap_ctrl, 1.);

//  //============= delta_adj-delta_adj row ===============================
//  KK->matrix_set_off_diagonal_values_blocked(l2GMap_adj, l2GMap_adj, 1.);

 //============= delta_mu-delta_ctrl row  - BEGIN ===============================
 for (unsigned i = 0; i < sol_actflag_dofs[kdim].size(); i++) if (sol_actflag_dofs[kdim][i] != 0 ) sol_actflag_dofs[kdim][i] = ineq_flag * c_compl;

  if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_ctrl, sol_actflag_dofs[kdim]); }
 //============= delta_mu-delta_ctrl row  - END ===============================

 //============= delta_mu-delta_mu row - BEGIN ===============================
  for (unsigned i = 0; i < sol_actflag_dofs[kdim].size(); i++) sol_actflag_dofs[kdim][i] =   ineq_flag * (1 - sol_actflag_dofs[kdim][i]/c_compl)  + (1-ineq_flag) * 1.;  //can do better to avoid division, maybe use modulo operator

  if (assembleMatrix) {    KK->matrix_set_off_diagonal_values_blocked(l2GMap_mu, l2GMap_mu, sol_actflag_dofs[kdim] );  }
 //============= delta_mu-delta_mu row - END ===============================


       }//kdim


 }


 
 static void node_insertion_bdry(const unsigned int iel,
                         const unsigned int iface,
                         const   Mesh* msh,
                         const     std::vector < std::vector < int > > & L2G_dofmap,
                                                             const std::vector<unsigned int> pos_mu,
                                                             const std::vector<unsigned int> pos_ctrl,
                         const std::vector < std::vector < double > > & sol_eldofs,
                         const std::vector < unsigned int > & Sol_n_el_dofs,
                         std::vector < std::vector < double > > & sol_actflag_dofs,
                         const std::vector < std::vector < double > > & ctrl_lower_dofs,
                         const std::vector < std::vector < double > > & ctrl_upper_dofs,
                         const double ineq_flag,
                         const double c_compl,
                         SparseMatrix*             KK,
                         NumericVector* RES,
                          const bool assembleMatrix
                        ) {


   const unsigned int   n_components_ctrl = pos_ctrl.size();

       for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {



// Create the L2G boundary maps from the volume ones
  std::vector < int > L2G_dofmap_mu_bdry(sol_actflag_dofs[kdim].size());
  std::vector < int > L2G_dofmap_ctrl_bdry(sol_actflag_dofs[kdim].size());

      for (int i_bdry = 0; i_bdry < sol_actflag_dofs[kdim].size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
  L2G_dofmap_mu_bdry[i_bdry]   = L2G_dofmap[ pos_mu[kdim]  ][i_vol];
  L2G_dofmap_ctrl_bdry[i_bdry] = L2G_dofmap[ pos_ctrl[kdim] ][i_vol];
      }

 //============= delta_mu row - BEGIN ===============================
      std::vector<double> Res_mu_bdry (sol_actflag_dofs[kdim].size());     std::fill(Res_mu_bdry.begin(),Res_mu_bdry.end(), 0.);
//       std::vector<double> Res_mu (Sol_n_el_dofs[pos_mu]);       std::fill(Res_mu.begin(),Res_mu.end(), 0.);

      for (int i_bdry = 0; i_bdry < sol_actflag_dofs[kdim].size(); i_bdry++)  {
	    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);

      if (sol_actflag_dofs[kdim][i_bdry] == 0) {  //inactive
//          Res_mu [i_vol]      = - ineq_flag * ( 1. * sol_eldofs[pos_mu[kdim]][i_vol] - 0. );
         Res_mu_bdry[i_bdry] = - ineq_flag * ( 1. * sol_eldofs[pos_mu[kdim]][i_vol] - 0. );
      }
      else if (sol_actflag_dofs[kdim][i_bdry] == 1) {  //active_a
// 	 Res_mu [i_vol]      = - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl[kdim]][i_vol] - c_compl * ctrl_lower_dofs[kdim][i_bdry]);
     Res_mu_bdry[i_bdry] = - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl[kdim]][i_vol] - c_compl * ctrl_lower_dofs[kdim][i_bdry]);

    }
      else if (sol_actflag_dofs[kdim][i_bdry] == 2) {  //active_b
// 	Res_mu [i_vol]      =  - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl[kdim]][i_vol] - c_compl * ctrl_upper_dofs[kdim][i_bdry]);
    Res_mu_bdry[i_bdry] =  - ineq_flag * ( c_compl *  sol_eldofs[pos_ctrl[kdim]][i_vol] - c_compl * ctrl_upper_dofs[kdim][i_bdry]);
      }
    }


//     RES->insert(Res_mu,  L2G_dofmap[pos_mu]);
    RES->insert(Res_mu_bdry,  L2G_dofmap_mu_bdry);
 //============= delta_mu row - END ===============================

 //============= delta_mu-delta_ctrl row  - BEGIN ===============================
 //auxiliary volume vector for act flag
//  unsigned nDof_actflag_vol  = msh->GetElementDofNumber(iel, solFEType_act_flag);
//  std::vector<double> sol_actflag_dofs_vol(nDof_actflag_vol);


 for (unsigned i_bdry = 0; i_bdry < sol_actflag_dofs[kdim].size(); i_bdry++) if (sol_actflag_dofs[kdim][i_bdry] != 0 ) sol_actflag_dofs[kdim][i_bdry] = ineq_flag * c_compl;

//  std::fill(sol_actflag_dofs_vol.begin(), sol_actflag_dofs_vol.end(), 0.);
//     for (int i_bdry = 0; i_bdry < sol_actflag_dofs.size(); i_bdry++)  {
//        unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
//        sol_actflag_dofs_vol[i_vol] = sol_actflag_dofs[i_bdry];
//     }

//  KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_ctrl], sol_actflag_dofs_vol);
 if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap_mu_bdry, L2G_dofmap_ctrl_bdry, sol_actflag_dofs[kdim]); }
 //============= delta_mu-delta_ctrl row - END ===============================

 //============= delta_mu-delta_mu row - BEGIN ===============================
 // Attention: this equation goes in contrast with \mu = 0 on \Omega \setminus \Gamma_c
 // In fact, here we shouldn't insert all VOLUME values, but only the BOUNDARY ones
 // The best way is to then do a L2G_map ON THE BOUNDARY only

  for (unsigned i_bdry = 0; i_bdry < sol_actflag_dofs[kdim].size(); i_bdry++) sol_actflag_dofs[kdim][i_bdry] =  ineq_flag * (1 - sol_actflag_dofs[kdim][i_bdry]/c_compl)  + (1-ineq_flag) * 1.;  //can do better to avoid division, maybe use modulo operator

//  std::fill(sol_actflag_dofs_vol.begin(), sol_actflag_dofs_vol.end(), 0.);
//     for (int i_bdry = 0; i_bdry < sol_actflag_dofs.size(); i_bdry++)  {
//        unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
//        sol_actflag_dofs_vol[i_vol] = sol_actflag_dofs[i_bdry];
//     }

//   KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap[pos_mu], L2G_dofmap[pos_mu], sol_actflag_dofs_vol );
  if (assembleMatrix) { KK->matrix_set_off_diagonal_values_blocked(L2G_dofmap_mu_bdry, L2G_dofmap_mu_bdry, sol_actflag_dofs[kdim]);  }
 //============= delta_mu-delta_mu row - END ===============================

       }//kdim


}




///@todo This is being added to a weak form?
///this is the same for volume or boundary
static void add_one_times_mu_res_ctrl(const unsigned iproc,
                         const double ineq_flag,
                         const std::vector<unsigned int> pos_ctrl_in_Mat,
                         const std::vector<unsigned int> pos_mu_in_Mat,
                         const std::vector < unsigned > & SolIndex,
                         const Solution*                sol,
                         const femus::NonLinearImplicitSystemWithPrimalDualActiveSetMethod * mlPdeSys,
                         const  LinearEquationSolver* pdeSys,
                         NumericVector* RES) {

     assert(pos_ctrl_in_Mat.size() == pos_mu_in_Mat.size());

    const unsigned int   n_components_ctrl = pos_ctrl_in_Mat.size();

       for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {

 const unsigned int ctrl_index_in_Mat = pos_ctrl_in_Mat[kdim];
 const unsigned int mu_index_in_Mat   = pos_mu_in_Mat[kdim];
 const unsigned int mu_index_in_Sol   = pos_mu_in_Mat[kdim];

  unsigned int ctrl_size_iproc_in_Mat = pdeSys->KKoffset[ctrl_index_in_Mat + 1][iproc] - pdeSys->KKoffset[ctrl_index_in_Mat][iproc];
  unsigned int mu_size_iproc_in_Sol = (*sol->_Sol[ SolIndex[mu_index_in_Sol] ]).last_local_index() - (*sol->_Sol[ SolIndex[mu_index_in_Sol] ]).first_local_index(); // pdeSys->KKoffset[mu_index_in_Mat + 1][iproc] - pdeSys->KKoffset[mu_index_in_Mat][iproc];

  assert(ctrl_size_iproc_in_Mat == mu_size_iproc_in_Sol);

  std::vector<double>  one_times_mu(ctrl_size_iproc_in_Mat, 0.);
  std::vector<int>    positions_ctrl_in_Res(ctrl_size_iproc_in_Mat);
  std::vector<int>    positions_mu_in_Sol(ctrl_size_iproc_in_Mat);

  for (unsigned i = 0; i < positions_ctrl_in_Res.size(); i++) {
    positions_ctrl_in_Res[i] = pdeSys->KKoffset[ctrl_index_in_Mat][iproc] + i;
    positions_mu_in_Sol[i]   = (*sol->_Sol[ SolIndex[mu_index_in_Sol] ]).first_local_index()/*pdeSys->KKoffset[mu_index_in_Mat][iproc]*/ + i;
    //this should not come from pdeSys but from Sol//actually I can take it from the Numeric Vector!
    ///@todo put the Dof range for Sol
//          unsigned solDof = msh->GetSolutionDof(i, iel, SolFEType[k]);  //this needs iel, in fact it is only in an ELEMENT loop, but here I am in a NODE loop

    one_times_mu[i] = ineq_flag * 1. * (*sol->_Sol[ SolIndex[mu_index_in_Sol] ])(positions_mu_in_Sol[i]) ;
    }

    RES->add_vector_blocked(one_times_mu, positions_ctrl_in_Res);


       }


 }



 static void store_act_flag_in_old(  const NonLinearImplicitSystemWithPrimalDualActiveSetMethod* mlPdeSys,
                               const MultiLevelSolution *    ml_sol,
                              Solution *                sol,
                             std::vector<unsigned int>  & solIndex_act_flag) {

     const unsigned int   n_components_ctrl = mlPdeSys->GetActiveSetFlagName().size();

          for (unsigned kdim = 0; kdim < n_components_ctrl; kdim++) {

  const std::string act_flag_name = mlPdeSys->GetActiveSetFlagName()[kdim];

  solIndex_act_flag[kdim]  = ml_sol->GetIndex(act_flag_name.c_str());

     if(sol->GetSolutionTimeOrder(solIndex_act_flag[kdim]) == 2) {
       *(sol->_SolOld[solIndex_act_flag[kdim]]) = *(sol->_Sol[solIndex_act_flag[kdim]]);
     }

  }

  }

};


//*********************** Domain and Mesh Independent - END *****************************************

}


}



#endif

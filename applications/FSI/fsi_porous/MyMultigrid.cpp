#include "MyMultigrid.hpp"
#include "main.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelMesh.hpp"
#include "ElemType.hpp"
#include "LinearEquationSolver.hpp"
#include "Elem.hpp"
#include "NumericVector.hpp"
#include "PetscVector.hpp"
#include "SparseMatrix.hpp"
#include "PetscMatrix.hpp"
#include "Parameter.hpp"

// // // MyMultiGrid::MyMultiGrid(MultiLevelMesh * mlmesh_in ) 
// // // :   MultiLevelProblem(mlmesh_in) {  
// // //       
// // //       _qty_idx[0] = IDX_P;
// // //       _qty_idx[1] = IDX_D;
// // //       _qty_idx[2] = IDX_VEL;
// // //       
// // //      _qty_ncomps[0] = NVAR_P; 
// // //      _qty_ncomps[1] = NVAR_D; 
// // //      _qty_ncomps[2] = NVAR_VEL; 
// // //       
// // //     }
// // // 
// // // 
// // //     
// // //   double MyMultiGrid::IncreaseStep(unsigned i)  {
// // //     
// // //     // Now time_step is no more a variable of MultiLevelProblem but of TransientSystem
// // // //    _time_step += i; 
// // //         
// // //   }
// // //     
// // //     
// // //  void MyMultiGrid::AddParameters(FemusInputParser<double>* runtimein) {
// // //      
// // //      _runtime_double = runtimein;
// // //      return;
// // //   }
// // //  
// // //  
// // //     
// // // //ok, I need to compute the integral of some function on some domain.
// // // //The point is how to specify the function
// // // //i need to loop over the boundary elements of face 1
// // // // then i need to compute the gradient of pressure at that face
// // // 
// // // //i might either build a node field gradP from the node field p,
// // // // and then use this, or do a volume loop, compute gradP at the volume nodes,
// // // //which are also boundary nodes, and compute at the surface gauss points
// // // //ok, so we must compute at the boundary gauss values 
// // // //The point is: GIVEN A LAGRANGIAN QUADRATIC FIELD,
// // // // HOW CAN I COMPUTE THE GRADIENT OF THAT FIELD AT THE NODES
// // // //IF I HAVE THAT THE NODE SHAPE FUNCTIONS HAVE DISCONTINUOUS DERIVATIVES?!
// // // // THE DOFS OF THAT FIELD ARE ALSO COMPUTED IN WEAK FORM
// // // // FROM THE RESPECTIVE EQUATION!
// // // 
// // // // AAA this routine needs the ID of the FACE, so it depends on the MESH FILE!!!
// // // 
// // // double MyMultiGrid::ComputeProductivityIndexNum(int bd) {
// // // 
// // // // // //   double Integral = 0.;
// // // // // //   double Integral_measure = 0.;
// // // // // // 
// // // // // //   PetscErrorCode ierr;
// // // // // // 
// // // // // // // Physics  
// // // // // //   const double rho  = _fluid->get_density();
// // // // // //   const double mu   = _fluid->get_viscosity();
// // // // // //   const double uref = _fluid->_parameter.Get_reference_velocity();
// // // // // // 
// // // // // // //VARIABLES involved  
// // // // // //   std::string varnames_u[NVAR_VEL];
// // // // // //   varnames_u[0] = "UX";
// // // // // //   varnames_u[1] = "UY";
// // // // // //   varnames_u[2] = "UZ";
// // // // // //   unsigned ku=GetIndex("UX");
// // // // // // 
// // // // // //    int indexSolD[NVAR_D];
// // // // // //   std::string varnamesD[NVAR_D];
// // // // // //   varnamesD[0] = "DX";
// // // // // //   varnamesD[1] = "DY";
// // // // // //   varnamesD[2] = "DZ";
// // // // // // 
// // // // // //    for (int i=0; i<NVAR_D; ++i)    indexSolD[i] = GetIndex(varnamesD[i].c_str());  
// // // // // //     
// // // // // // // parameters
// // // // // //    int moving_dom_ints = _runtime_double->get("moving_dom_ints");
// // // // // // 
// // // // // //     
// // // // // // // face types  
// // // // // //   const short unsigned NV1[6][2] = {{9,9},{6,6},{9,6},{3,3},{3,3},{1,1}}; //this is only for the faces, it doesnt work in 3D
// // // // // //   const unsigned      FELT[6][2] = {{3,3},{4,4},{3,4},{5,5},{5,5}}; //its length is 5 not 6!
// // // // // // 
// // // // // // // real geometry
// // // // // //   double vx[SPACEDIM][27];
// // // // // //   double normal[3];
// // // // // //   PetscInt node_geom[27];
// // // // // // 
// // // // // // //FE  
// // // // // //   double phi2[27],gradphi2[27][SPACEDIM];
// // // // // // 
// // // // // // //quadrature  
// // // // // //   double Weight2XJac;
// // // // // //   
// // // // // // //boundary integral
// // // // // //   for (unsigned lev=0; lev<gridn; lev++) {  //levels //here you have to loop over all levels, because you have coarser or finer elements
// // // // // // 
// // // // // //   LinearSolverM*    lsyspde_lev = Lin_Solver_[lev];
// // // // // //   unsigned order_ind2 = lsyspde_lev->SolType[ku];
// // // // // //   
// // // // // // // vectors for getting dofs   ====================
// // // // // // //here i am getting all the dofs of all the variables  
// // // // // //   int sol_size = lsyspde_lev->Sol_.size();
// // // // // //   PetscScalar **vec_sol     = new PetscScalar*[sol_size];
// // // // // //   vector<PetscVectorM*> petsc_vec_sol;
// // // // // //   petsc_vec_sol.resize(sol_size);
// // // // // // 
// // // // // //   for (int k=0; k<sol_size; k++) {
// // // // // //     petsc_vec_sol[k]  = static_cast<PetscVectorM*>(lsyspde_lev->Sol_[k]);
// // // // // //     ierr = VecGetArray(petsc_vec_sol[k]->vec(),&vec_sol[k]);    CHKERRQ(ierr);
// // // // // //   }
// // // // // // // vectors for getting dofs   ====================  
// // // // // // 
// // // // // // 
// // // // // //     for (unsigned iel=0; iel< lsyspde_lev->GetElementNumber(); iel++) {
// // // // // //       if ( lsyspde_lev->el->GetRefinedElementIndex(iel)==0 || lev==gridn-1u ) { //do i need this? this one is ok to avoid spanning twice domain portions //is this gridn or gridn - 1?
// // // // // // 
// // // // // //  
// // // // // //         short unsigned kelt = lsyspde_lev->el->GetElementType(iel);
// // // // // //         unsigned     nfaces = lsyspde_lev->el->GetElementFaceNumber(iel);
// // // // // // 
// // // // // //         for (unsigned jface=0; jface<nfaces; jface++) {
// // // // // // 	  
// // // // // // 	  int face_elem_index = lsyspde_lev->el->GetFaceElementIndex(iel,jface);
// // // // // // 	  
// // // // // //           if ( face_elem_index < 0 ) {  //what is dis?
// // // // // // 
// // // // // //             int gambit_face_number = - ( face_elem_index +1 ) ;
// // // // // // 	    
// // // // // //             if (gambit_face_number == bd) {  //EXTRACTING MY BOUNDARY FACES
// // // // // // 
// // // // // //           unsigned is_face_type_zero_or_one = ( jface<(lsyspde_lev->el->GetElementFaceNumber(iel,0)) );
// // // // // // 
// // // // // //               unsigned               nve =  NV1[kelt][is_face_type_zero_or_one];
// // // // // //               unsigned face_element_type = FELT[kelt][is_face_type_zero_or_one];
// // // // // // 
// // // // // //     for (unsigned i=0;i<nve;i++) {
// // // // // //       unsigned inode = lsyspde_lev->el->GetFaceVertexIndex(iel,jface,i)-1u; //WAS VOLUME: myel->GetElementDofIndex(kel,i)-1u;
// // // // // //       node_geom[i] = inode;
// // // // // //     for (unsigned j=0; j<SPACEDIM; j++)    {  vx[j][i]=vt[j][inode] + moving_dom_ints*vec_sol[ indexSolD[j] ][inode]; }
// // // // // //     }
// // // // // //        
// // // // // //               for (unsigned igs=0; igs < type_elem[face_element_type][order_ind2]->GetGaussPointNumber(); igs++) {
// // // // // // 		
// // // // // //                 (type_elem[face_element_type][order_ind2]->*type_elem[face_element_type][order_ind2]->Jacobian_sur_ptr)(vx,igs,Weight2XJac,phi2,gradphi2,normal);
// // // // // // 
// // // // // //        double SolU[NVAR_VEL]={0.,0.,0.};
// // // // // //         for (unsigned i=0; i<nve; i++) {
// // // // // //               for (unsigned d=0; d<NVAR_VEL; d++) {
// // // // // // 	       unsigned k = GetIndex( varnames_u[d].c_str() );
// // // // // //           double soli = vec_sol[k][node_geom[i]];
// // // // // //           SolU[d]+=phi2[i]*soli;
// // // // // // //           for (unsigned j=0; j<SPACEDIM; j++)  SolU[j]+= -gradphi2[i][j]*soli;  // drag term (- alpha*1/alpha grad P)
// // // // // //         }
// // // // // // 	  
// // // // // // 	}
// // // // // // 	
// // // // // // 	double udotn = 0.;
// // // // // // 	for (unsigned j=0; j<SPACEDIM; j++) udotn +=SolU[j]*normal[j];
// // // // // //                   Integral += udotn*Weight2XJac;
// // // // // //                   Integral_measure += 1.*Weight2XJac;
// // // // // // 
// // // // // //               } //end loop on gauss point
// // // // // // 
// // // // // //             }
// // // // // //           }
// // // // // //         }
// // // // // // 
// // // // // //       }
// // // // // //     }
// // // // // // 
// // // // // // 
// // // // // //     //free memory
// // // // // //   for (int k=0; k<lsyspde_lev->Sol_.size(); k++) {
// // // // // //     ierr = VecRestoreArray(petsc_vec_sol[k]->vec(),&vec_sol[k]);
// // // // // //     CHKERRQ(ierr);
// // // // // //     delete [] vec_sol[k];
// // // // // //   }
// // // // // //   delete [] vec_sol;
// // // // // //   
// // // // // // 
// // // // // //   } //level loop
// // // // // //   
// // // // // //   
// // // // // //    std::cout << "The integral in the numerator is " << Integral <<std::endl;
// // // // // //    std::cout << "The measure of the well is " << Integral_measure <<std::endl;
// // // // // // 
// // // // // //   return Integral;
// // // 
// // // }
// // //     
// // // 
// // // //list of involved variables, list of operators acting on those variables, operations that join the variables (sum,dot,cross...)
// // // 
// // // //quadrature rule to be used
// // // //type of domain (volume or boundary), regions in that domain
// // // //the point is that you might have CONSTANTS inside the integrand function
// // // 
// // // std::vector<double> MyMultiGrid::ComputeProductivityIndexDen(unsigned n_integrals) {
// // // 
// // // // // //   std::vector<double> Integrals(n_integrals);
// // // // // //  for (unsigned i=0;i<n_integrals;i++) { Integrals[i]=0.;  }
// // // // // //   
// // // // // // //VARIABLES involved  
// // // // // //   std::string varnames_p[NVAR_P];
// // // // // //   varnames_p[0] = "p";
// // // // // //   unsigned kp = GetIndex(varnames_p[0].c_str());
// // // // // // 
// // // // // //   int indexSolD[NVAR_D];
// // // // // //   std::string varnamesD[NVAR_D];
// // // // // //   varnamesD[0] = "DX";
// // // // // //   varnamesD[1] = "DY";
// // // // // //   varnamesD[2] = "DZ";
// // // // // // 
// // // // // //    for (int i=0; i<NVAR_D; ++i)    indexSolD[i] = GetIndex(varnamesD[i].c_str());  
// // // // // // 
// // // // // //   int moving_dom_ints = _runtime_double->get("moving_dom_ints");
// // // // // //    
// // // // // //   PetscErrorCode ierr;
// // // // // //   
// // // // // // // real geometry
// // // // // //   double vx[SPACEDIM][27];
// // // // // //   PetscInt node_geom[27];
// // // // // // 
// // // // // // //FE  
// // // // // //   double phi2[27],gradphi2[27][SPACEDIM];
// // // // // // 
// // // // // // //quadrature  
// // // // // //   int order_pick_quadrature = 0;
// // // // // //   double Weight2XJac;
// // // // // //   
// // // // // // 
// // // // // //   for (unsigned lev=0; lev<gridn; lev++) {  //levels //here you have to loop over all levels, because you have coarser or finer elements
// // // // // // 
// // // // // //   LinearSolverM*    lsyspde_lev = Lin_Solver_[lev];
// // // // // //   elem*    myel  = lsyspde_lev->el;
// // // // // //   unsigned order_ind2 = lsyspde_lev->SolType[kp];
// // // // // //   unsigned end_ind2   = lsyspde_lev->END_IND[order_ind2];
// // // // // // 
// // // // // // // vectors for getting dofs   ====================
// // // // // // //here i am getting all the dofs of all the variables  
// // // // // //   int sol_size = lsyspde_lev->Sol_.size();
// // // // // //   PetscScalar **vec_sol     = new PetscScalar*[sol_size];
// // // // // //   vector<PetscVectorM*> petsc_vec_sol;
// // // // // //   petsc_vec_sol.resize(sol_size);
// // // // // // 
// // // // // //   for (int k=0; k<sol_size; k++) {
// // // // // //     petsc_vec_sol[k]  = static_cast<PetscVectorM*>(lsyspde_lev->Sol_[k]);
// // // // // //     ierr = VecGetArray(petsc_vec_sol[k]->vec(),&vec_sol[k]);    CHKERRABORT(MPI_COMM_WORLD,ierr);
// // // // // //   }
// // // // // // // vectors for getting dofs   ====================  
// // // // // // 
// // // // // // 
// // // // // //     for (unsigned iel=0; iel< lsyspde_lev->GetElementNumber(); iel++) {
// // // // // //       
// // // // // //       if ( myel->GetRefinedElementIndex(iel)==0 || lev==gridn-1u ) { //do i need this? this one is ok to avoid spanning twice domain portions //is this gridn or gridn - 1?
// // // // // //  
// // // // // //      short unsigned kelt = myel->GetElementType(iel);
// // // // // //     unsigned         nve = myel->GetElementDofNumber(iel,end_ind2);
// // // // // // 
// // // // // //     for (unsigned i=0;i<nve;i++) {
// // // // // //       unsigned inode = myel->GetElementDofIndex(iel,i)-1u;
// // // // // //       node_geom[i] = inode;
// // // // // //     for (unsigned j=0; j<SPACEDIM; j++)    {  vx[j][i]=vt[j][inode] + moving_dom_ints*vec_sol[ indexSolD[j] ][inode]; }
// // // // // //     }
// // // // // //        
// // // // // //      for (unsigned igs=0; igs < type_elem[kelt][order_pick_quadrature]->GetGaussPointNumber(); igs++) {
// // // // // // 	  
// // // // // //             (type_elem[kelt][order_ind2]->*type_elem[kelt][order_ind2]->Jacobian_ptr)(vx,igs,Weight2XJac,phi2,gradphi2);  //i don't like that it is different from vol and surface, it is SUR or WITHOUT SUR
// // // // // // 
// // // // // //       double Solp[NVAR_P]={0.};
// // // // // //         for (unsigned i=0; i<nve; i++) {
// // // // // //               for (unsigned d=0; d<NVAR_P; d++) {
// // // // // // 	       unsigned k = GetIndex( varnames_p[d].c_str() );
// // // // // //           double soli = vec_sol[k][node_geom[i]];
// // // // // //           Solp[d]+=phi2[i]*soli;
// // // // // //         }
// // // // // // 	  
// // // // // // 	}
// // // // // // 	
// // // // // //                   Integrals[0]     += Solp[0]*Weight2XJac;  //here you should subtract p_at_the_well, which is ZERO for me
// // // // // //                   Integrals[1]     += Weight2XJac;
// // // // // // 
// // // // // //               } //end loop on gauss point
// // // // // // 
// // // // // //       }
// // // // // //     }
// // // // // // 
// // // // // // //  After the element loop is over, you have to subtract the pressure at the well surface, p_in... it is ZERO for us, at least so far...
// // // // // //     Integrals[0] -= _runtime_double->get("p_in");
// // // // // // 
// // // // // //     //free memory ====================
// // // // // //   for (int k=0; k<lsyspde_lev->Sol_.size(); k++) {
// // // // // //     ierr = VecRestoreArray(petsc_vec_sol[k]->vec(),&vec_sol[k]);
// // // // // //     CHKERRABORT(MPI_COMM_WORLD,ierr);
// // // // // //     delete [] vec_sol[k];
// // // // // //   }
// // // // // //   delete [] vec_sol;
// // // // // //     //free memory ===================
// // // // // //   
// // // // // // 
// // // // // //   } //level loop
// // // // // //   
// // // // // //   
// // // // // //   for (unsigned i=0;i<n_integrals;i++) {   std::cout << "The integral in the denominator is " << Integrals[i] <<std::endl;  }
// // // // // // 
// // // // // //   return Integrals;
// // // 
// // // }    
// // //     
// // //     
// // //     
// // //     
// // // double MyMultiGrid::Error() {
// // // // // //   double _Td=2;
// // // // // //   PetscInt node2[27];
// // // // // //   double vx[3][27];
// // // // // //   double phi2[27],gradphi2[27][3],Weight2;
// // // // // // 
// // // // // //   unsigned order_ind2 = Lin_Solver_[0]->SolType[GetIndex("T")];
// // // // // //   unsigned end_ind2   = Lin_Solver_[0]->END_IND[order_ind2];
// // // // // // 
// // // // // //   double Error=0;
// // // // // //   for (unsigned ig=0;ig<gridn;ig++) {
// // // // // //     PetscScalar *MYSOL[2];
// // // // // //     unsigned kT=GetIndex("T");
// // // // // //     PetscVectorM* petsc_vec_solT = static_cast<PetscVectorM*>(Lin_Solver_[ig]->Sol_[kT]);
// // // // // //     PetscErrorCode ierr;
// // // // // //     ierr = VecGetArray(petsc_vec_solT->vec(),&MYSOL[0]);
// // // // // //     CHKERRQ(ierr);
// // // // // // 
// // // // // //     unsigned kT0=GetIndex("T0");
// // // // // //     PetscVectorM* petsc_vec_solT0 = static_cast<PetscVectorM*>(Lin_Solver_[ig]->Sol_[kT0]);
// // // // // //     ierr = VecGetArray(petsc_vec_solT0->vec(),&MYSOL[1]);
// // // // // //     CHKERRQ(ierr);
// // // // // // 
// // // // // //     for (unsigned iel=0;iel<Lin_Solver_[ig]->GetElementNumber();iel++) {
// // // // // //       if (Lin_Solver_[ig]->el->GetRefinedElementIndex(iel)==0 || ig==gridn-1u) {
// // // // // //         unsigned nve2=Lin_Solver_[ig]->el->GetElementDofNumber(iel,end_ind2);
// // // // // //         double xe=0.,ye=0.,ze=0.;
// // // // // //         for (unsigned i=0;i<nve2;i++) {
// // // // // //           unsigned inode=Lin_Solver_[ig]->el->GetElementDofIndex(iel,i)-1u;
// // // // // //           node2[i]=inode;
// // // // // //           vx[0][i]=vt[0][inode];
// // // // // //           vx[1][i]=vt[1][inode];
// // // // // //           vx[2][i]=vt[2][inode];
// // // // // // 
// // // // // //           xe+=vx[0][i];
// // // // // //           ye+=vx[1][i];
// // // // // //           ze+=vx[2][i];
// // // // // //         }
// // // // // //         xe/=nve2;
// // // // // //         ye/=nve2;
// // // // // //         ze/=nve2;
// // // // // //         double testc=(xe>0.25 && xe<.75 && ye>0.25 && ye<.75)?1.:0.;
// // // // // //         if (testc) {
// // // // // //           short unsigned ielt=Lin_Solver_[ig]->el->GetElementType(iel);
// // // // // //           for (unsigned ig=0;ig < type_elem[ielt][order_ind2]->GetGaussPointNumber(); ig++) {
// // // // // //             // *** get Jacobian and test function and test function derivatives ***
// // // // // //             (type_elem[ielt][order_ind2]->*type_elem[ielt][order_ind2]->Jacobian_ptr)(vx,ig,Weight2,phi2,gradphi2);
// // // // // // 
// // // // // //             double SolT=0;
// // // // // //             double gradSolT[3]={0.,0.,0.};
// // // // // //             for (unsigned i=0; i<nve2; i++) {
// // // // // //               double soli=MYSOL[0][node2[i]];
// // // // // //               SolT+=phi2[i]*soli;
// // // // // //               gradSolT[0]+=gradphi2[i][0]*soli;
// // // // // //               gradSolT[1]+=gradphi2[i][1]*soli;
// // // // // //               gradSolT[2]+=gradphi2[i][2]*soli;
// // // // // //             }
// // // // // //             double SolT0=0;
// // // // // //             double gradSolT0[3]={0.,0.,0.};
// // // // // //             for (unsigned i=0; i<nve2; i++) {
// // // // // //               double soli=MYSOL[1][node2[i]];
// // // // // //               SolT0+=phi2[i]*soli;
// // // // // //               gradSolT0[0]+=gradphi2[i][0]*soli;
// // // // // //               gradSolT0[1]+=gradphi2[i][1]*soli;
// // // // // //               gradSolT0[2]+=gradphi2[i][2]*soli;
// // // // // //             }
// // // // // //             Error+=(_Td-(SolT+SolT0))*(_Td-(SolT+SolT0))*Weight2;
// // // // // //           }
// // // // // //         }
// // // // // //       }
// // // // // //     }
// // // // // // 
// // // // // //     ierr = VecRestoreArray(petsc_vec_solT->vec(),&MYSOL[0]);
// // // // // //     CHKERRQ(ierr);
// // // // // //     ierr = VecRestoreArray(petsc_vec_solT0->vec(),&MYSOL[1]);
// // // // // //     CHKERRQ(ierr);
// // // // // //   }
// // // // // // // 	for(unsigned j=0;j<Lin_Solver_[ig]->el->GetElementDofNumber(iel,index);j++) {
// // // // // // // 	  iconn = Lin_Solver_[ig]->el->GetElementDofIndex(iel,j)-1u;
// // // // // // //
// // // // // // 
// // // // // //   return Error;
// // // }
// // // 
// // // 
// // // 
// // // 
// // // //=================================
// // // //This function builds the sparsity pattern for the BIG MATRIX of the Problem.
// // // //I will loop over the single equations, make a "system element matrix", no need to do the same for the RHS,
// // // // and add all zeros
// // // 
// // // //how do I do that in my code? 
// // // //I first fill a Graph object, and then I pass it to the update sparsity pattern function
// // // //Now I want to do that by doing a "fake assembly loop"
// // // 
// // //  //for a given element, loop over the equations in that element, and add each of them to
// // // //every equation has several quantities and each of them has a different number of dofs 
// // // //but here we are using the number of dofs as a GEOMETRIC THING, and filling the node_dof with only the dofs we actually need
// // // //if you loop with the quantities inside the elements you'll put the elements outside and you'll pick an element only once
// // // //but in this way you'll have a different element matrix for every variable so you'll need to allocate these every time inside an element.
// // // //looping first with elements and then with quantities or the opposite does not change the number of operations,
// // // //but if SOME OPERATIONS are SLOWER than others, then THE ORDER MATTERS. So the order does not matter for the NUMBER OF OPERATIONS which is invariant, but for the SPEED.
// // // //putting some operations INSIDE or OUTSIDE might change things a lot...
// // // //yes, it is much better to put the QUANTITY LOOP OUTSIDE, so that ALLOCATIONS are performed OUTSIDE
// // // //maybe one could loop over the LIST OF ALL VARIABLES, and there are things which might be done in a shot.
// // // //One can make it later starting from this structure  
// // // 
// // // //now we are doing a QUANTITY LOOP but we should do an EQUATION LOOP, 
// // // // because the quantities are coupled so you have all the blocks to consider for the coupling
// // // 
// // // //Also, the SPARSITY PATTERN not only depends on the EQUATION, but also on the OPERATORS 
// // // //acting on that equation!!!
// // // //If you have a laplacian you don't have the same sparsity pattern as with the curl x curl !!!
// // // //So, one more reason to make the sparsity pattern be dependent on the EQUATION, and in particular the OPERATORS of the EQUATION!!!
// // // 
// // // //LU does not need the sparsity pattern because both L and U will be dense
// // // //ILU uses the SAME SPARSITY PATTERN as the initial matrix,
// // // //so you need the sparsity
// // // 
// // // void   MyMultiGrid::BuildSparsityPattern()  {
// // // 
// // // // // //   PetscErrorCode ierr;
// // // // // //   unsigned level_pick_order = 0;
// // // // // //   unsigned component_pick_SolIndex = 0;
// // // // // //  
// // // // // //   for (unsigned level=0; level<gridn; level++) {  //gridn belongs to the multigrid class
// // // // // //   
// // // // // // 
// // // // // //   LinearSolverM*     lsyspde_lev = Lin_Solver_[level];
// // // // // //   LinearSolverM* lsyspdemesh_lev = Lin_Solver_[level];
// // // // // //   
// // // // // //   Mat&           myKK      = lsyspde_lev->KK;
// // // // // //   vector <int>& myKKIndex  = lsyspde_lev->KKIndex;   //offset of the variables in the global matrix and rhs 
// // // // // //   Vec&           myRES     = lsyspde_lev->RES;
// // // // // // 
// // // // // //   ierr = MatZeroEntries(myKK);  CHKERRABORT(PETSC_COMM_WORLD,ierr);
// // // // // //   ierr = VecZeroEntries(myRES);  CHKERRABORT(PETSC_COMM_WORLD,ierr);
// // // // // //   
// // // // // // //         PetscViewer viewer;
// // // // // // //         ierr=PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);CHKERRQ(ierr);
// // // // // // //         ierr= MatView(myKK,viewer);CHKERRQ(ierr);
// // // // // // // //         ierr= VecView(myRES,viewer);CHKERRQ(ierr);
// // // // // // //         double ff;
// // // // // // //         std::cin>>ff;
// // // // // // 
// // // // // // //   PetscBool assem;
// // // // // // //   MatAssembled(myKK,&assem);
// // // // // // 
// // // // // //   // Element, geometry
// // // // // //   elem*    myel = lsyspdemesh_lev->el;
// // // // // //   unsigned nel  = lsyspdemesh_lev->GetElementNumber();
// // // // // //   PetscInt node_geom_el[MAX_EL_NODES]; 
// // // // // //   
// // // // // //   //begin equation loop
// // // // // //   // quantity loop
// // // // // //   // variables   ==============
// // // // // //   for (unsigned qty=0; qty< N_QTIES; qty++)  { 
// // // // // //   unsigned qty_idx_current = qty;
// // // // // //   
// // // // // //   unsigned ncomps = _qty_ncomps[qty_idx_current];
// // // // // // 
// // // // // //   std::vector<std::string> qty_comps_names(ncomps);
// // // // // // 
// // // // // //   if (qty==IDX_P) {
// // // // // //   qty_comps_names[0] = "p";
// // // // // //   }  
// // // // // //   if (qty==IDX_D) {
// // // // // //   qty_comps_names[0] = "DX";
// // // // // //   qty_comps_names[1] = "DY";
// // // // // //   qty_comps_names[2] = "DZ";
// // // // // //   }  
// // // // // //   if (qty==IDX_VEL) {
// // // // // //   qty_comps_names[0] = "UX";
// // // // // //   qty_comps_names[1] = "UY";
// // // // // //   qty_comps_names[2] = "UZ";
// // // // // //   }  
// // // // // //   
// // // // // // 
// // // // // //   unsigned qty_offset_prev_qties = 0;
// // // // // //   for (unsigned i=0; i< qty_idx_current; i++) { qty_offset_prev_qties += _qty_ncomps[i];  }
// // // // // //   
// // // // // //   PetscInt **  node_qty; node_qty = new PetscInt*[ncomps];
// // // // // //   for (unsigned i=0; i< ncomps; i++) node_qty[i] = new PetscInt[MAX_EL_NODES];
// // // // // //   std::vector<unsigned> index_qty_glob(ncomps);
// // // // // //   std::vector<unsigned> index_qty_localBlocks(ncomps);
// // // // // //   
// // // // // //   for (int i=0; i<ncomps; ++i) {   
// // // // // //     index_qty_glob[i]           = GetMGIndex(qty_comps_names[i].c_str());
// // // // // //     index_qty_localBlocks[i] = index_qty_glob[i] - qty_offset_prev_qties;
// // // // // //   }
// // // // // // 
// // // // // //   //element Quantity matrix
// // // // // //   double*** B;
// // // // // //   B = new double**[ncomps];
// // // // // //   for (unsigned i=0; i< ncomps; i++) {
// // // // // //     B[i] = new double*[ncomps];
// // // // // //     for (unsigned j=0; j< ncomps; j++) B[i][j] = new double[MAX_EL_NODES*MAX_EL_NODES];
// // // // // //   }
// // // // // //   double** Rhs;
// // // // // //   Rhs = new double*[ncomps];
// // // // // //   for (unsigned i=0; i< ncomps; i++)   Rhs[i] = new double[MAX_EL_NODES];
// // // // // // 
// // // // // //   unsigned SolIndex_qty      = GetIndex(qty_comps_names[component_pick_SolIndex].c_str());
// // // // // //   unsigned order_ind_qty     = Lin_Solver_[level_pick_order]->SolType[SolIndex_qty];  //the function here is INDEPENDENT of the LEVEL where I am taking it from!
// // // // // //   unsigned end_ind_var_qty   = Lin_Solver_[level_pick_order]->END_IND[order_ind_qty];         //the function here is INDEPENDENT of the LEVEL where I am taking it from!
// // // // // //  // end variables  ========================
// // // // // //   
// // // // // //   
// // // // // //      
// // // // // //   for (unsigned kel=0;kel<nel;kel++) {
// // // // // // 
// // // // // //    unsigned nve_qty = myel->GetElementDofNumber(kel,end_ind_var_qty);     //this should depend on the 
// // // // // //    
// // // // // //     short unsigned kelt = myel->GetElementType(kel);
// // // // // //     
// // // // // //     for (unsigned i=0; i<ncomps; i++) {
// // // // // //       memset(Rhs[index_qty_localBlocks[i]],0,nve_qty*sizeof(double));
// // // // // //       for (unsigned j=0; j<ncomps; j++) {
// // // // // //         memset(B[index_qty_localBlocks[i]][index_qty_localBlocks[j]],0,nve_qty*nve_qty*sizeof(double));
// // // // // //       }
// // // // // //     }
// // // // // // 
// // // // // //     for (unsigned i=0;i<nve_qty;i++) {
// // // // // //       unsigned inode = myel->GetElementDofIndex(kel,i)-1u;
// // // // // //       node_geom_el[i] = inode;
// // // // // //       for (unsigned j=0; j<ncomps; j++)       node_qty[ index_qty_localBlocks[j] ][i] = node_geom_el[i] + myKKIndex[index_qty_glob[j]];   //global
// // // // // //     }
// // // // // // 
// // // // // // 	
// // // // // //     for (unsigned i=0;i<ncomps;i++) {
// // // // // //       ierr = VecSetValues(myRES,nve_qty,node_qty[ index_qty_localBlocks[i] ],Rhs[ index_qty_localBlocks[i] ],ADD_VALUES);    CHKERRABORT(PETSC_COMM_WORLD,ierr);
// // // // // //       for (unsigned j=0;j<ncomps;j++) {
// // // // // //         ierr = MatSetValuesBlocked(myKK,nve_qty,node_qty[ index_qty_localBlocks[i] ],nve_qty,node_qty[ index_qty_localBlocks[j] ],B[ index_qty_localBlocks[i] ][ index_qty_localBlocks[j] ],ADD_VALUES);     CHKERRABORT(PETSC_COMM_WORLD,ierr);
// // // // // //         }
// // // // // //      }
// // // // // //     
// // // // // //     } //end element loop
// // // // // //    
// // // // // //    
// // // // // //   }//end quantity/equation loop
// // // // // //    
// // // // // //    
// // // // // //   //once you call the AssembleBegin and AssembleEnd functions, you cannot change the sparsity pattern any longer! 
// // // // // //   //BEGIN MATRIX ASSEMBLY ==========
// // // // // //   ierr = MatAssemblyBegin(myKK,MAT_FINAL_ASSEMBLY);     CHKERRABORT(PETSC_COMM_WORLD,ierr);
// // // // // //   ierr = MatAssemblyEnd(myKK,MAT_FINAL_ASSEMBLY);     CHKERRABORT(PETSC_COMM_WORLD,ierr);
// // // // // //   //END MATRIX ASSEMBLY ============
// // // // // // 
// // // // // //   //BEGIN RESIDUAL ASSEMBLY ==========
// // // // // //   ierr = VecAssemblyBegin(myRES);     CHKERRABORT(PETSC_COMM_WORLD,ierr);
// // // // // //   ierr = VecAssemblyEnd(myRES);     CHKERRABORT(PETSC_COMM_WORLD,ierr);
// // // // // //   //END RESIDUAL ASSEMBLY ============
// // // // // //   
// // // // // //   
// // // // // // //         PetscViewer viewer;
// // // // // // //         ierr=PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);CHKERRABORT(PETSC_COMM_WORLD,ierr);
// // // // // // //         ierr= MatView(myKK,viewer);CHKERRABORT(PETSC_COMM_WORLD,ierr);
// // // // // // // //         ierr= VecView(myRES,viewer);CHKERRQ(ierr);
// // // // // // //         double ff;
// // // // // // //         std::cin>>ff;
// // // // // // 
// // // // // // 
// // // // // //     } //end level  //----------------------------------------------------------------------------------------------
// // // 
// // //   return;
// // // 
// // //   }
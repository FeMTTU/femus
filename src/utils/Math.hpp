#ifndef __femus_utils_Math_hpp__
#define __femus_utils_Math_hpp__


#include "Typedefs.hpp"

#include "CurrentQuantity.hpp"
#include "ElemType.hpp"
#include "MultiLevelProblem.hpp"
#include "MultiLevelSolution.hpp"

namespace femus {

// remember that you have to declare all these functions "inline", otherwise you get "multiple definition" in linking

 // Operations ---------------------------------

namespace Math {
  
   inline void zeroN(double* x,const uint N); 
   inline double dotN(const double* x,const double* y,const uint N);  //TODO this will be deleted 
   inline double dot(const double* x,const double* y, const uint spacedim);
   inline  void cross(const double* a,const double* b, double* res);
   inline void extend(const double* a, double* a3D, const uint spacedim);
   inline void extend_nds(const uint,const double*, double*, const uint spacedim);
   inline void normalize(double* x,const double fac, const uint spacedim);



// ============== inline functions ==========================   
   
/// set to zero - n components
inline void zeroN(double* x,const uint N)  {

  for (uint i=0; i< N; i++)  x[i]=0.;
  return;
}


/// dot product
inline void normalize(double* x,const double fac, const uint spacedim)  {
  for (uint idim=0; idim< spacedim; idim++)  x[idim] /= fac;
  return;
}

/// dot product - n components
inline double dotN(const double* x,const double* y,const uint N)  {
  double dotprod=0.;
  for (uint idim=0; idim< N; idim++)  dotprod += x[idim]*y[idim];
  return dotprod;
}

/// dot product
inline double dot(const double* x,const double* y, const uint spacedim) {
  double dotprod=0.;
  for (uint idim=0; idim < spacedim; idim++)  dotprod += x[idim]*y[idim];
  return dotprod;
}

/// Cross product
inline void cross(const double* a,const double* b, double* res) {
//a,b,res are 3D vectors
//clean then fill
  for (uint i=0; i<3; i++) res[i]=0.;
  for (uint i=0; i<3; i++) res[i] = (a[(i+1)%3]*b[(i+2)%3] - a[(i+2)%3]*b[(i+1)%3]);
  return;
}

/// extend to a 3D vector a vector with dimension 
inline void extend(const double* a, double* a3D, const uint spacedim)  {
  for (uint i=0; i<3; i++) a3D[i]=0.;
  for (uint i=0; i< spacedim; i++)  a3D[i]=a[i];
  return;
}

/// extend to 3D an element dof vector
inline void extend_nds(const uint el_ndofs,const double* a_nds, double* a_nds3D, const uint spacedim)  {

//AAA: valid from spacedim to 3

//set to zero
  for (uint eln=0; eln<el_ndofs; eln++)  {
    for (uint i=0; i<3; i++) {
      a_nds3D[eln+i*el_ndofs]=0.;
    }
  }
//extend
  for (uint eln=0; eln<el_ndofs; eln++)    {
    for (uint idim=0; idim<spacedim; idim++) {
      a_nds3D[eln+idim*el_ndofs] = a_nds[eln+idim*el_ndofs];
    }
  }

  return;
}


//=================================================================

inline double FunctionIntegral (const uint vb, MultiLevelProblem & ml_prob, double (*pt2func)(double, const std::vector<double> ) )  {

  const uint mesh_vb = vb;
  
  
const uint Level  = ml_prob.GetMeshTwo()._NoLevels - 1;
const uint myproc = ml_prob.GetMeshTwo()._iproc;
  double time = 0.;
  
  Mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(Level);
  

  double integral = 0.;
  

//parallel sum
    const uint nel_e = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level+1];
    const uint nel_b = ml_prob.GetMeshTwo()._off_el[mesh_vb][ml_prob.GetMeshTwo()._NoLevels*myproc+Level];
  
    for (uint iel=0; iel < (nel_e - nel_b); iel++) {
      
    CurrentElem       currelem(iel,myproc,Level,vb,NULL,ml_prob.GetMeshTwo(),ml_prob.GetElemType(),mymsh); //element without equation
    CurrentGaussPointBase & currgp = CurrentGaussPointBase::build(currelem,ml_prob.GetQrule(currelem.GetDim()));

    const uint el_ngauss = ml_prob.GetQrule(currelem.GetDim()).GetGaussPointsNumber();

//========= DOMAIN MAPPING
    CurrentQuantity xyz(currgp);
    xyz._dim      = ml_prob.GetMeshTwo().get_dim();
    xyz._FEord    = MESH_MAPPING_FE;
    xyz._ndof     = currelem.GetElemType(xyz._FEord)->GetNDofs();
    xyz.Allocate();

    currelem.SetDofobjConnCoords();
    
    currelem.ConvertElemCoordsToMappingOrd(xyz);

     
    for (uint qp = 0; qp < el_ngauss; qp++) {

       for (uint fe = 0; fe < QL; fe++)     {  
             currgp.SetPhiElDofsFEVB_g (fe,qp); 
             currgp.SetDPhiDxezetaElDofsFEVB_g (fe,qp);
        }  
     
double  Jac_g=0.;
          if (vb==0)   Jac_g = currgp.JacVectVV_g(xyz);  //not xyz_refbox!      
     else if (vb==1)   Jac_g = currgp.JacVectBB_g(xyz);  //not xyz_refbox!      

   const double  wgt_g = ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);


 xyz.val_g();
double myval_g = pt2func(time,xyz._val_g); 

 
  integral += wgt_g*Jac_g*myval_g;

   
    }//gauss loop
     
    }//element loop
    
         std::cout << std::endl  << " ^^^^^^^^^^^^^^^^^L'integrale sul processore "<< myproc << " vale: " << integral << std::endl;

//     double weights_sum = 0.;
//     for (uint qp = 0; qp < el_ngauss; qp++)  weights_sum += ml_prob.GetQrule(currelem.GetDim()).GetGaussWeight(qp);
//        std::cout << std::endl << " ^^^^^^^^^^^^^^^^^ La somma dei pesi  vale: " << weights_sum << std::endl;

       double J=0.;
   #ifdef HAVE_MPI
//       MPI_Reduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );  //This one gives J only to processor 0 !
      MPI_Allreduce( &integral, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );  //THIS IS THE RIGHT ONE!!
   #else
   J = integral;
   #endif
    
     std::cout << std::endl << " ^^^^^^^^^^^^^^^^^L'integrale totale vale: " << J << std::endl;

  return J;  
  
}



} //end namespace Math



namespace FE_convergence {
 
 
template < class type >
  class Function {
 
  public:
      
 virtual type value(const std::vector < type >& x) const = 0;

 virtual vector < type >  gradient(const std::vector < type >& x) const  = 0;

 virtual type laplacian(const std::vector < type >& x) const = 0;
 
  type helmholtz(const std::vector < type >& x) const { return ( - laplacian(x) + value(x) ); };

};
  


//this is based on the AddSolution function in MLSol
 class Unknowns_definition {
     
 public:
     
     std::string _name;
     FEFamily _fe_family;     
     FEOrder  _fe_order;     
     
 };
    
 
template < class type >
 inline   std::vector < std::vector < std::vector < type > > >  initialize_vector_of_norms(const unsigned unknowns_size, 
                                                                                             const unsigned max_number_of_meshes, 
                                                                                             const unsigned norm_flag) {
   
       //how many Unknowns, how many mesh levels, how many norms
       
   std::vector < std::vector < std::vector < type > > > norms( unknowns_size );
  
     for (unsigned int u = 0; u < unknowns_size; u++) {
              norms[u].resize( max_number_of_meshes );
       for (int i = 0; i < norms[u].size(); i++) {   // loop on the mesh level
               norms[u][i].resize(norm_flag + 1);
           }   
       }
 //Error norm definition  ==================
    
    return norms; 
       
}

    
   
  inline const MultiLevelSolution  initialize_convergence_study(const std::vector< FE_convergence::Unknowns_definition > &  unknowns,  
                                                                MultiLevelMesh & ml_mesh_all_levels, 
                                                                const unsigned max_number_of_meshes, 
                                                                const MultiLevelSolution::BoundaryFunc SetBoundaryCondition)  {

   //Mesh: construct all levels  ==================
        unsigned numberOfUniformLevels_finest = max_number_of_meshes;
        ml_mesh_all_levels.RefineMesh(numberOfUniformLevels_finest, numberOfUniformLevels_finest, NULL);
//      ml_mesh_all_levels.EraseCoarseLevels(numberOfUniformLevels - 2);  // need to keep at least two levels to send u_(i-1) projected(prolongated) into next refinement

 
 //Solution ==================
//         std::vector < MultiLevelSolution * >   ml_sol_all_levels(unknowns.size());
//                ml_sol_all_levels[u] = new MultiLevelSolution (& ml_mesh_all_levels);  //with the declaration outside and a "new" inside it persists outside the loop scopes
               MultiLevelSolution ml_sol_all_levels(& ml_mesh_all_levels);
               
     for (unsigned int u = 0; u < unknowns.size(); u++) {
               ml_sol_all_levels.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order);  //We have to do so to avoid buildup of AddSolution with different FE families
               ml_sol_all_levels.Initialize("All");
               ml_sol_all_levels.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
               ml_sol_all_levels.GenerateBdc("All");
                }
 //Solution  ==================

 return ml_sol_all_levels;
} 
    


//   print the error and the order of convergence between different levels
 template <class type>
inline void output_convergence_order(const std::vector < std::vector < std::vector < type > > > &  norm,
                                    const unsigned int u,
                                    const unsigned int i,
                                    const unsigned int n) {

   if(i < norm[u].size() - 2)  {
//   std::cout << norm_name << " ERROR and ORDER OF CONVERGENCE: " << fam << " " << ord << "\n\n";

    std::cout << i + 1 << "\t\t";
    std::cout.precision(14);

    std::cout << norm[u][i][n] << "\t";

    std::cout << std::endl;
  
      std::cout.precision(3);
      std::cout << "\t\t";

        std::cout << log( norm[u][i][n] / norm[u][i + 1][n] ) / log(2.) << "  \t\t\t\t";

      std::cout << std::endl;
    }
    
  
  
}


 template <class type>
inline void output_convergence_order_all(const std::vector< FE_convergence::Unknowns_definition > &  unknowns,
                                        const std::vector < std::vector < std::vector < type > > > &  norms, 
                                        const unsigned norm_flag, 
                                        const unsigned max_number_of_meshes) {
    
    assert( unknowns.size() == norms.size() );
    
    const std::vector< std::string > norm_names = {"L2-NORM","H1-SEMINORM"};
  
     for (unsigned int u = 0; u < unknowns.size(); u++) {
       for (int n = 0; n < norm_flag + 1; n++) {
            std::cout << unknowns[u]._name << " : " << norm_names[n] << " ERROR and ORDER OF CONVERGENCE"  << std::endl;
         for (int i = 0; i < max_number_of_meshes; i++) {
                output_convergence_order(norms,u,i,n);
            }
            std::cout << std::endl;
         }
      }
      
 }
     


 
 

 template <class type>
inline std::vector< type > compute_error_norms(const MultiLevelSolution* ml_sol, 
                                              const MultiLevelSolution* ml_sol_all_levels,
                                              const std::string & unknown,
                                              const unsigned current_level,
                                              const unsigned norm_flag,
                                              const unsigned conv_order_flag,
                                              const Function< type > * ex_sol_in = NULL
                                             ) {
     
  // (//0 = only L2: //1 = L2 + H1)
  
// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv 
//i.e. ||prol_(u_(i-1)) - u_(i)|| = err(i) => err(i-1)/err(i) = 2^alpha ,implemented as log(err(i)/err(i+1))/log2

   if (conv_order_flag == 1 && ex_sol_in == NULL) { std::cout << "Please provide analytical solution" << std::endl; abort(); }
   
    
  const unsigned num_norms = norm_flag + 1;
  //norms that we are computing here //first L2, then H1 ============
  std::vector< type > norms(num_norms);                  std::fill(norms.begin(), norms.end(), 0.);   
  std::vector< type > norms_exact_dofs(num_norms);       std::fill(norms_exact_dofs.begin(), norms_exact_dofs.end(), 0.);
  std::vector< type > norms_inexact_dofs(num_norms);     std::fill(norms_inexact_dofs.begin(), norms_inexact_dofs.end(), 0.);
  //norms that we are computing here //first L2, then H1 ============
  
  
  unsigned level = ml_sol->_mlMesh->GetNumberOfLevels() - 1u;
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = ml_sol->_mlMesh->GetLevel(level);
  elem*     el  = msh->el;
  const Solution* sol = ml_sol->GetSolutionLevel(level);

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

 unsigned iproc = msh->processor_id();

  //solution variable
  unsigned soluIndex = ml_sol->GetIndex(unknown.c_str()); // ml_sol->GetIndex("u");    // get the position of "u" in the ml_sol object
  unsigned soluType = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  
  // ======================================
  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  
  
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)
  vector < vector < type > > x(dim);    // local coordinates
  for (unsigned i = 0; i < dim; i++)   x[i].reserve(maxSize);

//-----------------  
  vector < type > phi_coords;
  vector < type > phi_coords_x;
  vector < type > phi_coords_xx;

  phi_coords.reserve(maxSize);
  phi_coords_x.reserve(maxSize * dim);
  phi_coords_xx.reserve(maxSize * dim2);

  vector < type > phi;
  vector < type > phi_x;
  vector < type > phi_xx;
  
  type weight; // gauss point weight


  vector < type >  solu;                               solu.reserve(maxSize);
  vector < type >  solu_exact_at_dofs;   solu_exact_at_dofs.reserve(maxSize);
  vector < type >  solu_coarser_prol;     solu_coarser_prol.reserve(maxSize);


  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);
  phi_xx.reserve(maxSize * dim2);


  
  // element loop: each process loops only on the elements that owns
  for (int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {

    
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);

    // resize local arrays
    solu.resize(nDofu);
    solu_exact_at_dofs.resize(nDofu);
    solu_coarser_prol.resize(nDofu);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->_topology->_Sol[jdim])(xDof);  // global extraction and local storage for the element coordinates
      }
    }
    
    

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
        std::vector< type > x_at_node(dim,0.);
        for (unsigned jdim = 0; jdim < dim; jdim++) x_at_node[jdim] = x[jdim][i];
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);
                   solu[i]  =                                        (*sol->_Sol[soluIndex])(solDof);
      solu_coarser_prol[i]  = (*ml_sol_all_levels->GetSolutionLevel(current_level)->_Sol[soluIndex])(solDof);
      if (ex_sol_in != NULL) solu_exact_at_dofs[i] = ex_sol_in->value(x_at_node);
    }


    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
        
      // *** get gauss point weight, test function and test function partial derivatives ***
     static_cast<const elem_type_2D*>( msh->_finiteElement[ielGeom][soluType] )
                                         ->Jacobian_type_non_isoparametric< type >( static_cast<const elem_type_2D*>( msh->_finiteElement[ielGeom][xType] ), x, ig, weight, phi, phi_x, phi_xx);
//       msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, phi_xx);
      msh->_finiteElement[ielGeom][xType]->Jacobian(x, ig, weight, phi_coords, phi_coords_x, phi_coords_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      type solu_gss = 0.;
      type exactSol_from_dofs_gss = 0.;
      type solu_coarser_prol_gss = 0.;
      vector < type > gradSolu_gss(dim, 0.);
      vector < type > gradSolu_exact_at_dofs_gss(dim, 0.);
      vector < type > gradSolu_coarser_prol_gss(dim, 0.);
      vector < type > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss                += phi[i] * solu[i];
        exactSol_from_dofs_gss  += phi[i] * solu_exact_at_dofs[i];
        solu_coarser_prol_gss   += phi[i] * solu_coarser_prol[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim]                += phi_x[i * dim + jdim] * solu[i];
          gradSolu_exact_at_dofs_gss[jdim]  += phi_x[i * dim + jdim] * solu_exact_at_dofs[i];
          gradSolu_coarser_prol_gss[jdim]   += phi_x[i * dim + jdim] * solu_coarser_prol[i];
          x_gss[jdim] += x[jdim][i] * phi_coords[i];
        }
      }

// H^0 ==============      
//     if (norm_flag == 0) {
      type exactSol = 0.; if (ex_sol_in != NULL) exactSol = ex_sol_in->value(x_gss);
      norms[0]               += (solu_gss - exactSol)                * (solu_gss - exactSol)       * weight;
      norms_exact_dofs[0]    += (solu_gss - exactSol_from_dofs_gss)  * (solu_gss - exactSol_from_dofs_gss) * weight;
      norms_inexact_dofs[0]  += (solu_gss - solu_coarser_prol_gss)   * (solu_gss - solu_coarser_prol_gss)  * weight;
//     }
    
// H^1 ==============      
    /*else*/ if (norm_flag == 1) {
      vector < type > exactGradSol(dim,0.);    if (ex_sol_in != NULL) exactGradSol = ex_sol_in->gradient(x_gss);

      for (unsigned j = 0; j < dim ; j++) {
        norms[1]               += ((gradSolu_gss[j] - exactGradSol[j])               * (gradSolu_gss[j]  - exactGradSol[j])) * weight;
        norms_exact_dofs[1]    += ((gradSolu_gss[j] - gradSolu_exact_at_dofs_gss[j]) * (gradSolu_gss[j] - gradSolu_exact_at_dofs_gss[j])) * weight;
        norms_inexact_dofs[1]  += ((gradSolu_gss[j] - gradSolu_coarser_prol_gss[j])  * (gradSolu_gss[j] - gradSolu_coarser_prol_gss[j]))  * weight;
      }
   }
      
   } // end gauss point loop
   
  } //end element loop for each process

  // add the norms of all processes
  NumericVector* norm_vec;
                 norm_vec = NumericVector::build().release();
                 norm_vec->init(msh->n_processors(), 1 , false, AUTOMATIC);

         /*if (norm_flag == 0) {*/ norm_vec->set(iproc, norms[0]/*.value()*/);  norm_vec->close();  norms[0] = norm_vec->l1_norm(); /*}*/
    /*else*/ if (norm_flag == 1) { norm_vec->set(iproc, norms[1]/*.value()*/);  norm_vec->close();  norms[1] = norm_vec->l1_norm(); }

          delete norm_vec;

   // add the norms of all processes
  NumericVector* norm_vec_exact_dofs;
                 norm_vec_exact_dofs = NumericVector::build().release();
                 norm_vec_exact_dofs->init(msh->n_processors(), 1 , false, AUTOMATIC);

         /*if (norm_flag == 0) {*/ norm_vec_exact_dofs->set(iproc, norms_exact_dofs[0]/*.value()*/);  norm_vec_exact_dofs->close();  norms_exact_dofs[0] = norm_vec_exact_dofs->l1_norm(); /*}*/
    /*else*/ if (norm_flag == 1) { norm_vec_exact_dofs->set(iproc, norms_exact_dofs[1]/*.value()*/);  norm_vec_exact_dofs->close();  norms_exact_dofs[1] = norm_vec_exact_dofs->l1_norm(); }

          delete norm_vec_exact_dofs;

  // add the norms of all processes
  NumericVector* norm_vec_inexact;
                 norm_vec_inexact = NumericVector::build().release();
                 norm_vec_inexact->init(msh->n_processors(), 1 , false, AUTOMATIC);

         /*if (norm_flag == 0) {*/ norm_vec_inexact->set(iproc, norms_inexact_dofs[0]/*.value()*/);  norm_vec_inexact->close();  norms_inexact_dofs[0] = norm_vec_inexact->l1_norm(); /*}*/
    /*else*/ if (norm_flag == 1) { norm_vec_inexact->set(iproc, norms_inexact_dofs[1]/*.value()*/);  norm_vec_inexact->close();  norms_inexact_dofs[1] = norm_vec_inexact->l1_norm(); }

          delete norm_vec_inexact;

          
    for (int n = 0; n < norms.size(); n++) { 
                  norms[n] = sqrt(norms[n]);                                            
       norms_exact_dofs[n] = sqrt(norms_exact_dofs[n]);                                            
     norms_inexact_dofs[n] = sqrt(norms_inexact_dofs[n]);
    }
    
    
if (conv_order_flag == 0)  return norms_inexact_dofs;
if (conv_order_flag == 1)  return norms;
 //   return norms_exact_dofs;

 
} 
 
 
 

     
template <class type>
 inline void compute_error_norms_per_unknown_per_level(const MultiLevelSolution* ml_sol_single_level, 
                                          MultiLevelSolution* ml_sol_all_levels, 
                                          const std::vector< FE_convergence::Unknowns_definition > &  unknowns, 
                                          const unsigned i,
                                          const unsigned norm_flag, 
                                          std::vector < std::vector < std::vector < type > > > &  norms,
                                          const unsigned conv_order_flag,
                                          const Function< type > * ex_sol_in = NULL
                                         ) {
     
     
        if ( i > 0 ) {

            // ======= prolongate to the current level (i) from the coarser level (i-1) (so that you can compare the two) ========================
            ml_sol_all_levels->RefineSolution(i);
            
            // =======  compute the error norm at the current level (i) ========================
            for (unsigned int u = 0; u < unknowns.size(); u++) {  //this loop could be inside the below function
                
            const std::vector< type > norm_out = FE_convergence::compute_error_norms< type >(ml_sol_single_level, ml_sol_all_levels, unknowns[u]._name, i, norm_flag, conv_order_flag, ex_sol_in);

              for (int n = 0; n < norms[u][i-1].size(); n++)      norms[u][i-1][n] = norm_out[n];
                                       
                   }
                   
        }
                
                 
              // ======= store the last computed solution to prepare the next iteration (the current level i is now overwritten) ========================
            ml_sol_all_levels->fill_at_level_from_level(i, ml_sol_single_level->_mlMesh->GetNumberOfLevels() - 1, *ml_sol_single_level);
        
                 
                 
                 
 } 
 


 
    
} //end FE_convergence




} //end namespace femus



#endif

#ifndef __femus_utils_FE_convergence_hpp__
#define __femus_utils_FE_convergence_hpp__

#include "solution_generation_single_level.hpp"


#include "MultiLevelSolution.hpp"
#include "VectorOperations.hpp"
#include "Function.hpp"
#include "Assemble_unknown.hpp"

#include <vector>
#include <iostream>


namespace femus {

class MultiLevelProblem;
class MultiLevelMesh;
class MultiLevelSolution;





/// The idea of this class is to generate a certain amount of MultilevelSolutions and to compute the convergence rate.
/// It may be used either to study Approximation Theory or to study Convergence of FE Solution of PDEs
template < class real_num = double >
class FE_convergence {
 
    
    
public: 
 

static  void  convergence_study(MultiLevelProblem & ml_prob,
                          MultiLevelMesh & ml_mesh,
                          MultiLevelMesh * ml_mesh_all_levels_needed_for_incremental,   //auxiliary
                          const unsigned max_number_of_meshes,
                          const std::vector < bool > convergence_rate_computation_method_Flag,
                          const std::vector < bool > volume_or_boundary_Flag,
                          const std::vector < bool > sobolev_norms_Flag,
                          const bool                               main_in_has_equation_solve,                                          // true only if I have a System inside
                          const Solution_generation_single_level & main_in,
                          const std::vector< Unknown > & unknowns,                            //vector of Solutions
                          const std::vector< Math::Function< double > * >  & exact_sol,       // analytical function to Approximate
                          const MultiLevelSolution::InitFuncMLProb SetInitialCondition,       // routine that does the Approximation - needed in all cases
                          const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition  // routine for the BC - needed only if I have a System inside
                         );    //vector of Exact Solutions, if available

    
//   print the error and the order of convergence between different levels
static  void output_convergence_order(const std::vector < std::vector < real_num > >  &  norm,
                                      const unsigned int u,
                                      const unsigned int i);


private: 
  
static       const std::vector< std::string > norm_names;
static       const std::vector< std::string > convergence_computation_names;
static       const std::vector< std::string > volume_or_boundary_names;


static   const MultiLevelSolution  initialize_convergence_study(const MultiLevelProblem & ml_prob,
                                                                const std::vector< Unknown > &  unknowns,  
                                                                const std::vector< Math::Function< double > * >  & exact_sol,
                                                                MultiLevelMesh & ml_mesh_all_levels_needed_for_incremental, 
                                                                const unsigned max_number_of_meshes, 
                                                                const MultiLevelSolution::InitFuncMLProb SetInitialCondition_in,
                                                                const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition,
                                                                const bool main_in_has_equation_solve
                                                               );


static   std::vector < std::vector < real_num > >    initialize_vector_of_norms(const unsigned unknowns_size, 
                                                                                const unsigned max_number_of_meshes);

    
   



static  void output_convergence_order_all(const std::vector< Unknown > &  unknowns,
                                        const std::vector < std::vector < real_num > >   &  norms, 
                                        const unsigned sobolev_norms, 
                                        const unsigned volume_or_boundary,
                                        const unsigned convergence_computation_method);


static  void output_convergence_rate( double norm_i, double norm_ip1, std::string norm_name, unsigned maxNumberOfMeshes , int loop_i); 


static  std::vector< real_num > compute_error_norms_volume_or_boundary_with_analytical_sol_or_not(const std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, double> *  > > & elem_all,
                                                const std::vector<Gauss> & quad_rules,
                                                const MultiLevelSolution* ml_sol, 
                                                const MultiLevelSolution* ml_sol_all_levels_needed_for_incremental,
                                                const std::string & unknown,
                                                const Math::Function< real_num > * ex_sol_in,
                                                const unsigned current_level,
                                                const unsigned sobolev_norms,
                                                const unsigned convergence_rate_computation_method,
                                                const unsigned volume_or_boundary
                                             );


static  void compute_error_norms_per_unknown_per_level(const std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, double> *  > > & elem_all,
                                                       const std::vector<Gauss> & quad_rules,
                                                       const MultiLevelSolution* ml_sol_single_level,
                                                       MultiLevelSolution* ml_sol_all_levels_needed_for_incremental,
                                                       const std::vector< Unknown > &  unknowns,
                                                       const std::vector< Math::Function< real_num > * > & ex_sol_in,
                                                       std::vector < std::vector < real_num > >   &  norms,
                                                       const unsigned i,
                                                       const unsigned sobolev_norms,
                                                       const unsigned volume_or_boundary,
                                                       const unsigned convergence_rate_computation_method
                                                       );


 
    
};



template < class real_num>
/*static*/       const std::vector< std::string > FE_convergence< real_num >::norm_names = {"L2-NORM", "H1-SEMINORM"};
template < class real_num>
/*static*/       const std::vector< std::string > FE_convergence< real_num >::convergence_computation_names = {"Incremental", "Absolute, with analytical sol"};
template < class real_num>
/*static*/       const std::vector< std::string > FE_convergence< real_num >::volume_or_boundary_names = {"Volume", "Boundary"};



    

template < class real_num>
/*static*/  void  FE_convergence< real_num >::convergence_study(MultiLevelProblem & ml_prob,
                                                  MultiLevelMesh & ml_mesh,
                                                  MultiLevelMesh * ml_mesh_all_levels_needed_for_incremental,
                                                  const unsigned max_number_of_meshes,
                                                  const std::vector < bool > convergence_rate_computation_method_Flag,
                                                  const std::vector < bool > volume_or_boundary_Flag,
                                                  const std::vector < bool > sobolev_norms_Flag,
                                                  const bool main_in_has_equation_solve,
                                                  const Solution_generation_single_level & main_in,
                                                  const std::vector< Unknown > & unknowns,
                                                  const std::vector< Math::Function< double > * >  & exact_sol,
                                                  const MultiLevelSolution::InitFuncMLProb SetInitialCondition,
                                                  const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition
                                                 ) {
   
      // Controls - BEGIN 
      if ( main_in_has_equation_solve == true &&  SetBoundaryCondition == NULL)  { std::cout << "Must provide a BC function" << std::endl; abort(); }
      
       if (convergence_rate_computation_method_Flag[0] == true &&  ml_mesh_all_levels_needed_for_incremental == NULL)  {
         std::cout << "Incremental method requires auxiliary mesh for all levels" << std::endl; abort(); }
         
       if (convergence_rate_computation_method_Flag[1] == true && exact_sol.size() == 0)  { std::cout << "Must provide exact sol" << std::endl; abort(); }

    
    if ( convergence_rate_computation_method_Flag.size() !=   convergence_computation_names.size() ) { std::cout << "Not implemented" << std::endl; abort(); }
    if ( volume_or_boundary_Flag.size() != volume_or_boundary_names.size() )                  { std::cout << "Not implemented" << std::endl; abort(); }
    if ( sobolev_norms_Flag.size() != norm_names.size() )                       { std::cout << "Not implemented" << std::endl; abort(); }
       // Controls - END
  
    
    // ================= initialize norms - BEGIN 
    
     std::vector < std::vector < std::vector < std::vector < std::vector < double > > > > >   norms;
     
     norms.resize(convergence_rate_computation_method_Flag.size());
     
      for (unsigned convergence_rate_computation_method = 0; 
                    convergence_rate_computation_method < convergence_rate_computation_method_Flag.size(); 
                    convergence_rate_computation_method++) {
        //0: incremental 1: absolute (with analytical sol)  2: absolute (with projection of finest sol)...
        
     norms[convergence_rate_computation_method].resize( volume_or_boundary_Flag.size() );  

      for (unsigned volume_or_boundary = 0;
                    volume_or_boundary < volume_or_boundary_Flag.size(); 
                    volume_or_boundary++) {// (//0 = only V: //1 = only B) 

        norms[convergence_rate_computation_method][volume_or_boundary].resize( sobolev_norms_Flag.size() );   
     
      for (unsigned sobolev_norms = 0; 
                    sobolev_norms < sobolev_norms_Flag.size();
                    sobolev_norms++) {  //  what Sobolev norms to compute //0 = only L2: //1 = only H1)

     // norms[convergence_rate_computation_method][volume_or_boundary][sobolev_norms].resize( unknowns.size() );

     norms[convergence_rate_computation_method][volume_or_boundary][sobolev_norms] = FE_convergence::initialize_vector_of_norms (unknowns.size(), max_number_of_meshes);
     
           }
        }
      }
     
    // ================= initialize norms - END
     
    
    // ================= initialize convergence study, ONLY FOR INCREMENTAL - BEGIN 
     MultiLevelSolution         ml_sol_all_levels_needed_for_incremental = FE_convergence::initialize_convergence_study(ml_prob,
                                                                                                 unknowns,
                                                                                                 exact_sol,
                                                                                                 *ml_mesh_all_levels_needed_for_incremental, 
                                                                                                 max_number_of_meshes, 
                                                                                                 SetInitialCondition, 
                                                                                                 SetBoundaryCondition,
                                                                                                 main_in_has_equation_solve);
    // ================= initialize convergence study, ONLY FOR INCREMENTAL - END


    // ================= Error computation at all levels - BEGIN =============

     
    // ================= FE Evaluations at Quadrature - BEGIN 
  //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, double> *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
    // ================= FE Evaluations at Quadrature - END
  
  
            
       for (int lev = 0; lev < max_number_of_meshes; lev++) {
         
                  
    // ================= Compute solution - BEGIN 
            const MultiLevelSolution ml_sol_single_level = main_in.run_on_single_level(ml_prob,
                                                                                       ml_mesh,
                                                                                       lev,
                                                                                       unknowns,
                                                                                       exact_sol,
                                                                                       SetInitialCondition, 
                                                                                       SetBoundaryCondition,
                                                                                       main_in_has_equation_solve
                                                                                      );
    // ================= Compute solution - END
            
            
    // ================= Compute error - BEGIN 
      for (unsigned convergence_rate_computation_method = 0; 
                    convergence_rate_computation_method < convergence_rate_computation_method_Flag.size();
                    convergence_rate_computation_method++) {
        if( convergence_rate_computation_method_Flag[convergence_rate_computation_method] ) {
        
        for (unsigned volume_or_boundary = 0;
                      volume_or_boundary < volume_or_boundary_Flag.size();
                      volume_or_boundary++) { 
          if( volume_or_boundary_Flag[volume_or_boundary] ) {

     
           for (unsigned sobolev_norms = 0; 
                         sobolev_norms < sobolev_norms_Flag.size();
                         sobolev_norms++) {
          if( sobolev_norms_Flag[sobolev_norms] ) {
        

            FE_convergence::compute_error_norms_per_unknown_per_level ( elem_all,
                                                                        ml_prob.GetQuadratureRuleAllGeomElems(),
                                                                        & ml_sol_single_level,
                                                                        & ml_sol_all_levels_needed_for_incremental,
                                                                        unknowns,
                                                                        exact_sol,
                                                                        norms[convergence_rate_computation_method][volume_or_boundary][sobolev_norms],
                                                                        lev,
                                                                        sobolev_norms,
                                                                        volume_or_boundary,
                                                                        convergence_rate_computation_method);
           }}
        }}
       }}
    // ================= Compute error - END 
        
        
      }
      
    // ================= Error computation at all levels - END =============
      
   
    // ================= Error output - BEGIN 
      for (unsigned convergence_rate_computation_method = 0; 
                    convergence_rate_computation_method < convergence_rate_computation_method_Flag.size();
                    convergence_rate_computation_method++) {
        if( convergence_rate_computation_method_Flag[convergence_rate_computation_method] ) {
        
        for (unsigned volume_or_boundary = 0;
                      volume_or_boundary < volume_or_boundary_Flag.size();
                      volume_or_boundary++) { 
          if( volume_or_boundary_Flag[volume_or_boundary] ) {

     
           for (unsigned sobolev_norms = 0; 
                         sobolev_norms < sobolev_norms_Flag.size();
                         sobolev_norms++) {
          if( sobolev_norms_Flag[sobolev_norms] ) {
                     
       FE_convergence::output_convergence_order_all(unknowns,
                                                    norms[convergence_rate_computation_method][volume_or_boundary][sobolev_norms],
                                                    sobolev_norms,
                                                    volume_or_boundary,
                                                    convergence_rate_computation_method);
       
             }}
          }}
        }}
    // ================= Error output - END 
 
 
}

    
    
   
template < class real_num>
/*static*/   const MultiLevelSolution  FE_convergence< real_num >::initialize_convergence_study(const MultiLevelProblem & ml_prob, 
                                                                                            const std::vector< Unknown > &  unknowns,
                                                                                            const std::vector< Math::Function< double > * > &  exact_sol,
                                                                                            MultiLevelMesh & ml_mesh_all_levels_needed_for_incremental,
                                                                                            const unsigned max_number_of_meshes,
                                                                                            const MultiLevelSolution::InitFuncMLProb SetInitialCondition_in,
                                                                                            const MultiLevelSolution::BoundaryFuncMLProb SetBoundaryCondition_in,
                                                                                            const bool main_in_has_equation_solve
                                                                                           )  {

 //Problem - BEGIN ==================
  MultiLevelProblem  ml_prob_aux(ml_prob);
 //Problem - END ==================
  
  
 //Mesh: construct all levels - BEGIN  ==================
        unsigned numberOfUniformLevels_finest = max_number_of_meshes;
        ml_mesh_all_levels_needed_for_incremental.RefineMesh(numberOfUniformLevels_finest, numberOfUniformLevels_finest, NULL);
//      ml_mesh_all_levels_needed_for_incremental.EraseCoarseLevels(numberOfUniformLevels - 2);  // need to keep at least two levels to send u_(i-1) projected(prolongated) into next refinement
 //Mesh: construct all levels - END ==================

 
 //Solution - BEGIN ==================
// MultiLevelSolution *   ml_sol_all_levels_needed_for_incremental = new MultiLevelSolution (& ml_mesh_all_levels_needed_for_incremental);  
          //with the declaration outside and a "new" inside it persists outside the loop scopes
               MultiLevelSolution ml_sol_all_levels_needed_for_incremental(& ml_mesh_all_levels_needed_for_incremental);
 //Solution - END ==================

 //Problem - BEGIN ==================
  ml_prob_aux.SetMultiLevelMeshAndSolution(& ml_sol_all_levels_needed_for_incremental);
 //Problem - END ==================
               
               
  //Solution, Initialize - BEGIN ==================
               for (unsigned int u = 0; u < unknowns.size(); u++) {
               ml_sol_all_levels_needed_for_incremental.AddSolution(unknowns[u]._name.c_str(), unknowns[u]._fe_family, unknowns[u]._fe_order);  //We have to do so to avoid buildup of AddSolution with different FE families
               ml_sol_all_levels_needed_for_incremental.set_analytical_function(unknowns[u]._name.c_str(), exact_sol[u]);   
               ml_sol_all_levels_needed_for_incremental.Initialize(unknowns[u]._name.c_str(), SetInitialCondition_in, & ml_prob_aux);
               }
  //Solution, Initialize - END ==================
               
               
 // Only for an EQUATION - BEGIN  ==================
         if (main_in_has_equation_solve) {
             
         ml_sol_all_levels_needed_for_incremental.AttachSetBoundaryConditionFunction(SetBoundaryCondition_in);
               for (unsigned int u = 0; u < unknowns.size(); u++) {
               ml_sol_all_levels_needed_for_incremental.GenerateBdc(unknowns[u]._name.c_str(), "Steady", & ml_prob_aux);
               }
               
         }
 // Only for an EQUATION - END  ==================

 return ml_sol_all_levels_needed_for_incremental;
} 
    


template < class real_num>
/*static*/    std::vector < std::vector < real_num > >  
         FE_convergence< real_num >::initialize_vector_of_norms(const unsigned unknowns_size,
                                                                const unsigned max_number_of_meshes) {
   
           
   std::vector < std::vector < real_num > >   norms( unknowns_size );
   
  
      for (unsigned int u = 0; u < norms.size(); u++) {
              norms[u].resize( max_number_of_meshes );
               std::fill(norms[u].begin(), norms[u].end(), 0.);
       }

    
    return norms; 
       
}

    
    
    

//   print the error and the order of convergence between different levels
template < class real_num>
/*static*/  void FE_convergence< real_num >::output_convergence_order(const std::vector < std::vector < real_num > >  &  norm,
                                                                      const unsigned int solution,
                                                                      const unsigned int mesh_lev
                                                                      ) {

   if(mesh_lev < norm[solution].size() - 2)  {

    std::cout << mesh_lev + 1;
    
    //error value
    std::cout << "\t\t";
    std::cout.precision(14);

    std::cout << norm[solution][mesh_lev] << "\t";

    std::cout << std::endl;
  
    //convergence rate value
    std::cout << "\t\t" /*<<  std::fixed*/;
    std::cout.precision(3);

    std::cout << log( norm[solution][mesh_lev] / norm[solution][mesh_lev + 1] ) / log(2.);

    std::cout << std::endl;
      
    }
    
  
  
}


template < class real_num>
/*static */ void FE_convergence< real_num >::output_convergence_order_all(const std::vector< Unknown > &  unknowns,
                                        const std::vector < std::vector < real_num > >  &  norms, 
                                        const unsigned sobolev_norms,
                                        const unsigned volume_or_boundary,
                                        const unsigned convergence_computation_method) {
    
    
    assert( unknowns.size() == norms.size() );
    
    std::cout << std::endl;
      
        
     std::cout << "==== Convergence computation method: " << convergence_computation_names[convergence_computation_method] << std::endl;
     
     std::cout << "==== Volume = 0, boundary = 1: here we have " << volume_or_boundary_names[volume_or_boundary] << std::endl;

     std::cout << "==== ERROR and ORDER OF CONVERGENCE for the norm " << norm_names[sobolev_norms] << std::endl;

    
     for (unsigned int u = 0; u < unknowns.size(); u++) {
         
            std::cout << unknowns[u]._name          << " , " << 
            " FE Family " << unknowns[u]._fe_family << " , " <<
            " FE Order "  << unknowns[u]._fe_order  << std::endl;
            
            for (int lev = 0; lev < norms[u].size(); lev++) {
                output_convergence_order(norms, u, lev);
            }
            std::cout << std::endl;
            
       }

      
 }
     


template < class real_num>
/*static */ void FE_convergence< real_num >::output_convergence_rate(double norm_i,
                                                                     double norm_ip1,
                                                                     std::string norm_name,
                                                                     unsigned maxNumberOfMeshes,
                                                                     int loop_i) {

    std::cout << loop_i + 1 << "\t\t" <<  std::setw(11) << std::setprecision(10) << norm_i << "\t\t\t\t" ;
  
    if (loop_i < maxNumberOfMeshes/*norm.size()*/ - 2) {
      std::cout << std::setprecision(3) << log( norm_i/ norm_ip1 ) / log(2.) << std::endl;
    }
  
} 

 
 

template < class real_num>
/*static*/  std::vector< real_num > FE_convergence< real_num >::compute_error_norms_volume_or_boundary_with_analytical_sol_or_not(
                                                                            const std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, double> *  > > & elem_all,
                                                                            const std::vector<Gauss> & quad_rules,
                                                                            const MultiLevelSolution* ml_sol_single_level,
                                                                            const MultiLevelSolution* ml_sol_all_levels_needed_for_incremental,
                                                                            const std::string & unknown,
                                                                            const Math::Function< real_num > * ex_sol_in,
                                                                            const unsigned current_level,
                                                                            const unsigned sobolev_norms,
                                                                            const unsigned convergence_rate_computation_method,
                                                                            const unsigned volume_or_boundary
                                             ) {
     

  
// ||u_h - u_(h/2)||/||u_(h/2)-u_(h/4)|| = 2^alpha, alpha is order of conv 
//i.e. ||prol_(u_(i-1)) - u_(i)|| = err(i) => err(i-1)/err(i) = 2^alpha ,implemented as log(err(i)/err(i+1))/log2

   if (convergence_rate_computation_method == 1 && ex_sol_in == NULL) { std::cout << "Please provide analytical solution" << std::endl; abort(); }
   
    
  const unsigned num_norms = norm_names.size();

  //norms that we are computing here //first L2, then H1 - BEGIN ============
  std::vector< real_num > norms_exact_function_at_qp(num_norms);                  std::fill(norms_exact_function_at_qp.begin(), norms_exact_function_at_qp.end(), 0.);   
  std::vector< real_num > norms_exact_function_at_dofs(num_norms);       std::fill(norms_exact_function_at_dofs.begin(), norms_exact_function_at_dofs.end(), 0.);
  std::vector< real_num > norms_inexact_dofs(num_norms);     std::fill(norms_inexact_dofs.begin(), norms_inexact_dofs.end(), 0.);
  //norms that we are computing here //first L2, then H1 - END ============
  
  
  unsigned level = ml_sol_single_level->GetMLMesh()->GetNumberOfLevels() - 1u; //this is supposed to be zero because a single level is used
  
  //  extract pointers to the several objects that we are going to use
  const Mesh*     msh = ml_sol_single_level->GetMLMesh()->GetLevel(level);
  const Solution* sol = ml_sol_single_level->GetSolutionLevel(level);

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)

 unsigned iproc = msh->processor_id();

  
 //***************************************************  
  constexpr unsigned int space_dim = 3;
  const unsigned int dim_offset_grad = 3;

  // ======================================
  // reserve memory for the local standar vectors
  const unsigned max_size = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  
   
  //solution variable
  unsigned sol_uIndex = ml_sol_single_level->GetIndex(unknown.c_str());
  unsigned sol_uType  = ml_sol_single_level->GetSolutionType(sol_uIndex);    // get the finite element type for "u"

  

  

  CurrentElem < double > geom_element_iel(dim, msh);

  
//-----------------  VOLUME, common to both - BEGIN
  std::vector < real_num >  sol_u;                               sol_u.reserve(max_size);
  std::vector < real_num >  sol_u_exact_at_dofs;   sol_u_exact_at_dofs.reserve(max_size);
  std::vector < real_num >  sol_u_inexact_prolongated_from_coarser_level;     sol_u_inexact_prolongated_from_coarser_level.reserve(max_size);

  
  
  unsigned solType_coords = CONTINUOUS_BIQUADRATIC;
  std::vector < std::vector < real_num > > x(dim_offset_grad);    // local coordinates
  for (unsigned i = 0; i < x.size(); i++)   x[i].reserve(max_size);

//-----------------  VOLUME, common to both - END

  
//-----------------  VOLUME initialization - BEGIN
 
//***************************************************  
   std::vector < std::vector < double > >  JacI_qp(space_dim);
   std::vector < std::vector < double > >  Jac_qp(dim);
    for (unsigned d = 0; d < dim; d++) { Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < space_dim; d++) { JacI_qp[d].resize(dim); }
    
  real_num detJac_qp = 0.;
  real_num weight = 0.; // gauss point weight
 //***************************************************  


  std::vector < real_num > phi;
  std::vector < real_num > phi_x;
  
  phi.reserve(max_size);
  phi_x.reserve(max_size * dim_offset_grad);
  
  
  
  std::vector < real_num > phi_coords;
  std::vector < real_num > phi_coords_x;

  phi_coords.reserve(max_size);
  phi_coords_x.reserve(max_size * dim_offset_grad);


 

//-----------------  VOLUME initialization - END

  
  
  
//-----------------  BOUNDARY initialization - BEGIN
    
 //***************************************************  
     std::vector < std::vector < double > >  JacI_qp_bdry(space_dim);
     std::vector < std::vector < double > >  Jac_qp_bdry(dim-1);
    for (unsigned d = 0; d < Jac_qp_bdry.size(); d++) {   Jac_qp_bdry[d].resize(space_dim); }
    for (unsigned d = 0; d < JacI_qp_bdry.size(); d++) { JacI_qp_bdry[d].resize(dim-1); }
    
    real_num detJac_iqp_bdry = 0.;
  real_num weight_iqp_bdry = 0.;
 //***************************************************  
    
  std::vector <double> phi_u_bdry;  
  std::vector <double> phi_u_x_bdry; 

  phi_u_bdry.reserve(max_size);
  phi_u_x_bdry.reserve(max_size * dim_offset_grad);
  
  
  std::vector < real_num > phi_coords_bdry;
  std::vector < real_num > phi_coords_x_bdry;

  phi_coords_bdry.reserve(max_size);
  phi_coords_x_bdry.reserve(max_size * dim_offset_grad);

//-----------------  BOUNDARY initialization - END

  
  
  
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

//---------------- Geometry - BEGIN
      geom_element_iel.set_coords_at_dofs_and_geom_type(iel, solType_coords);
        
    unsigned nDofx = msh->GetElementDofNumber(iel, solType_coords);
    
    for (int i = 0; i < x.size(); i++) {
      x[i].resize(nDofx);
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, solType_coords);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim_offset_grad/*dim*/; jdim++) {
        x[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);  // global extraction and local storage for the element coordinates
      }
    }
    
//---------------- Geometry - END
    
    
//---------------- Solution - BEGIN
    unsigned nDofu  = msh->GetElementDofNumber(iel, sol_uType);

    // resize local arrays
    sol_u.resize(nDofu);
    sol_u_exact_at_dofs.resize(nDofu);
    sol_u_inexact_prolongated_from_coarser_level.resize(nDofu);

    
    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
        
        std::vector< real_num > x_at_dof(dim_offset_grad, 0.);
        for (unsigned jdim = 0; jdim < x_at_dof.size(); jdim++) x_at_dof[jdim] = x[jdim][i];
        
      unsigned solDof = msh->GetSolutionDof(i, iel, sol_uType);
                   sol_u[i]  =                                        (*sol->_Sol[sol_uIndex])(solDof);
      sol_u_inexact_prolongated_from_coarser_level[i]  = (*ml_sol_all_levels_needed_for_incremental->GetSolutionLevel(current_level)->_Sol[sol_uIndex])(solDof);
      if (ex_sol_in != NULL) sol_u_exact_at_dofs[i] = ex_sol_in->value(x_at_dof);
    }
//---------------- Solution - END


// *** VOLUME - BEGIN ***
 if (volume_or_boundary == 0) {  
      
      
    const short unsigned ielGeom = geom_element_iel.geom_type();
    
   
    for (unsigned ig = 0; ig < quad_rules[ielGeom].GetGaussPointsNumber(); ig++) {

      // *** get gauss point weight, test function and test function partial derivatives ***
	elem_all[ielGeom][solType_coords]->JacJacInv(x, ig, Jac_qp, JacI_qp, detJac_qp, space_dim);
    elem_all[ielGeom][sol_uType]->shape_funcs_current_elem(ig, JacI_qp, phi, phi_x, boost::none, space_dim);
    
    elem_all[ielGeom][solType_coords]->shape_funcs_current_elem(ig, JacI_qp, phi_coords, phi_coords_x,  boost::none, space_dim);
    
    weight = detJac_qp * quad_rules[ielGeom].GetGaussWeightsPointer()[ig];

    
		 //========== QP eval - BEGIN ===============================================
      real_num sol_u_gss = 0.;
      real_num exactSol_from_dofs_gss = 0.;
      real_num sol_u_inexact_prolongated_from_coarser_level_gss = 0.;
      
      std::vector < real_num > gradSolu_gss(dim_offset_grad, 0.);
      std::vector < real_num > gradSolu_exact_at_dofs_gss(dim_offset_grad, 0.);
      std::vector < real_num > gradSolu_inexact_prolongated_from_coarser_level_gss(dim_offset_grad, 0.);


      for (unsigned i = 0; i < phi.size(); i++) {
        sol_u_gss                += phi[i] * sol_u[i];
        exactSol_from_dofs_gss  += phi[i] * sol_u_exact_at_dofs[i];
        sol_u_inexact_prolongated_from_coarser_level_gss   += phi[i] * sol_u_inexact_prolongated_from_coarser_level[i];

        for (unsigned jdim = 0; jdim < dim_offset_grad; jdim++) {
          gradSolu_gss[jdim]                += phi_x[i * dim_offset_grad + jdim] * sol_u[i];
          gradSolu_exact_at_dofs_gss[jdim]  += phi_x[i * dim_offset_grad + jdim] * sol_u_exact_at_dofs[i];
          gradSolu_inexact_prolongated_from_coarser_level_gss[jdim]   += phi_x[i * dim_offset_grad + jdim] * sol_u_inexact_prolongated_from_coarser_level[i];
        }
        
 
      }
      
      
      std::vector < real_num > x_gss(dim_offset_grad, 0.);
      
       for (unsigned i = 0; i < phi_coords.size(); i++) {
          for (unsigned jdim = 0; jdim < dim_offset_grad; jdim++) {
           x_gss[jdim] += x[jdim][i] * phi_coords[i];
          }
       }
		 //========== QP eval - END ===============================================
      
      
      
      

// H^0 ==============      
    if (sobolev_norms == 0) {
      real_num exactSol = 0.; if (ex_sol_in != NULL) exactSol = ex_sol_in->value(x_gss);

      norms_exact_function_at_qp[0] += (sol_u_gss - exactSol)                * (sol_u_gss - exactSol)       * weight;
      norms_exact_function_at_dofs[0]     += (sol_u_gss - exactSol_from_dofs_gss)  * (sol_u_gss - exactSol_from_dofs_gss) * weight;
      norms_inexact_dofs[0]   += (sol_u_gss - sol_u_inexact_prolongated_from_coarser_level_gss)   * (sol_u_gss - sol_u_inexact_prolongated_from_coarser_level_gss)  * weight;
    }
    
// H^1 ==============
    else if (sobolev_norms == 1) {
      std::vector < real_num > exactGradSol(dim_offset_grad,0.);    if (ex_sol_in != NULL) exactGradSol = ex_sol_in->gradient(x_gss);

      for (unsigned d = 0; d < dim_offset_grad ; d++) {
        norms_exact_function_at_qp[1] += ((gradSolu_gss[d] - exactGradSol[d])               * (gradSolu_gss[d]  - exactGradSol[d])) * weight;
        norms_exact_function_at_dofs[1]     += ((gradSolu_gss[d] - gradSolu_exact_at_dofs_gss[d]) * (gradSolu_gss[d] - gradSolu_exact_at_dofs_gss[d])) * weight;
        norms_inexact_dofs[1]   += ((gradSolu_gss[d] - gradSolu_inexact_prolongated_from_coarser_level_gss[d])  * (gradSolu_gss[d] - gradSolu_inexact_prolongated_from_coarser_level_gss[d]))  * weight;
      }
   }
      
   } // end gauss point loop
  }
// *** VOLUME - END ***
  

  
// *** BOUNDARY - BEGIN ***
else if (volume_or_boundary == 1 )	{
  
	  std::vector<double> normal_at_qp_bdry(space_dim, 0.);
	       
	  for(unsigned iface = 0; iface < msh->GetElementFaceNumber(iel); iface++) {
          
       const unsigned ielGeom_bdry = msh->GetElementFaceType(iel, iface);    
       
// ----------
       geom_element_iel.set_coords_at_dofs_bdry_3d(iel, iface, solType_coords);
 
// ----------

		            const int bdry_index_j = msh->GetMeshElements()->GetFaceElementIndex(iel, iface);

	    if(  bdry_index_j < 0 ) {

	
		
		for(unsigned ig_bdry = 0; ig_bdry < quad_rules[ielGeom_bdry].GetGaussPointsNumber(); ig_bdry++) {
		  
    elem_all[ielGeom_bdry][solType_coords]->JacJacInv(geom_element_iel.get_coords_at_dofs_bdry_3d(), ig_bdry, Jac_qp_bdry, JacI_qp_bdry, detJac_iqp_bdry, space_dim);
    weight_iqp_bdry = detJac_iqp_bdry * quad_rules[ielGeom_bdry].GetGaussWeightsPointer()[ig_bdry];
    
    elem_all[ielGeom_bdry][sol_uType]     ->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_u_bdry, phi_u_x_bdry, boost::none, space_dim);
    elem_all[ielGeom_bdry][solType_coords]->shape_funcs_current_elem(ig_bdry, JacI_qp_bdry, phi_coords_bdry, phi_coords_x_bdry,  boost::none, space_dim);
	elem_all[ielGeom_bdry][solType_coords]->compute_normal(Jac_qp_bdry, normal_at_qp_bdry);

		  
		 //========== QP eval - BEGIN ===============================================
      real_num		  sol_u_bdry_gss = 0.;
      real_num exactSol_from_dofs_bdry_gss = 0.;
      real_num sol_u_inexact_prolongated_from_coarser_level_gss = 0.;
          
      std::vector < real_num >  sol_u_x_bdry_gss(space_dim, 0.);
      std::vector < real_num >  gradSolu_exact_at_dofs_bdry_gss(dim_offset_grad, 0.);
      std::vector < real_num >  gradSolu_inexact_prolongated_from_coarser_level_bdry_gss(dim_offset_grad, 0.);

                  
		      for (int i_bdry = 0; i_bdry < phi_u_bdry.size(); i_bdry++)  {
		    unsigned int i_vol = msh->GetLocalFaceVertexIndex(iel, iface, i_bdry);
			
			sol_u_bdry_gss                                    +=  sol_u[i_vol]              * phi_u_bdry[i_bdry];
            exactSol_from_dofs_bdry_gss                       += sol_u_exact_at_dofs[i_vol] * phi_u_bdry[i_bdry];
            
            sol_u_inexact_prolongated_from_coarser_level_gss  += sol_u_inexact_prolongated_from_coarser_level[i_vol]              * phi_u_bdry[i_bdry];
        
            for (int d = 0; d < space_dim; d++) {
			      sol_u_x_bdry_gss[d] += sol_u[i_vol] * phi_u_x_bdry[i_bdry * space_dim + d];
                  gradSolu_exact_at_dofs_bdry_gss[d] += sol_u_exact_at_dofs[i_vol] * phi_u_x_bdry[i_bdry * space_dim + d];
			      gradSolu_inexact_prolongated_from_coarser_level_bdry_gss[d] += sol_u_inexact_prolongated_from_coarser_level[i_vol] * phi_u_x_bdry[i_bdry * space_dim + d];
			    }
		    }
		      
		      double laplace_u_bdry = 0.;  for (int d = 0; d < space_dim; d++) { laplace_u_bdry += sol_u_x_bdry_gss[d] * sol_u_x_bdry_gss[d]; }


       std::vector < real_num > x_gss_bdry(dim_offset_grad, 0.);
       
    for (unsigned i = 0; i < phi_coords_bdry.size(); i++) {
          for (unsigned jdim = 0; jdim < dim_offset_grad; jdim++) {
           x_gss_bdry[jdim] += geom_element_iel.get_coords_at_dofs_bdry_3d()[jdim][i] * phi_coords_bdry[i];
          }
       }
		 //========== QP eval - END ===============================================

// H^0 - BEGIN ============== 
 if (sobolev_norms == 0) {
      real_num exactSol_bdry = 0.; if (ex_sol_in != NULL) exactSol_bdry = ex_sol_in->value(x_gss_bdry);
      
      norms_exact_function_at_qp[0] += (sol_u_bdry_gss - exactSol_bdry) * (sol_u_bdry_gss - exactSol_bdry) * weight_iqp_bdry;
      
      norms_exact_function_at_dofs[0]  += (sol_u_bdry_gss - exactSol_from_dofs_bdry_gss)  * 
                                          (sol_u_bdry_gss - exactSol_from_dofs_bdry_gss) * weight;

      norms_inexact_dofs[0]   += (sol_u_bdry_gss - sol_u_inexact_prolongated_from_coarser_level_gss)   * (sol_u_bdry_gss - sol_u_inexact_prolongated_from_coarser_level_gss)  * weight_iqp_bdry;
  }
// H^0 - END ==============      
                  
                  
                  
// H^1 - BEGIN ==============
    else if (sobolev_norms == 1) {
      
      std::cout << "Here you have to ROTATE the boundary element appropriately, I think" << std::endl;
      
      std::vector < real_num > exactGradSol_bdry(dim_offset_grad, 0.);    if (ex_sol_in != NULL) exactGradSol_bdry = ex_sol_in->gradient(x_gss_bdry);  
      ///@todo THIS IS WRONG! It is NOT the SURFACE GRADIENTTTT, so we have to be careful when we have more than one boundary face!!!
      
      std::vector < real_num > exactGradSol_bdry_dot_tangent =  Math::tangent_vector_from_normal( exactGradSol_bdry, normal_at_qp_bdry, exactGradSol_bdry.size() );
      
      std::vector < real_num > sol_u_x_bdry_gss_dot_tangent =  Math::tangent_vector_from_normal( sol_u_x_bdry_gss, normal_at_qp_bdry, sol_u_x_bdry_gss.size() );

      std::vector < real_num > gradSolu_exact_at_dofs_bdry_gss_dot_tangent =  Math::tangent_vector_from_normal( gradSolu_exact_at_dofs_bdry_gss, normal_at_qp_bdry, gradSolu_exact_at_dofs_bdry_gss.size() );
      
      std::vector < real_num > gradSolu_inexact_prolongated_from_coarser_level_bdry_gss_dot_tangent =  Math::tangent_vector_from_normal( gradSolu_inexact_prolongated_from_coarser_level_bdry_gss, normal_at_qp_bdry, gradSolu_inexact_prolongated_from_coarser_level_bdry_gss.size() );

      
      for (unsigned d = 0; d < dim_offset_grad ; d++) {
        norms_exact_function_at_qp[1] += ( (sol_u_x_bdry_gss_dot_tangent[d] - exactGradSol_bdry_dot_tangent[d])  *  
                                     (sol_u_x_bdry_gss_dot_tangent[d] - exactGradSol_bdry_dot_tangent[d]) ) * weight_iqp_bdry;
        
        norms_exact_function_at_dofs[1]   += (( sol_u_x_bdry_gss_dot_tangent[d] - gradSolu_exact_at_dofs_bdry_gss_dot_tangent[d]) * 
                                              ( sol_u_x_bdry_gss_dot_tangent[d] - gradSolu_exact_at_dofs_bdry_gss_dot_tangent[d])) * weight_iqp_bdry;
                                     
        norms_inexact_dofs[1]   += ( (sol_u_x_bdry_gss_dot_tangent[d] - gradSolu_inexact_prolongated_from_coarser_level_bdry_gss_dot_tangent[d])   *
                                     (sol_u_x_bdry_gss_dot_tangent[d] - gradSolu_inexact_prolongated_from_coarser_level_bdry_gss_dot_tangent[d]) ) * weight_iqp_bdry;
      }
   }
// H^1 - END ==============
                  
        }
            
        
	      } //end face
	      
	  }  // loop over element faces   
	  

//=====================================================================================================================  
//=====================================================================================================================  
//=====================================================================================================================  
     
    
}   
// *** BOUNDARY - END ***

  
  
  } //end element loop for each process

  
// *** Add norms to all processes - BEGIN ***

  NumericVector* norm_vec_exact_using_qp;
                 norm_vec_exact_using_qp = NumericVector::build().release();
                 norm_vec_exact_using_qp->init(msh->n_processors(), 1 , false, AUTOMATIC);

         if (sobolev_norms == 0) { norm_vec_exact_using_qp->set(iproc, norms_exact_function_at_qp[0]);  norm_vec_exact_using_qp->close();  norms_exact_function_at_qp[0] = norm_vec_exact_using_qp->l1_norm(); }
    else if (sobolev_norms == 1) { norm_vec_exact_using_qp->set(iproc, norms_exact_function_at_qp[1]);  norm_vec_exact_using_qp->close();  norms_exact_function_at_qp[1] = norm_vec_exact_using_qp->l1_norm(); }

          delete norm_vec_exact_using_qp;

   // add the norms of all processes
  NumericVector* norm_vec_exact_using_dofs;
                 norm_vec_exact_using_dofs = NumericVector::build().release();
                 norm_vec_exact_using_dofs->init(msh->n_processors(), 1 , false, AUTOMATIC);

         if (sobolev_norms == 0) { norm_vec_exact_using_dofs->set(iproc, norms_exact_function_at_dofs[0]);  norm_vec_exact_using_dofs->close();  norms_exact_function_at_dofs[0] = norm_vec_exact_using_dofs->l1_norm(); }
    else if (sobolev_norms == 1) { norm_vec_exact_using_dofs->set(iproc, norms_exact_function_at_dofs[1]);  norm_vec_exact_using_dofs->close();  norms_exact_function_at_dofs[1] = norm_vec_exact_using_dofs->l1_norm(); }

          delete norm_vec_exact_using_dofs;

  // add the norms of all processes
  NumericVector* norm_vec_inexact;
                 norm_vec_inexact = NumericVector::build().release();
                 norm_vec_inexact->init(msh->n_processors(), 1 , false, AUTOMATIC);

         if (sobolev_norms == 0) { norm_vec_inexact->set(iproc, norms_inexact_dofs[0]);  norm_vec_inexact->close();  norms_inexact_dofs[0] = norm_vec_inexact->l1_norm(); }
    else if (sobolev_norms == 1) { norm_vec_inexact->set(iproc, norms_inexact_dofs[1]);  norm_vec_inexact->close();  norms_inexact_dofs[1] = norm_vec_inexact->l1_norm(); }

          delete norm_vec_inexact;

          
    for (int n = 0; n < norms_exact_function_at_qp.size(); n++) { 
     norms_exact_function_at_qp[n] = sqrt(norms_exact_function_at_qp[n]);                                            
   norms_exact_function_at_dofs[n] = sqrt(norms_exact_function_at_dofs[n]);                                            
             norms_inexact_dofs[n] = sqrt(norms_inexact_dofs[n]);
    }

// *** Add norms to all processes - END ***
    
    
if (convergence_rate_computation_method == 0)  return norms_inexact_dofs;
else if (convergence_rate_computation_method == 1)  return norms_exact_function_at_qp;
// if (convergence_rate_computation_method == 1)  return   norms_exact_function_at_dofs; //@todo this one does not work in the ABSOLUTE case but only in the INCREMENTAL case
else  abort();
 
} 
 
 
 

// 
// If the method uses the exact solution, in the vector of norms I will put
// || u_h - u_{true} ||       coarsest level
// || u_{h/2} - u_{true} ||
// || u_{h/4} - u_{true} ||
// 
//  Every u_h is the C^0 LAGRANGE INTERPOLANT at that level of the original function
//  If I evaluate the analytical true solution at the real Quadrature Point, it is like using a BETTER INTERPOLANT,
// that is the interpolant with ORTHOGONAL POLYNOMIALS that is under the construction of Gauss quadrature rules 
// So actually every level is of the type 
// || u_h - u_{true} || = || I_{h, Lagrange} u  -  I_{h, Orthogonal Polynomials} u  ||
// 
// These observations tell us that if we compared the ALGEBRAIC VECTORS (which only contain DOFS of LAGRANGE INTERPOLANT)
// we would probably get ALL ZEROs, so in that case we could only do it in INCREMENTAL way
///  @todo do it  || x_h - x_{h/2} || where the vectors are algebraic

// If the method is incremental, in the vector of norms I will put
// || u_h     - u_{h/2} ||       coarsest level and 1st finer
// || u_{h/2} - u_{h/4} ||
// || u_{h/4} - u_{h/8} ||
// So, it is clear that at each stage I involve 2 different levels instead of one


template < class real_num>
/*static*/  void FE_convergence< real_num >::compute_error_norms_per_unknown_per_level(
                                     const std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, double> *  > > & elem_all,
                                                                                   const std::vector<Gauss> & quad_rules,
                                                                                   const MultiLevelSolution* ml_sol_single_level,
                                                                                   MultiLevelSolution* ml_sol_all_levels_needed_for_incremental,
                                                                                   const std::vector< Unknown > &  unknowns,
                                                                                   const  std::vector< Math::Function< real_num > * > &  ex_sol_in,
                                                                                   std::vector < std::vector < real_num > >  &  norms,
                                                                                   const unsigned lev,
                                                                                   const unsigned sobolev_norms,
                                                                                   const unsigned volume_or_boundary,
                                                                                   const unsigned convergence_rate_computation_method
                                         ) {
     
    
        // ======= Control - BEGIN
            if (volume_or_boundary == 0) {
            }
            else if (volume_or_boundary == 1) {
            }
            else { std::cout << "Boundary of boundary not implemented in function " << __func__; abort(); }
        // ======= Control - END
            
            
   
        if ( lev > 0 ) {

            // ======= prolongate to the current level (lev) from the coarser level (lev-1) (so that you can compare the two) ========================
            ml_sol_all_levels_needed_for_incremental->RefineSolution(lev);
            
            // =======  compute the error norm at the current level (lev) ========================
      
            for (unsigned int u = 0; u < unknowns.size(); u++) {  //this loop could be inside the below function
                
            std::vector< real_num > norm_out;
            
            norm_out = FE_convergence::compute_error_norms_volume_or_boundary_with_analytical_sol_or_not (elem_all,
                                                                                                   quad_rules, 
                                                                                                   ml_sol_single_level,
                                                                                                   ml_sol_all_levels_needed_for_incremental, 
                                                                                                   unknowns[u]._name,
                                                                                                   ex_sol_in[u], 
                                                                                                   lev,
                                                                                                   sobolev_norms, 
                                                                                                   convergence_rate_computation_method, 
                                                                                                   volume_or_boundary);
            
            
            norms[u][lev-1] = norm_out[sobolev_norms];
                                       
                   }
                   

        }         
                 
              // ======= store the last computed solution to prepare the next iteration (the current level i is now overwritten) ========================
              const unsigned level_to_pick_from = ml_sol_single_level->GetMLMesh()->GetNumberOfLevels() - 1;
            ml_sol_all_levels_needed_for_incremental->fill_at_level_from_level(lev, level_to_pick_from, *ml_sol_single_level);
        
                 
                 
                 
 } 
 



 
template < class real_num >
void compute_L2_norm_of_errors_of_unknowns_with_analytical_sol(MultiLevelProblem& ml_prob,
                                                              const unsigned int level_in,
                                                              const unsigned int system_size_in ) {


//***************** Check Number of Unknowns - BEGIN **********************************  
  const unsigned int n_vars = system_size_in;
  if (n_vars != 1)  { std::cout << "Function written for only 1 scalar unknown";  abort();  }
//***************** Check Number of Unknowns - END **********************************  

//***************** Check that true solution is provided - BEGIN **********************************  
  if (ml_prob.get_app_specs_pointer()->_true_solution_function == NULL) {  std::cout << "No true solution provided";  abort();  }
//***************** Check that true solution is provided - END **********************************  


//***************** Level - BEGIN **********************************  
  const unsigned level = level_in;
//***************** Level - END **********************************  
  

  Mesh*                    msh = ml_prob._ml_msh->GetLevel(level);
  elem*                     el = msh->GetMeshElements();

  MultiLevelSolution*    ml_sol = ml_prob._ml_sol;
  Solution*                sol = ml_prob._ml_sol->GetSolutionLevel(level);

  const unsigned  dim = msh->GetDimension();
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));

  unsigned    iproc = msh->processor_id();



  //=============== Geometry ========================================
  unsigned xType = CONTINUOUS_BIQUADRATIC; // the FE for the domain need not be biquadratic
  
  CurrentElem < double > geom_element(dim, msh);            // must be adept if the domain is moving, otherwise double
    
  constexpr unsigned int space_dim = 3;
//***************************************************  


 //******************** quadrature *******************************  
  double jacXweight_qp; 

 //********************* quadrature, geometry *********************** 
 //***************************************************  
  std::vector <double> phi_coords;
  std::vector <double> phi_coords_x; 
  std::vector <double> phi_coords_xx;

  phi_coords.reserve(maxSize);
  phi_coords_x.reserve(maxSize * space_dim);
  phi_coords_xx.reserve(maxSize * dim2);
  
  
 //********************* quadrature, unknowns *********************** 
 //***************************************************  
  std::vector <double> phi_u;
  std::vector <double> phi_u_x; 
  std::vector <double> phi_u_xx;

  phi_u.reserve(maxSize);
  phi_u_x.reserve(maxSize * space_dim);
  phi_u_xx.reserve(maxSize * dim2);
  
  const std::string solname_u = ml_sol->GetSolName_string_vec()[0];
  unsigned solIndex_u = ml_sol->GetIndex(solname_u.c_str()); 
  unsigned solFEType_u = ml_sol->GetSolutionType(solIndex_u); 


  std::vector < double >  sol_u;     sol_u.reserve(maxSize);
  std::vector < double >  sol_u_true;     sol_u_true.reserve(maxSize);
 //***************************************************  
 //***************************************************  


  
 //***************************************************  
     std::vector < std::vector < double /*real_num_mov*/ > >  JacI_qp(space_dim);
     std::vector < std::vector < double /*real_num_mov*/ > >  Jac_qp(dim);
    for (unsigned d = 0; d < dim; d++) { Jac_qp[d].resize(space_dim); }
    for (unsigned d = 0; d < space_dim; d++) { JacI_qp[d].resize(dim); }
    
    double /*real_num_mov*/ detJac_qp;
    

  //prepare Abstract quantities for all fe fams for all geom elems: all quadrature evaluations are performed beforehand in the main function
  std::vector < std::vector < /*const*/ elem_type_templ_base<real_num, double/*real_num_mov*/ > *  > > elem_all;
  ml_prob.get_all_abstract_fe(elem_all);
 //***************************************************  
  
  double  norm_dof_based = 0.;
  double  norm_qp_based = 0.;
  

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {
      

    geom_element.set_coords_at_dofs_and_geom_type(iel, xType);
        
    const short unsigned ielGeom = geom_element.geom_type();

 //**************** state **************************** 
    unsigned nDof_u     = msh->GetElementDofNumber(iel, solFEType_u);
    sol_u    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
     unsigned solDof_u = msh->GetSolutionDof(i, iel, solFEType_u);
      sol_u[i] = (*sol->_Sol[solIndex_u])(solDof_u);
    }
    //***************************************************  
    
 //**************** true sol **************************** 
    sol_u_true    .resize(nDof_u);
   // local storage of global mapping and solution
    for (unsigned i = 0; i < sol_u.size(); i++) {
        std::vector< double > xyz_i(dim);
             	for (unsigned d = 0; d < xyz_i.size(); d++) {
                   xyz_i[d] = geom_element.get_coords_at_dofs_3d()[d][i];
                }
               sol_u_true[i]  =   ml_prob.get_app_specs_pointer()->_true_solution_function->value( xyz_i);  
    }

    //***************************************************  
 
    
 //*************************************************** 
    

 
 //========= VOLUME ==================   
   
 //========= gauss value quantities ==================   
	std::vector<double> sol_u_x_gss(space_dim);     std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
 //===================================================   
    
    
    
      // *** Quadrature point loop ***
      for (unsigned i_qp = 0; i_qp < ml_prob.GetQuadratureRule(ielGeom).GetGaussPointsNumber(); i_qp++) {
          
        // *** get gauss point weight, test function and test function partial derivatives ***
// 	msh->_finiteElement[ielGeom][solFEType_u]->Jacobian(geom_element.get_coords_at_dofs_3d(),    i_qp, weight,    phi_u,    phi_u_x,    boost::none /*phi_u_xx*/);
          
	elem_all[ielGeom][xType]->JacJacInv(geom_element.get_coords_at_dofs_3d(), i_qp, Jac_qp, JacI_qp, detJac_qp, space_dim);
    jacXweight_qp = detJac_qp * ml_prob.GetQuadratureRule(ielGeom).GetGaussWeightsPointer()[i_qp];
    elem_all[ielGeom][solFEType_u]->shape_funcs_current_elem(i_qp, JacI_qp, phi_u, phi_u_x, boost::none /*phi_u_xx*/, space_dim);
    elem_all[ielGeom][xType]->shape_funcs_current_elem(i_qp, JacI_qp, phi_coords, phi_coords_x, boost::none /*phi_coords_xx*/, space_dim);


//--------------    
	std::fill(sol_u_x_gss.begin(), sol_u_x_gss.end(), 0.);
	
	for (unsigned i = 0; i < sol_u.size(); i++) {
// 	                                                sol_u_gss      += sol_u[i] * phi_u[i];
                   for (unsigned d = 0; d < sol_u_x_gss.size(); d++)   sol_u_x_gss[d] += sol_u[i] * phi_u_x[i * space_dim + d];
          }
//--------------    

//--------------    
 /// @assignment You need to evaluate your manufactured right hand side at the quadrature point qp.
 /// Hence, you need to compute the coordinates of the quadrature point. Let us call them x_qp.
 /// These are obtained just like every quantity at a quadrature point, i.e., by interpolating the values of the quantity at the element nodes.
 /// The interpolation is performed by using the shape functions.
 /// In other words, 
 ///         (x_qp) = summation of (x_nodes) * (shape function of that node, evaluated at qp)
 /// 
 ///   (x_nodes) are obtained from   geom_element.get_coords_at_dofs_3d()  (this is a  std::vector < vector >,  where the outer index is the dimension and the inner index ranges over the nodes) 
 ///   (shape function of that node, evaluated at qp)  is obtained from phi_u  (this is a vector, whose index ranges over the nodes)
 
 std::vector<double> x_qp(dim, 0.);
          
        for (unsigned i = 0; i < phi_coords.size(); i++) {
          	for (unsigned d = 0; d < dim; d++) {
	                                                x_qp[d]    += geom_element.get_coords_at_dofs_3d()[d][i] * phi_coords[i]; // fetch of coordinate points
            }
        }
 
//--------------    

          
	// *** phi_i loop ***
	double u_qp = 0.;        for (unsigned i = 0; i < phi_u.size(); i++) {            u_qp +=  sol_u[i] * phi_u[i];        }

 	double u_true_qp = 0.;   for (unsigned i = 0; i < phi_u.size(); i++) {            u_true_qp +=  sol_u_true[i] * phi_u[i];        }
       
  
	      
       norm_dof_based +=  jacXweight_qp * (  u_qp - u_true_qp) * (  u_qp - u_true_qp) ;
       norm_qp_based  +=  jacXweight_qp * (  u_qp -  ml_prob.get_app_specs_pointer()->_true_solution_function->value( x_qp )) * (  u_qp - ml_prob.get_app_specs_pointer()->_true_solution_function->value( x_qp )) ;
//======================Residuals=======================
	      
// // //           if (assembleMatrix) {
// // // 	    
// // //             // *** phi_j loop ***
// // //             for (unsigned j = 0; j < nDof_max; j++) {
// // // 
// // // 
// // // //--------------    
// // //               double laplace_mat_du_u_i_j = 0.;
// // // 
// // //               if ( i < nDof_u && j < nDof_u ) {
// // //                   for (unsigned kdim = 0; kdim < space_dim; kdim++) {
// // //                          laplace_mat_du_u_i_j           += phi_u_x   [i * space_dim + kdim] *
// // //                                                            phi_u_x   [j * space_dim + kdim];
// // // 	              }
// // //              }
// // // //--------------    
// // // 
// // // 
// // // 
// // //               //============ delta_state row ============================
// // //               //DIAG BLOCK delta_state - state
// // // 		  if ( i < nDof_u && j < nDof_u )       Jac[ (0 + i) * nDof_AllVars   + 	(0 + j) ]  += jacXweight_qp * laplace_mat_du_u_i_j;
// // // // 		  if ( i < nDof_u && j < nDof_u )       Jac[ (0 + i) * nDof_AllVars   + 	(0 + j) ]  += jacXweight_qp * laplace_beltrami_mat_du_u_i_j; ///@todo On a flat domain, this must coincide with the standard Laplacian, so we can do a double check with this
// // //             } // end phi_j loop
// // //           } // endif assemble_matrix

        
      } // end gauss point loop


   
  } //end element loop for each process

   double norm_dof_based_parallel = 0.; MPI_Allreduce( &norm_dof_based, &norm_dof_based_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   double norm_qp_based_parallel = 0.; MPI_Allreduce( &norm_qp_based, &norm_qp_based_parallel, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
 
  std::cout << std::scientific << std::setw(20) << std::setprecision(15) << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& L2-Norm of the error, dof-based: " << norm_dof_based_parallel << std::endl;
  
  std::cout << std::scientific << std::setw(20) << std::setprecision(15) << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& L2-Norm of the error, qp-based " << norm_qp_based_parallel << std::endl;



  return;
}




//This function does not require an equation
std::pair < double, double > GetErrorNorm_L2_H1_with_analytical_sol(const MultiLevelSolution* ml_sol, 
                                                                    const std::string solution_name,
                                                                    const Math::Function< double > * function_in,
                                                                    double    (* function_value )  (const std::vector<double> & ),
                                                                    void    (* function_gradient)  (const std::vector < double > & , std::vector < double >&  )
                                                                    ) {      
  if (function_in == NULL) {
     if (function_value == NULL) { abort(); }
     if (function_gradient == NULL) { abort(); }
  }
  else 
  {
     if (function_value != NULL) { abort(); }
     if (function_gradient != NULL) { abort(); }
  }
  
  
  
  
  const unsigned level = ml_sol->GetMLMesh()->GetNumberOfLevels() - 1u;
  
  //  extract pointers to the several objects that we are going to use
  const Mesh*     msh = ml_sol->GetMLMesh()->GetLevel(level);    // pointer to the mesh (level) object
  const elem*     el  = msh->GetMeshElements();  // pointer to the elem object in msh (level)
  const Solution* sol = ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  const unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  const unsigned soluIndex = ml_sol->GetIndex( solution_name.c_str() );    // get the position of "u" in the ml_sol object
  const unsigned soluType = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  std::vector < double >  solu; // local solution

  std::vector < std::vector < double > > x(dim);    // local coordinates
  const unsigned xType = CONTINUOUS_BIQUADRATIC; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives

  double weight = 0.; // gauss point weight

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);

  double seminorm = 0.;
  double l2norm = 0.;

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    
    const short unsigned ielGeom = msh->GetElementType(iel);
    const unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    const unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    solu.resize(nDofu);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
      const unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      const unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);  // global extraction and local storage for the element coordinates
      }
    }

    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
        
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, boost::none);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double solu_gss = 0.;
      std::vector < double > gradSolu_gss(dim, 0.);
      std::vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss += phi[i] * solu[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }

      std::vector <double> exactGradSol(dim);

      if (function_in != NULL) { exactGradSol = function_in->gradient(x_gss); }
      else {      function_gradient(x_gss, exactGradSol); }

      for (unsigned j = 0; j < dim ; j++) {
        seminorm   += ((gradSolu_gss[j] - exactGradSol[j]) * (gradSolu_gss[j] - exactGradSol[j])) * weight;
      }

      const double exactSol = (function_in) ? function_in->value(x_gss) : function_value(x_gss);
      
      l2norm += (exactSol - solu_gss) * (exactSol - solu_gss) * weight;
    } // end gauss point loop
  } //end element loop for each process

  // add the norms of all processes
  NumericVector* norm_vec;
  norm_vec = NumericVector::build().release();
  norm_vec->init(msh->n_processors(), 1 , false, AUTOMATIC);

  norm_vec->set(iproc, l2norm);
  norm_vec->close();
  l2norm = norm_vec->l1_norm();

  norm_vec->set(iproc, seminorm);
  norm_vec->close();
  seminorm = norm_vec->l1_norm();

  delete norm_vec;

  return std::pair < double, double > (sqrt(l2norm), sqrt(seminorm));

}


std::pair < double, double > GetErrorNorm_L2_H1_with_analytical_sol(const MultiLevelSolution* ml_sol, 
                                                                    const std::vector< Unknown > & unknowns_vec,
                                                                    double    (* function_value )  (const std::vector<double> & ),
                                                                    void    (* function_gradient)  (const std::vector < double > & , std::vector < double >&  )
                                                                    ){
   if (unknowns_vec.size() != 1) abort();
                                                                     
 return   GetErrorNorm_L2_H1_with_analytical_sol(ml_sol, unknowns_vec[0]._name, NULL, function_value, function_gradient) ;                                                                       
                                                                      
                                                                    }


std::pair < double, double > GetErrorNorm_L2_H1_with_analytical_sol(const MultiLevelSolution* ml_sol, 
                                                                    const std::string solution_name,
                                                                    const Math::Function< double > * function_in
                                                                    ) {
                                                                      
  return GetErrorNorm_L2_H1_with_analytical_sol(ml_sol, solution_name, function_in, NULL, NULL);
  
}                                                                    
                                                                    
                                                                    

std::pair < double, double > GetErrorNorm_L2_H1_with_analytical_sol(const MultiLevelSolution* ml_sol, 
                                                                    const std::string solution_name,
                                                                    double    (* function_value )  (const std::vector<double> & ),
                                                                    void    (* function_gradient)  (const std::vector < double > & , std::vector < double >&  )
                                                                    ) {

  return GetErrorNorm_L2_H1_with_analytical_sol(ml_sol, solution_name, NULL, function_value, function_gradient);
  
}




// ||u_i - u_h||/||u_i-u_(h/2)|| = 2^alpha, alpha is order of conv 
std::pair < double, double > GetErrorNorm_L2_H1_multiple_methods(MultiLevelSolution* ml_sol,
                                          Solution* sol_finer,
                                          std::vector< Unknown > & unknowns_vec,
                                                                    double    (* function_value )  (const std::vector<double> & ),
                                                                    void    (* function_gradient)  (const std::vector < double > & , std::vector < double >&  )
                                          ) {
  
  if (unknowns_vec.size() != 1) abort();
  
  unsigned level = ml_sol->GetMLMesh()->GetNumberOfLevels() - 1u;
  
  //  extract pointers to the several objects that we are going to use
  Mesh*     msh = ml_sol->GetMLMesh()->GetLevel(level);    // pointer to the mesh (level) object
  elem*     el  = msh->GetMeshElements();  // pointer to the elem object in msh (level)
  Solution* sol = ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned iproc = msh->processor_id(); // get the process_id (for parallel computation)

  //solution variable
  unsigned soluIndex = ml_sol->GetIndex( unknowns_vec[0]._name.c_str() );    // get the position of "u" in the ml_sol object
  unsigned soluType = ml_sol->GetSolutionType(soluIndex);    // get the finite element type for "u"

  std::vector < double >  solu; // local solution

  std::vector < std::vector < double > > x(dim);    // local coordinates
  unsigned xType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  std::vector <double> phi;  // local test function
  std::vector <double> phi_x; // local test function first order partial derivatives

  double weight; // gauss point weight

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27
  solu.reserve(maxSize);

  std::vector < double >  solu_finer;   solu_finer.reserve(maxSize);
  
  std::vector < double >  solu_exact_at_dofs;   solu_exact_at_dofs.reserve(maxSize);

  for (unsigned i = 0; i < dim; i++)
    x[i].reserve(maxSize);

  phi.reserve(maxSize);
  phi_x.reserve(maxSize * dim);

  double seminorm = 0.;
  double l2norm = 0.;
  double seminorm_exact_dofs = 0.;
  double l2norm_exact_dofs = 0.;
  double seminorm_inexact = 0.;
  double l2norm_inexact = 0.;

  // element loop: each process loops only on the elements that owns
  for (int iel = msh->GetElementOffset(iproc); iel < msh->GetElementOffset(iproc + 1); iel++) {

    
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofu  = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofx = msh->GetElementDofNumber(iel, xType);    // number of coordinate element dofs

    // resize local arrays
    solu.resize(nDofu);
    solu_finer.resize(nDofu);
    solu_exact_at_dofs.resize(nDofu);

    for (int i = 0; i < dim; i++) {
      x[i].resize(nDofx);
    }

    // local storage of coordinates
    for (unsigned i = 0; i < nDofx; i++) {
      unsigned xDof  = msh->GetSolutionDof(i, iel, xType);    // global to global mapping between coordinates node and coordinate dof

      for (unsigned jdim = 0; jdim < dim; jdim++) {
        x[jdim][i] = (*msh->GetTopology()->_Sol[jdim])(xDof);  // global extraction and local storage for the element coordinates
      }
    }
    
    const double weird_multigrid_factor = 0.25;  //don't know!
    
         std::vector <double> x_at_node(dim,0.);
      for (unsigned i = 0; i < nDofu; i++) {
         for (unsigned jdim = 0; jdim < dim; jdim++) {
             x_at_node[jdim] = x[jdim][i];
         }
      }

    // local storage of global mapping and solution
    for (unsigned i = 0; i < nDofu; i++) {
//         std::vector<double> x_at_node(dim,0.);
//         for (unsigned jdim = 0; jdim < dim; jdim++) x_at_node[jdim] = x[jdim][i];
      unsigned solDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i]       = (*sol->_Sol[soluIndex])(solDof);      // global extraction and local storage for the solution
      solu_finer[i] = /*weird_multigrid_factor * */(*sol_finer->_Sol[soluIndex])(solDof);
      solu_exact_at_dofs[i] = function_value(x_at_node);
    }


    // *** Gauss point loop ***
    for (unsigned ig = 0; ig < msh->_finiteElement[ielGeom][soluType]->GetGaussPointNumber(); ig++) {
        
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][soluType]->Jacobian(x, ig, weight, phi, phi_x, boost::none);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      double solu_gss = 0.;
      double solu_finer_gss = 0.;
      double exactSol_from_dofs_gss = 0.;
      std::vector < double > gradSolu_gss(dim, 0.);
      std::vector < double > gradSolu_exact_at_dofs_gss(dim, 0.);
      std::vector < double > gradSolu_finer_gss(dim, 0.);
      std::vector < double > x_gss(dim, 0.);

      for (unsigned i = 0; i < nDofu; i++) {
        solu_gss       += phi[i] * solu[i];
        solu_finer_gss += phi[i] * solu_finer[i];
        exactSol_from_dofs_gss += phi[i] * solu_exact_at_dofs[i];

        for (unsigned jdim = 0; jdim < dim; jdim++) {
          gradSolu_gss[jdim] += phi_x[i * dim + jdim] * solu[i];
          gradSolu_exact_at_dofs_gss[jdim] += phi_x[i * dim + jdim] * solu_exact_at_dofs[i];
          gradSolu_finer_gss[jdim] += phi_x[i * dim + jdim] * solu_finer[i];
          x_gss[jdim] += x[jdim][i] * phi[i];
        }
      }

      std::vector <double> exactGradSol(dim);
      function_gradient(x_gss, exactGradSol);

      for (unsigned j = 0; j < dim ; j++) {
        seminorm           += ((gradSolu_gss[j]       - exactGradSol[j]) * (gradSolu_gss[j]       - exactGradSol[j])) * weight;
        seminorm_inexact   += ((gradSolu_gss[j] - gradSolu_finer_gss[j]) * (gradSolu_gss[j] - gradSolu_finer_gss[j])) * weight;
        seminorm_exact_dofs      += ((gradSolu_gss[j]       - gradSolu_exact_at_dofs_gss[j]) * (gradSolu_gss[j]       - gradSolu_exact_at_dofs_gss[j])) * weight;
      }

      double exactSol = function_value(x_gss);
      l2norm              += (solu_gss - exactSol)                * (solu_gss - exactSol)       * weight;
      l2norm_exact_dofs   += (solu_gss - exactSol_from_dofs_gss)  * (solu_gss - exactSol_from_dofs_gss)       * weight;
      l2norm_inexact      += (solu_gss - solu_finer_gss)          * (solu_gss - solu_finer_gss) * weight;
   } // end gauss point loop
  } //end element loop for each process

  // add the norms of all processes
  NumericVector* norm_vec;
  norm_vec = NumericVector::build().release();
  norm_vec->init(msh->n_processors(), 1 , false, AUTOMATIC);

  norm_vec->set(iproc, l2norm);
  norm_vec->close();
  l2norm = norm_vec->l1_norm();

  norm_vec->set(iproc, seminorm);
  norm_vec->close();
  seminorm = norm_vec->l1_norm();

  delete norm_vec;

   // add the norms of all processes
  NumericVector* norm_vec_exact_using_dofs;
  norm_vec_exact_using_dofs = NumericVector::build().release();
  norm_vec_exact_using_dofs->init(msh->n_processors(), 1 , false, AUTOMATIC);

  norm_vec_exact_using_dofs->set(iproc, l2norm_exact_dofs);
  norm_vec_exact_using_dofs->close();
  l2norm_exact_dofs = norm_vec_exact_using_dofs->l1_norm();

  norm_vec_exact_using_dofs->set(iproc, seminorm_exact_dofs);
  norm_vec_exact_using_dofs->close();
  seminorm_exact_dofs = norm_vec_exact_using_dofs->l1_norm();

  delete norm_vec_exact_using_dofs;

  // add the norms of all processes
  NumericVector* norm_vec_inexact;
  norm_vec_inexact = NumericVector::build().release();
  norm_vec_inexact->init(msh->n_processors(), 1 , false, AUTOMATIC);

  norm_vec_inexact->set(iproc, l2norm_inexact);
  norm_vec_inexact->close();
  l2norm_inexact = norm_vec_inexact->l1_norm();

  norm_vec_inexact->set(iproc, seminorm_inexact);
  norm_vec_inexact->close();
  seminorm_inexact = norm_vec_inexact->l1_norm();

  delete norm_vec_inexact;

  std::pair < double, double > inexact_pair(sqrt(l2norm_inexact), sqrt(seminorm_inexact));
  
  return std::pair < double, double > (sqrt(l2norm), sqrt(seminorm));
  // return std::pair < double, double > (sqrt(l2norm_exact_dofs), sqrt(seminorm_exact_dofs));  ///@todo does not seem to work
  // return std::pair < double, double > (sqrt(l2norm_inexact), sqrt(seminorm_inexact));        ///@todo does not seem to work

}

 
 
 
 




}  //end namespace









#endif

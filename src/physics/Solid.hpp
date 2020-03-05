/*=========================================================================

 Program: FEMUS
 Module: Solid
 Authors: Simone Bnà, Giorgio Bornia

 Copyright (c) FEMTTU
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_physics_Solid_hpp__
#define __femus_physics_Solid_hpp__

//----------------------------------------------------------------------------
// includes :
//----------------------------------------------------------------------------
#include <vector>

#include "Material.hpp"


namespace femus {



//------------------------------------------------------------------------------
// Forward declarations
//------------------------------------------------------------------------------
class Parameter;

class Solid : public Material {

public:

    /** constructors */
    Solid(Parameter& par);

    Solid(Parameter& par,
          const double young_module, 
          const double poisson_coeff,
          const double density, 
          const char model[]= "Linear_elastic",
          const double k = 1., 
          const double cp = 1., 
          const double alpha = 1.e-06);

    Solid();

    /** destructor */
    ~Solid() {};

    /** To be Added */
    void set_young_module(const double young_module);

    /** To be Added */
    void set_poisson_coeff(const double poisson_coeff);

    /** To be Added */
    const double get_young_module() const;

    /** To be Added */
    const double get_poisson_coeff() const;

    /** To be Added */
    const double get_lame_lambda() const;

    /** To be Added */
    const double get_lame_shear_modulus() const;

    /** To be Added */
    const unsigned get_physical_model() const;
    
    /** To be Added */
    const bool get_if_penalty() const;
    
    /** To be Added */
    const bool get_if_mass_penalty() const;

    /** printing operator */
    friend std::ostream & operator << (std::ostream & os, const Solid & solid);

    /** overloading of the = operator */
    Solid & operator=(const Solid &solid);
    
template < class real_num_mov >
static std::vector < std::vector < real_num_mov > >  get_Cauchy_stress_tensor(const unsigned int solid_model,
                                                                              const double & mus,
                                                                              const double & lambda,
                                                                              const bool  incompressible,
                                                                              const unsigned int dim,
                                                                              const unsigned int sol_index_displ,
                                                                              const unsigned int sol_pde_index_press,
                                                                              const std::vector < std::vector < real_num_mov > > & gradSolVAR_hat_qp,
                                                                              const          std::vector < real_num_mov >   & SolVAR_qp,
                                                                              const          std::vector < unsigned int >   & SolPdeIndex,
                                                                              real_num_mov & J_hat,
                                                                              real_num_mov & trace_e_hat);
    


template < class real_num_mov >
static real_num_mov   get_mass_balance_reference_domain(const unsigned int solid_model,
                                                        const bool penalty,
                                                        const bool incompressible,
                                                        const double & lambda,
                                                        const real_num_mov & trace_e_hat,
                                                        const real_num_mov & J_hat,
                                                        const std::vector < real_num_mov >   & SolVAR_qp,
                                                        const std::vector < unsigned int >   & SolPdeIndex,
                                                        const unsigned int sol_pde_index_press);


template < class real_num_mov >
static real_num_mov   get_mass_balance_moving_domain(const std::vector < std::vector < real_num_mov > > & gradSolVAR_qp,
                                                     const std::vector < unsigned int >   & SolPdeIndex); 

    

private:

    double _young_module;

    double _poisson_coeff;

    double _lambda_lame;

    double _mu_lame;

    unsigned _model;
    
    bool _penalty;
    
    bool _mass_penalty;

};



template < class real_num_mov >
/*static*/ std::vector < std::vector < real_num_mov > >  Solid::get_Cauchy_stress_tensor(const unsigned int solid_model,
                                                             const double & mus,
                                                             const double & lambda,
                                                             const bool  incompressible,
                                                             const unsigned int dim,
                                                             const unsigned int sol_index_displ,
                                                             const unsigned int sol_pde_index_press,
                                                             const std::vector < std::vector < real_num_mov > > & gradSolVAR_hat_qp,
                                                             const          std::vector < real_num_mov >   & SolVAR_qp,
                                                             const          std::vector < unsigned int >   & SolPdeIndex,
                                                             real_num_mov & J_hat,
                                                             real_num_mov & trace_e_hat
) {
    
    
//     const unsigned int is_incompressible = 1;  //0 means compressible
          
          
          std::vector < std::vector < real_num_mov > > Cauchy(3);    for (int i = 0; i < Cauchy.size(); i++) Cauchy[i].resize(3);
          
          real_num_mov Identity[3][3] = {{ 1., 0., 0.}, { 0., 1., 0.}, { 0., 0., 1.}};

          real_num_mov I1_B = 0.;
          real_num_mov I2_B = 0.;

          if (solid_model == 0) { // Saint-Venant
            real_num_mov e[3][3];

            //computation of the stress tensor
            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                e[i][j] = 0.5 * (gradSolVAR_hat_qp[SolPdeIndex[i]][j] + gradSolVAR_hat_qp[SolPdeIndex[j]][i]);
              }
            }

          trace_e_hat = 0.;    //trace of deformation tensor
          
            for (int i = 0; i < dim; i++) {  trace_e_hat += e[i][i];   }

            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                //incompressible
                Cauchy[i][j] = 2. * mus *  e[i][j] -  (incompressible) * SolVAR_qp[SolPdeIndex[sol_pde_index_press]] * Identity[i][j];  ///@todo check that mus is multiplying everything or only the deformation tensor
                //+(penalty)*lambda*trace_e*Identity[i][j];
              }
            }
          }

          else { // hyperelastic non linear material
            real_num_mov F[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}};
            real_num_mov B[3][3];

            for (int i = 0; i < dim; i++) {
              for (int j = 0; j < dim; j++) {
                F[i][j] += gradSolVAR_hat_qp[SolPdeIndex[i]][j];
              }
            }

            J_hat  = F[0][0] * F[1][1] * F[2][2] + F[0][1] * F[1][2] * F[2][0] + F[0][2] * F[1][0] * F[2][1]
                   - F[2][0] * F[1][1] * F[0][2] - F[2][1] * F[1][2] * F[0][0] - F[2][2] * F[1][0] * F[0][1];  //Jacobian wrt fixed domain
                                

            for (int I = 0; I < 3; ++I) {
              for (int J = 0; J < 3; ++J) {
                B[I][J] = 0.;

                for (int K = 0; K < 3; ++K) {
                  //left Cauchy-Green deformation tensor or Finger tensor (b = F*F^T)
                  B[I][J] += F[I][K] * F[J][K];
                }
              }
            }

            if (solid_model <= 4) { // Neo-Hookean
              I1_B = B[0][0] + B[1][1] + B[2][2];

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  if (1  ==  solid_model) 
			Cauchy[I][J] = mus * B[I][J] - (incompressible) * mus * I1_B * SolVAR_qp[SolPdeIndex[sol_pde_index_press]] * Identity[I][J]; 	//Wood-Bonet J_hat  =1;   ///@todo check presence of mu_s here
                  else if (2  ==  solid_model) 
			Cauchy[I][J] = mus / J_hat * B[I][J] - (incompressible) * mus / J_hat * SolVAR_qp[SolPdeIndex[sol_pde_index_press]] * Identity[I][J]; //Wood-Bonet J_hat !=1;  ///@todo check presence of mu_s here
                  else if (3  ==  solid_model) 
			Cauchy[I][J] = mus * (B[I][J] - Identity[I][J]) / J_hat + lambda / J_hat * log(J_hat) * Identity[I][J]; 	//Wood-Bonet penalty
                  else if (4  ==  solid_model) 
			Cauchy[I][J] = mus * (B[I][J] - I1_B * Identity[I][J] / 3.) / pow(J_hat, 5. / 3.) + lambda * (J_hat - 1.) * Identity[I][J];  	 //Allan-Bower

                }
              }
            }
            else if (5  ==  solid_model) {  //Mooney-Rivlin
              real_num_mov detB =   	B[0][0] * (B[1][1] * B[2][2] - B[2][1] * B[1][2])
                                      - B[0][1] * (B[2][2] * B[1][0] - B[1][2] * B[2][0])
                                      + B[0][2] * (B[1][0] * B[2][1] - B[2][0] * B[1][1]);
              real_num_mov invdetB = 1. / detB;
              real_num_mov invB[3][3];

              invB[0][0] =  (B[1][1] * B[2][2] - B[2][1] * B[1][2]) * invdetB;
              invB[1][0] = -(B[0][1] * B[2][2] - B[0][2] * B[2][1]) * invdetB;
              invB[2][0] =  (B[0][1] * B[1][2] - B[0][2] * B[1][1]) * invdetB;
              invB[0][1] = -(B[1][0] * B[2][2] - B[1][2] * B[2][0]) * invdetB;
              invB[1][1] =  (B[0][0] * B[2][2] - B[0][2] * B[2][0]) * invdetB;
              invB[2][1] = -(B[0][0] * B[1][2] - B[1][0] * B[0][2]) * invdetB;
              invB[0][2] =  (B[1][0] * B[2][1] - B[2][0] * B[1][1]) * invdetB;
              invB[1][2] = -(B[0][0] * B[2][1] - B[2][0] * B[0][1]) * invdetB;
              invB[2][2] =  (B[0][0] * B[1][1] - B[1][0] * B[0][1]) * invdetB;

              I1_B = 	B[0][0] + B[1][1] + B[2][2];
              I2_B = 	B[0][0] * B[1][1] + B[1][1] * B[2][2] + B[2][2] * B[0][0]
                      - B[0][1] * B[1][0] - B[1][2] * B[2][1] - B[2][0] * B[0][2];

              double C1 = mus / 3.;
              double C2 = C1 / 2.;

              for (int I = 0; I < 3; ++I) {
                for (int J = 0; J < 3; ++J) {
                  Cauchy[I][J] =  2.*(C1 * B[I][J] - C2 * invB[I][J])
                                  //- (2. / 3.) * (C1 * I1_B - C2 * I2_B) * SolVAR_qp[SolPdeIndex[sol_pde_index_press]] * Identity[I][J];
                                  - (incompressible) * SolVAR_qp[SolPdeIndex[sol_pde_index_press]] * Identity[I][J];
               }
              }

            }
          }

          
          return Cauchy;
}


template < class real_num_mov >
/*static*/ real_num_mov   Solid::get_mass_balance_reference_domain(const unsigned int solid_model, 
                                const bool penalty,
                                const bool incompressible,
                                const double & lambda,
                                const real_num_mov & trace_e_hat,
                                const real_num_mov & J_hat,
                                const std::vector < real_num_mov >   & SolVAR_qp,
                                const std::vector < unsigned int >   & SolPdeIndex,
                                const unsigned int sol_pde_index_press) {
    
  real_num_mov  mass_balance = 0.;
  
              if (!penalty) {
                     if (0  ==  solid_model)                           mass_balance = trace_e_hat;
                else if (1  ==  solid_model || 5  ==  solid_model)     mass_balance = J_hat - 1.         + (!incompressible) / lambda * SolVAR_qp[SolPdeIndex[sol_pde_index_press]];
                else if (2  ==  solid_model)                           mass_balance = log(J_hat) / J_hat + (!incompressible) / lambda * SolVAR_qp[SolPdeIndex[sol_pde_index_press]];
              }
                else if (3  ==  solid_model || 4  ==  solid_model)     mass_balance = SolVAR_qp[SolPdeIndex[sol_pde_index_press]] ; // pressure = 0 in the solid
              
 return mass_balance;              
              
}


template < class real_num_mov >
/*static*/ real_num_mov   Solid::get_mass_balance_moving_domain(const std::vector < std::vector < real_num_mov > > & gradSolVAR_qp,
                                                     const std::vector < unsigned int >   & SolPdeIndex) {
    
 real_num_mov  mass_balance = 0.;
 
 real_num_mov div_displ = 0.;  for (int idim = 0; idim < gradSolVAR_qp[SolPdeIndex[idim]].size(); idim++) div_displ += gradSolVAR_qp[SolPdeIndex[idim]][idim];

 mass_balance = div_displ;
 
 
 return mass_balance;              
    
}

    



} //end namespace femus



#endif

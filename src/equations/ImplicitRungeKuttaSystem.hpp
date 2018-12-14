/*=========================================================================

 Program: FEMuS
 Module: ImplicitRungeKuttaSystem
 Authors: Eugenio Aulisa and Erdi Kara

 Copyright (c) FEMuS
 All rights reserved.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __femus_equations_ImplicitRungeKuttaSystem_hpp__
#define __femus_equations_ImplicitRungeKuttaSystem_hpp__

//------------------------------------------------------------------------------
// includes
//------------------------------------------------------------------------------
#include "TransientSystem.hpp"
#include <string>
#include <vector>
#include "assert.h"

namespace femus {
    
/**
 * This class provides a specific system class for the time integration of system PDE
 * using the Runge Kutta schemes.
 */

// ------------------------------------------------------------
// TransientSystem class definition
template <class Base> 
class ImplicitRungeKuttaSystem : public TransientSystem<Base> { 

public:

    /** Constructor.  Initializes required data structures. */
    ImplicitRungeKuttaSystem (MultiLevelProblem& ml_probl,
                            const std::string& name,
                            const unsigned int number,
                            const MgSmoother & smoother_type);

    /** Destructor. */
    virtual ~ImplicitRungeKuttaSystem ();

    virtual void clear ();
    
    void AddSolutionToSystemPDE(const char solname[]);
    
    inline void SetRungeKuttaStages(const unsigned & RK){
      _RK = RK;
    }
    
    inline unsigned GetRungeKuttaStages(){
      return _RK;
    }
    
    void UpdateSolution();
    
    /** calling the parent solve */
    void MLsolve();

    /** calling the parent solve */
    void MGsolve( const MgSmootherType& mgSmootherType = MULTIPLICATIVE );
    
    const std::vector < std::ostringstream > & GetSolkiNames(const char name[]);
    
    void GetIntermediateSolutions(const std::vector < double > &soluOld, 
        const std::vector < std::vector < adept::adouble > > & solk, 
        std::vector < std::vector < adept::adouble > > & solu);
    
private:
   unsigned _RK;
   static const double _b[5][5];
   static const double _a[5][5][5];
   std::vector < std::ostringstream > _solName;
   std::vector < std::vector < std::ostringstream > > _solKiName;
};


template <class Base>
const double ImplicitRungeKuttaSystem<Base>::_b[5][5]={
   {1.},
   {0.5,0.5},
   {5./18., 4./9., 5./18.}
}; 
   
template <class Base>
const double ImplicitRungeKuttaSystem<Base>::_a[5][5][5]={
    { {0.5}},
    {
      {0.25, 0.25 - sqrt(3.)/6.},
      {0.25 + sqrt(3.)/6., 0.25 }
    },
    {
      { 5./36.,                   2./9. - sqrt(15.)/15.,   5./36. - sqrt(15.)/30.},
      { 5./36. + sqrt(15.)/24.,   2./9.,                   5./36. - sqrt(15.)/24.},
      { 5./36. + sqrt(15.)/30.,   2./9. + sqrt(15.)/15.,   5./36.}
    }
};   


template <class Base>
ImplicitRungeKuttaSystem<Base>::ImplicitRungeKuttaSystem(
            MultiLevelProblem& ml_probl,
            const std::string& name,
            const unsigned int number, 
            const MgSmoother & smoother_type):
            TransientSystem<Base>(ml_probl, name, number, smoother_type),
		    _RK(1)
{
 
}

/** Destructor. */
template <class Base>
ImplicitRungeKuttaSystem<Base>::~ImplicitRungeKuttaSystem()
{
  this->clear();
}

template <class Base>
void ImplicitRungeKuttaSystem<Base>::clear()
{
  _RK = 1;  
  _solName.resize(0);
  _solKiName.resize(0);
  TransientSystem<Base>::clear();
}

template <class Base>
void ImplicitRungeKuttaSystem<Base>::AddSolutionToSystemPDE(const char solname[]){
    
    unsigned size = _solName.size();
   _solName.resize(size+1);
   _solName[size] << solname;    
   
   _solKiName.resize(size+1);
   _solKiName[size].resize(_RK);
   
  for(unsigned i = 0; i < _RK; i++){
    
    _solKiName[size][i] << solname << "k" << i+1;
    
    unsigned solIndex = this->_ml_sol->GetIndex(solname);
        
    this->_ml_sol->AddSolution(_solKiName[size][i].str().c_str(), this->_ml_sol->GetSolutionFamily(solIndex), this->_ml_sol->GetSolutionOrder(solIndex)); 
    
    this->_ml_sol->Initialize(_solKiName[size][i].str().c_str()); // Since this is time depend problem.
    
    this->_ml_sol->GenerateBdc(_solKiName[size][i].str().c_str());
     
    this->Base::AddSolutionToSystemPDE(_solKiName[size][i].str().c_str());
  }
  
}


template <class Base>
  void ImplicitRungeKuttaSystem<Base>::MLsolve(){
    TransientSystem<Base>::MLsolve(); 
    UpdateSolution();
}

template <class Base>
  void ImplicitRungeKuttaSystem<Base>::MGsolve(const MgSmootherType& mgSmootherType){
    TransientSystem<Base>::MGsolve(mgSmootherType);  
    UpdateSolution();    
  }

template <class Base>
  void ImplicitRungeKuttaSystem<Base>::UpdateSolution(){
    unsigned level = this->_solution.size() - 1u;  
    
    for(unsigned i = 0; i < _solName.size(); i++){
      unsigned solIndex = this->_ml_sol->GetIndex(_solName[i].str().c_str());  
      for( unsigned j = 0; j < _RK; j++ ){
         
         unsigned solkiIndex = this->_ml_sol->GetIndex( _solKiName[i][j].str().c_str() );   
                 
         this->_solution[level]->_Sol[solIndex]-> add( _b[_RK - 1][j] * this->_dt, *(this->_solution[level]->_Sol[solkiIndex]));
      }
      this->_solution[level]->_Sol[solIndex]->close();
    }
  }
  
template <class Base>
  const std::vector < std::ostringstream > & ImplicitRungeKuttaSystem<Base>::GetSolkiNames(const char name[]){
      
    unsigned index = 0;

    while(strcmp(_solName[index].str().c_str(), name)) {
      index++;

      if(index == _solName.size()) {
       std:: cout << "error! invalid solution name: in entry ImplicitRungeKuttaSystem.GetSolkiNames(\"" << name << "\")" << std::endl;
        abort();
      }
    }

    return _solKiName[index];
  }  
      
template <class Base>
  void ImplicitRungeKuttaSystem<Base>::GetIntermediateSolutions(const std::vector < double > &soluOld, 
    const std::vector < std::vector < adept::adouble > > & solk, 
    std::vector < std::vector < adept::adouble > > & solu){
      
    //local storage of global mapping and solution
    for (unsigned i = 0; i < soluOld.size(); i++) {
      for( unsigned j = 0; j < _RK; j++ ){
        solu[j][i] = soluOld[i];
        for( unsigned k = 0; k < _RK; k++ ){
          solu[j][i] += this->_dt * _a[_RK-1][j][k]  * solk[k][i]; // global extraction and local storage for the solution
        }
      }
    }
  }      

// -----------------------------------------------------------
// Useful typedefs
typedef ImplicitRungeKuttaSystem<LinearImplicitSystem> ImplicitRungeKuttaLinearImplicitSystem;
typedef ImplicitRungeKuttaSystem<NonLinearImplicitSystem> ImplicitRungeKuttaNonlinearImplicitSystem;
typedef ImplicitRungeKuttaSystem<MonolithicFSINonLinearImplicitSystem> ImplicitRungeKuttaMonolithicFSINonlinearImplicitSystem;

} //end namespace femus



#endif







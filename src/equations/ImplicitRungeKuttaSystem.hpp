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

#include "ImplicitRKEnum.hpp"

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
                                const LinearEquationSolverType & smoother_type);

      /** Destructor. */
      virtual ~ImplicitRungeKuttaSystem ();

      void AddSolutionToSystemPDE (const char solname[]);

      void SetImplicitRungeKuttaScheme(const ImplicitRKScheme & RKtype);
      
      ImplicitRKScheme GetImplicitRungeKuttaScheme(){
        return _RKScheme;
      }
      

      inline unsigned GetRungeKuttaStages() {
        return _RK;
      }

      void UpdateSolution();

      /** calling the parent solve */
      void MGsolve (const MgSmootherType& mgSmootherType = MULTIPLICATIVE);

      const std::vector < std::string > & GetSolkiNames (const char name[]);

      void GetIntermediateSolutions (const std::vector < double > &soluOld,
                                     const std::vector < std::vector < adept::adouble > > & solk,
                                     std::vector < std::vector < adept::adouble > > & solu);

      void SetIntermediateTimes();
      const std::vector < double > & GetIntermediateTimes();
      void SetRKVariableType (const char solname[], const bool &type);
      
      void GetbAi();
      
    private:
      unsigned _RK;
      ImplicitRKScheme _RKScheme;

      const double *_A, *_Ai, *_b, *_c;
      std::vector< double > _bAi;
      double _bAiSum;
      
      
      std::vector < std::string > _solName;
      std::vector < unsigned > _solIndex;
      std::vector < bool > _solRKType;
      std::vector < std::vector < std::string > > _solKiName;
      std::vector < std::vector < unsigned > > _solKiIndex;
      std::vector < double > _itime;
      double _time0;
  };

  template <class Base>
  ImplicitRungeKuttaSystem<Base>::ImplicitRungeKuttaSystem (
    MultiLevelProblem& ml_probl,
    const std::string& name,
    const unsigned int number,
    const LinearEquationSolverType & smoother_type) :
    TransientSystem<Base> (ml_probl, name, number, smoother_type),
    _RK (1),
    _RKScheme(LEGENDRE1),
    _c (cLEGENDRE1),
    _b (bLEGENDRE1),
    _A (aLEGENDRE1),
    _Ai(aiLEGENDRE1)
    {
      GetbAi();
  }

  /** Destructor. */
  template <class Base>
  ImplicitRungeKuttaSystem<Base>::~ImplicitRungeKuttaSystem() {
    _RK = 1;
    _solName.resize (0);
    _solKiName.resize (0);
    _solIndex.resize (0);
    _solKiIndex.resize (0);
    _solRKType.resize(0);
    
  }

  template <class Base>
  void ImplicitRungeKuttaSystem<Base>::AddSolutionToSystemPDE (const char solname[]) {

    std::ostringstream soliname;
    soliname << solname; 
    unsigned size = _solName.size();
    _solName.resize (size + 1u);
    _solName[size] = soliname.str();
   
    _solIndex.resize (size + 1u);
    _solIndex[size] = this->_ml_sol->GetIndex (solname);

    _solKiName.resize (size + 1u);
    _solKiName[size].resize (_RK);

    _solKiIndex.resize (size + 1u);
    _solKiIndex[size].resize (_RK);

    _solRKType.resize(size + 1u);
    _solRKType[size] = true;
    
    _itime.assign (_RK, 0.);
    _time0 = 0.;

    for (unsigned i = 0; i < _RK; i++) {

      std::ostringstream solkiname; 
      solkiname << solname << "k" << i + 1; 
      _solKiName[size][i] = solkiname.str();

      this->_ml_sol->AddSolution (_solKiName[size][i].c_str(), this->_ml_sol->GetSolutionFamily (_solIndex[size]), this->_ml_sol->GetSolutionOrder (_solIndex[size]));

      this->_ml_sol->Initialize (_solKiName[size][i].c_str()); // Since this is time depend problem.

      _solKiIndex[size][i] = this->_ml_sol->GetIndex (_solKiName[size][i].c_str());

      this->Base::AddSolutionToSystemPDE (_solKiName[size][i].c_str());
    }

    this->_ml_sol->GenerateRKBdc (_solIndex[size], _solKiIndex[size], 0, _itime, _time0, 1., _Ai);
  }

  template <class Base>
  void ImplicitRungeKuttaSystem<Base>::SetRKVariableType (const char solname[], const bool &type) {
    
    unsigned index = 0;
    
    while (strcmp (_solName[index].c_str(), solname)) {
      index++;
      
      if (index == _solName.size()) {
        std:: cout << "error! invalid solution name: in entry ImplicitRungeKuttaSystem.SetRKVariableType(\" ... \")" << std::endl;
        abort();
      }
    }
    _solRKType[index] = type;
    
  }

  template <class Base>
  void ImplicitRungeKuttaSystem<Base>::MGsolve (const MgSmootherType& mgSmootherType) {
    TransientSystem<Base>::SetUpForSolve();

    SetIntermediateTimes();

    std::cout << std::endl;
    for (unsigned i = 0; i < _solIndex.size(); i++) {
      if (!strcmp (this->_ml_sol->GetBdcType (_solIndex[i]), "Time_dependent")) {
        this->_ml_sol->GenerateRKBdc (_solIndex[i], _solKiIndex[i], 0, _itime, _time0, this->_dt, _Ai);
      }
    }

    Base::MGsolve (mgSmootherType);

    UpdateSolution();
  }

  template <class Base>
  void ImplicitRungeKuttaSystem<Base>::UpdateSolution() {
    unsigned level = this->_solution.size() - 1u;

    for (unsigned i = 0; i < _solName.size(); i++) {
      
      unsigned solIndex = _solIndex[i];
      if( _solRKType[i] ) {
        * (this->_solution[level]->_Sol[solIndex]) = * (this->_solution[level]->_SolOld[solIndex]);
        for (unsigned j = 0; j < _RK; j++) {
          unsigned solkiIndex = _solKiIndex[i][j]; 
          this->_solution[level]->_Sol[solIndex]-> add (_b[j] * this->_dt, * (this->_solution[level]->_Sol[solkiIndex]));
        }
        this->_solution[level]->_Sol[solIndex]->close();
      }
      else if( this->_ml_sol->GetSolutionTimeOrder(solIndex) == 2 ){
        * (this->_solution[level]->_Sol[solIndex]) = * (this->_solution[level]->_SolOld[solIndex]);
        this->_solution[level]->_Sol[solIndex]->scale (1. - _bAiSum);
        for (unsigned j = 0; j < _RK; j++) {
          unsigned solkiIndex = _solKiIndex[i][j]; 
          this->_solution[level]->_Sol[solIndex]-> add (_bAi[j], * (this->_solution[level]->_Sol[solkiIndex]));
        }
      }
      else{       
        unsigned solkiIndex = _solKiIndex[i][_RK-1]; 
        * (this->_solution[level]->_Sol[solIndex]) = * (this->_solution[level]->_Sol[solkiIndex]);
      }
    }
  }

  template <class Base>
  const std::vector < std::string > & ImplicitRungeKuttaSystem<Base>::GetSolkiNames (const char name[]) {

    unsigned index = 0;

    while (strcmp (_solName[index].c_str(), name)) {
      index++;

      if (index == _solName.size()) {
        std:: cout << "error! invalid solution name: in entry ImplicitRungeKuttaSystem.GetSolkiNames(\"" << name << "\")" << std::endl;
        abort();
      }
    }

    return _solKiName[index];
  }

  template <class Base>
  void ImplicitRungeKuttaSystem<Base>::GetIntermediateSolutions (
      const std::vector < double > &soluOld,
      const std::vector < std::vector < adept::adouble > > & solk,
      std::vector < std::vector < adept::adouble > > & solu) {

    //local storage of global mapping and solution
    for (unsigned i = 0; i < soluOld.size(); i++) {
      for (unsigned j = 0; j < _RK; j++) {
        solu[j][i] = soluOld[i];
        for (unsigned k = 0; k < _RK; k++) {
          solu[j][i] += this->_dt * _A[ j * _RK + k]  * solk[k][i]; // global extraction and local storage for the solution
        }
      }
    }
  }

  template <class Base>
  void ImplicitRungeKuttaSystem<Base>::SetIntermediateTimes() {

    _itime.resize (_RK);

    _time0 = this->_time - this->_dt;

    for (unsigned i = 0; i < _RK; i++) {
      _itime[i] = _time0 + this->_dt * _c[i];
    }
  }
  template <class Base>
  const std::vector < double > & ImplicitRungeKuttaSystem<Base>::GetIntermediateTimes() {
    return _itime;
  }

  template <class Base>
  void ImplicitRungeKuttaSystem<Base>::SetImplicitRungeKuttaScheme(const ImplicitRKScheme & RKscheme){
    _RKScheme = RKscheme;
    _RK = nRK[RKscheme]; 
    _c  = cIRK[RKscheme];
    _b  = bIRK[RKscheme];
    _A  = aIRK[RKscheme];
    _Ai = aiIRK[RKscheme];
    GetbAi();
  }
  
  // This idea is used to update the algebraic variable using RK in DAE problems.
  // Refer to http://www.ams.org/journals/mcom/2006-75-254/S0025-5718-05-01809-0/S0025-5718-05-01809-0.pdfS
  template <class Base>
  void ImplicitRungeKuttaSystem<Base>::GetbAi(){
    _bAi.assign(_RK, 0.);
    _bAiSum = 0.;
    for (unsigned i = 0; i <_RK; i++ ){
      for (unsigned j = 0; j <_RK; j++ ){
        _bAi[i] += _b[i] * _Ai[ i *_RK + j];
      }
      _bAiSum += _bAi[i];
    }
  }
  
// -----------------------------------------------------------
// Useful typedefs
  typedef ImplicitRungeKuttaSystem<LinearImplicitSystem> ImplicitRungeKuttaLinearImplicitSystem;
  typedef ImplicitRungeKuttaSystem<NonLinearImplicitSystem> ImplicitRungeKuttaNonlinearImplicitSystem;
  typedef ImplicitRungeKuttaSystem<MonolithicFSINonLinearImplicitSystem> ImplicitRungeKuttaMonolithicFSINonlinearImplicitSystem;

} //end namespace femus



#endif







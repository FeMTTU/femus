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

      void AddSolutionToSystemPDE (const char solname[]);

      inline void SetRungeKuttaStages (const unsigned & RK) {
        _RK = RK;
      }

      inline unsigned GetRungeKuttaStages() {
        return _RK;
      }

      void UpdateSolution();

      /** calling the parent solve */
      void MLsolve();

      /** calling the parent solve */
      void MGsolve (const MgSmootherType& mgSmootherType = MULTIPLICATIVE);

      const std::vector < std::string > & GetSolkiNames (const char name[]);

      void GetIntermediateSolutions (const std::vector < double > &soluOld,
                                     const std::vector < std::vector < adept::adouble > > & solk,
                                     std::vector < std::vector < adept::adouble > > & solu);

      void SetIntermediateTimes();
      const std::vector < double > & GetIntermediateTimes();
      void SetRKVariableType (const char solname[], const bool &type);
      
    private:
      unsigned _RK;

      static const double _a[5][5][5], _aI[5][5][5];
      static const double _b[5][5];
      static const double _c[5][5];
      std::vector < std::string > _solName;
      std::vector < unsigned > _solIndex;
      std::vector < bool > _solRKType;
      std::vector < std::vector < std::string > > _solKiName;
      std::vector < std::vector < unsigned > > _solKiIndex;
      std::vector < double > _itime;
      double _time0;
  };

  template <class Base>
  const double ImplicitRungeKuttaSystem<Base>::_c[5][5] = {
    {0.5},
    {0.5 - sqrt (3.) / 6., 0.5 + sqrt (3) / 6.},
    {0.5 - sqrt (15.) / 10., 0.5, 0.5 + sqrt (15.) / 10.},
    {1.06858, .5, -0.06858}
  };


  template <class Base>
  const double ImplicitRungeKuttaSystem<Base>::_b[5][5] = {
    {1.},
    {0.5, 0.5},
    {5. / 18., 4. / 9., 5. / 18.},
    {0.1288859568,.7422280863, .1288859568}
  };

  template <class Base>
  const double ImplicitRungeKuttaSystem<Base>::_a[5][5][5] = {
    { {0.5}},
    {
      {0.25, 0.25 - sqrt (3.) / 6.},
      {0.25 + sqrt (3.) / 6., 0.25 }
    },
    {
      { 5. / 36.,                   2. / 9. - sqrt (15.) / 15.,   5. / 36. - sqrt (15.) / 30.},
      { 5. / 36. + sqrt (15.) / 24.,   2. / 9.,                   5. / 36. - sqrt (15.) / 24.},
      { 5. / 36. + sqrt (15.) / 30.,   2. / 9. + sqrt (15.) / 15.,   5. / 36.}
    },
    {
        {1.06858, 0.0,0.},
        {-.56858, 1.06858,0.},
        {2.13716, -3.27432, 1.06858}, // Nørsett's three-stage, 4th order Diagonally IRK
    }
    
    
  };

  template <class Base>
  const double ImplicitRungeKuttaSystem<Base>::_aI[5][5][5] = {
    {
      {2.}
    },
    {
      {3., 0.4641016151377544},
      {-6.464101615137755, 3.},

    },
    {
      {5., 1.1639777949432233, -0.16397779494322232},
      {-5.727486121839513, 2.,  0.7274861218395138},
      {10.163977794943225, -9.163977794943223, 5.}
    },
    {
        {.9358213704, 0.,0.},
        {.4979405519, .9358213704, 0.},
        {-.3458639439, 2.867523845, .9358213704},
    }
  };


  template <class Base>
  ImplicitRungeKuttaSystem<Base>::ImplicitRungeKuttaSystem (
    MultiLevelProblem& ml_probl,
    const std::string& name,
    const unsigned int number,
    const MgSmoother & smoother_type) :
    TransientSystem<Base> (ml_probl, name, number, smoother_type),
    _RK (1) {

  }

  /** Destructor. */
  template <class Base>
  ImplicitRungeKuttaSystem<Base>::~ImplicitRungeKuttaSystem() {
    this->clear();
  }

  template <class Base>
  void ImplicitRungeKuttaSystem<Base>::clear() {
    _RK = 1;
    _solName.resize (0);
    _solKiName.resize (0);
    _solIndex.resize (0);
    _solKiIndex.resize (0);
    _solRKType.resize(0);
    TransientSystem<Base>::clear();
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

    this->_ml_sol->GenerateRKBdc (_solIndex[size], _solKiIndex[size], 0, _itime, _time0, 1., _aI[_RK - 1]);
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
  void ImplicitRungeKuttaSystem<Base>::MLsolve() {

    TransientSystem<Base>::SetUpForSolve();

    SetIntermediateTimes();


    for (unsigned i = 0; i < _solIndex.size(); i++) {
      if (!strcmp (this->_ml_sol->GetBdcType (_solIndex[i]), "Time_dependent")) {
        this->_ml_sol->GenerateRKBdc (_solIndex[i], _solKiIndex[i], 0, _itime, _time0, this->_dt, _aI[_RK - 1]);
      }
    }

    // call the parent MLsolver
    Base::_MLsolver = true;
    Base::_MGsolver = false;

    Base::solve();

    UpdateSolution();
  }

  template <class Base>
  void ImplicitRungeKuttaSystem<Base>::MGsolve (const MgSmootherType& mgSmootherType) {
    TransientSystem<Base>::SetUpForSolve();

    SetIntermediateTimes();

    std::cout << std::endl;
    for (unsigned i = 0; i < _solIndex.size(); i++) {
      if (!strcmp (this->_ml_sol->GetBdcType (_solIndex[i]), "Time_dependent")) {
        this->_ml_sol->GenerateRKBdc (_solIndex[i], _solKiIndex[i], 0, _itime, _time0, this->_dt, _aI[_RK - 1]);
      }
    }

    // call the parent MLsolver
    Base::_MLsolver = false;
    Base::_MGsolver = true;

    Base::solve (mgSmootherType);

    UpdateSolution();
  }

  template <class Base>
  void ImplicitRungeKuttaSystem<Base>::UpdateSolution() {
    unsigned level = this->_solution.size() - 1u;

    for (unsigned i = 0; i < _solName.size(); i++) {
      
      unsigned solIndex = _solIndex[i];//this->_ml_sol->GetIndex (_solName[i].str().c_str());
      if( _solRKType[i] ) {
        
        * (this->_solution[level]->_Sol[solIndex]) = * (this->_solution[level]->_SolOld[solIndex]);

        for (unsigned j = 0; j < _RK; j++) {

          unsigned solkiIndex = _solKiIndex[i][j]; // this->_ml_sol->GetIndex (_solKiName[i][j].str().c_str());

          this->_solution[level]->_Sol[solIndex]-> add (_b[_RK - 1][j] * this->_dt, * (this->_solution[level]->_Sol[solkiIndex]));
        }
        this->_solution[level]->_Sol[solIndex]->close();
      }
      else{
        unsigned solkiIndex = _solKiIndex[i][_RK-1]; //this->_ml_sol->GetIndex (_solKiName[i][_RK - 1].str().c_str());
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
          solu[j][i] += this->_dt * _a[_RK - 1][j][k]  * solk[k][i]; // global extraction and local storage for the solution
        }
      }
    }
  }

  template <class Base>
  void ImplicitRungeKuttaSystem<Base>::SetIntermediateTimes() {

    _itime.resize (_RK);

    _time0 = this->_time - this->_dt;

    for (unsigned i = 0; i < _RK; i++) {
      _itime[i] = _time0 + this->_dt * _c[_RK - 1][i];
    }
  }
  template <class Base>
  const std::vector < double > & ImplicitRungeKuttaSystem<Base>::GetIntermediateTimes() {
    return _itime;
  }

// -----------------------------------------------------------
// Useful typedefs
  typedef ImplicitRungeKuttaSystem<LinearImplicitSystem> ImplicitRungeKuttaLinearImplicitSystem;
  typedef ImplicitRungeKuttaSystem<NonLinearImplicitSystem> ImplicitRungeKuttaNonlinearImplicitSystem;
  typedef ImplicitRungeKuttaSystem<MonolithicFSINonLinearImplicitSystem> ImplicitRungeKuttaMonolithicFSINonlinearImplicitSystem;

} //end namespace femus



#endif







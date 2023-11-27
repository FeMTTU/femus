#ifndef __femus_LinearImplicitSystemForSSC_hpp__
#define __femus_LinearImplicitSystemForSSC_hpp__



#include "LinearImplicitSystem.hpp"

#include "MultiLevelProblem.hpp"

#include <string>


namespace femus {

    
  enum SSCType {
    SYMMETRIC1111 = 0,
    SYMMETRIC1001,
    SYMMETRIC0110,
    ASYMMETRIC0101,
    ASYMMETRIC1010
  };   
    
  
  
  class LinearImplicitSystemForSSC : public LinearImplicitSystem {
      
      
  public:

    LinearImplicitSystemForSSC (MultiLevelProblem& ml_probl, 
                                const std::string& name,
                                const unsigned int number, 
                                const LinearEquationSolverType & smoother_type) :
                                LinearImplicitSystem(ml_probl, name, number, smoother_type),
                                _factorTest( false ),
                                _scale( 1.0 ),
                                _sscLevelSmoother(true),
                                _sscType{1, 1, 1, 1}
                                {
     std::cout << std::endl;
     std::cout << std::endl;
     std::cout << "!!!!!!!!!!!!!!!!" << std::endl;
     std::cout << "Need to see from where I can initialize _gridr, the number of uniformly refined levels!!!!!!!!!!!!!!!!" << std::endl;                             
                                  
     std::cout << "Also, review the MGsolve function with the specific implementations below!!!!!!!!!!!!!!!!" << std::endl;                             
     std::cout << std::endl;
     std::cout << std::endl;
                                };
  

      void SetFactorAndScale(const bool &factorTest, const double &scale) {
        _factorTest = factorTest;
        _scale = scale;
      };

      void SetSscLevelSmoother(const bool &value) {
        _sscLevelSmoother = value;
      };

      void SetSSCType(const SSCType &sscType){
	if (sscType == SYMMETRIC1111 ){
	   _sscType[0]=1;
	   _sscType[1]=1;
	   _sscType[2]=1;
	   _sscType[3]=1;
	}
	else if (sscType == SYMMETRIC1001 ){
	   _sscType[0]=1;
	   _sscType[1]=0;
	   _sscType[2]=0;
	   _sscType[3]=1;
	}
	else if (sscType == SYMMETRIC0110 ){
	   _sscType[0]=0;
	   _sscType[1]=1;
	   _sscType[2]=1;
	   _sscType[3]=0;
	}
	else if (sscType == ASYMMETRIC0101 ){
	   _sscType[0]=0;
	   _sscType[1]=1;
	   _sscType[2]=0;
	   _sscType[3]=1;
	}
	else if (sscType == ASYMMETRIC1010 ){
	   _sscType[0]=1;
	   _sscType[1]=0;
	   _sscType[2]=1;
	   _sscType[3]=0;
	}
	else{
	  std::cout << "wrong SSC type"<<std::endl;
	  abort();
	}
      }
      


  // ********************************************

// The following routine is intended to replace this piece of code in the MGsolve(...) function in the class father:
// // // 
// // // 
// // // 
// // //  ======= BEGIN ==========================
// // // _LinSolver[igridn]->MGInit(mgSmootherType, igridn + 1, _mgOuterSolver);
// // // 
// // // for(unsigned i = 0; i < igridn + 1; i++) {
// // //   unsigned npre = (i == 0) ? _npre0 : _npre;
// // //   unsigned npost = (i == 0) ? 0 : _npost;
// // //   if(_RR[i])
// // //     _LinSolver[i]->MGSetLevel(_LinSolver[igridn], igridn, _VariablesToBeSolvedIndex, _PP[i], _RR[i], npre, npost);
// // //   else
// // //     _LinSolver[i]->MGSetLevel(_LinSolver[igridn], igridn, _VariablesToBeSolvedIndex, _PP[i], _PP[i], npre, npost);
// // // }
// // // 
// // // Vcycle(igridn, mgSmootherType);
// // // 
// // // _LinSolver[igridn]->MGClear();
// // //  ======== END =========================
// // // 
// // // 
// // //   Replace the WHOLE code above (not just the Vcycle call!) with just this call:
// // // 
// // //   VcycleForSSC( igridn );
// // // 
 
  
  bool /*LinearImplicitSystemForSSC::*/VcycleForSSC( const unsigned& gridn )
  {

    clock_t start_mg_time = clock();


    _LinSolver[gridn - 1u]->SetEpsZero();

    bool linearIsConverged;

    for ( unsigned linearIterator = 0; linearIterator < _n_max_linear_iterations; linearIterator++ ) { //linear cycle

      std::cout << std::endl << " *************** Linear iteration " << linearIterator + 1 << " ***********" << std::endl;

      bool ksp_clean = !linearIterator * _assembleMatrix;

      _solution[gridn - 1]->UpdateRes( _SolSystemPdeIndex, _LinSolver[gridn - 1]->_RES, _LinSolver[gridn - 1]->KKoffset );
      linearIsConverged = HasLinearConverged( gridn - 1 );

      for ( unsigned ig = gridn - 1u; ig > 0; ig-- ) {


        unsigned factor = 1;
        //BEGIN da commentare
        if (_factorTest) {
          if (_gridr == 1) {
            if (ig >= _gridr ) {
              factor = (1 + ig - _gridr );
            }
          }
          else {
            if (ig >= _gridr ) {
              factor = (2 + ig - _gridr); //da controllare qnd _gridr lo mettiamo >1
            }
          }
          if (_sscLevelSmoother == true) {
            factor *= factor;
            factor = factor + 1;
          }
          else {
            factor = pow(2, gridn - _gridr -  factor);
          }
        }
        //END da commentare
        //std::cout<< "ig = "<< ig << " factor =  " <<factor<<std::endl;
        std::cout << "pre-smoothing level = " << ig << " 1 + level^2 = " << factor << std::endl;

        // ============== Presmoothing ==============
        for ( unsigned k = 0; k < _npre * factor; k++ ) {
          if (!_sscLevelSmoother) { //JK and BPWX
            _LinSolver[ig]->SetRichardsonScaleFactor(_scale);
            _LinSolver[ig]->Solve( _VariablesToBeSolvedIndex, ksp_clean * ( !k ) );
          }
          else//our
            SolveForSSC(ig + 1, ksp_clean * ( !k ), _sscType[0], _sscType[1]);
        }

        // ============== Restriction ==============
        Restrictor( ig );
      }

      // ============== Direct Solver ==============
      _LinSolver[0]->Solve( _VariablesToBeSolvedIndex, ksp_clean );

      for ( unsigned ig = 1; ig < gridn; ig++ ) {
        // ============== Standard Prolongation ==============
        Prolongator( ig );

        unsigned factor = 1;
        if (_factorTest) {
          //BEGIN da commentare
          if (_gridr == 1) {
            if (ig >= _gridr ) {
              factor = (1 + ig - _gridr );
            }
          }
          else {
            if (ig >= _gridr ) {
              factor = (2 + ig - _gridr);
            }
          }
          if (_sscLevelSmoother == true) {
            factor *= factor;
            factor = factor + 1;
          }
          else {
            factor = pow(2, gridn - _gridr - factor);
          }
        }
        //END da commentare
        //std::cout<< "ig = "<< ig << " factor =  " <<factor<<std::endl;

        std::cout << "post-smoothing level = " << ig << " 1 + level^2 = " << factor << std::endl;
        // ============== PostSmoothing ==============
        for ( unsigned k = 0; k < _npost * factor; k++ ) {
          if (!_sscLevelSmoother) { //JK and BPWX
            _LinSolver[ig]->SetRichardsonScaleFactor(_scale);
            _LinSolver[ig]->Solve( _VariablesToBeSolvedIndex, ksp_clean * ( !_npre ) * ( !k ) );
          }
          else { //our
	    //Solve(ig + 1, ksp_clean * ( !(_npre * factor) ) * ( !k ), 1, 1);
            SolveForSSC(ig + 1, ksp_clean * ( !(_npre * factor) ) * ( !k ), _sscType[2], _sscType[3]);
          }

        }
      }

      // ============== Update AMR Solution and Residual ( _gridr-1 <= ig <= gridn-2 ) ==============
      for ( unsigned ig = _gridr - 1; ig < gridn - 1; ig++ ) { // _gridr
        _solution[ig]->UpdateSol( _SolSystemPdeIndex, _LinSolver[ig]->_EPS, _LinSolver[ig]->KKoffset );
        _solution[ig]->UpdateRes( _SolSystemPdeIndex, _LinSolver[ig]->_RES, _LinSolver[ig]->KKoffset );
      }

      // ============== Update Fine Residual ==============
      _solution[gridn - 1]->UpdateRes( _SolSystemPdeIndex, _LinSolver[gridn - 1]->_RES, _LinSolver[gridn - 1]->KKoffset );
      linearIsConverged = HasLinearConverged( gridn - 1 );

      if ( linearIsConverged ) break;

    }

    // ============== Update Fine Solution ==============
    _solution[gridn - 1]->UpdateSol( _SolSystemPdeIndex, _LinSolver[gridn - 1]->_EPS, _LinSolver[gridn - 1]->KKoffset );


    std::cout << "\n ************ Linear-Cycle TIME:\t" << std::setw( 11 ) << std::setprecision( 6 ) << std::fixed
              << static_cast<double>( ( clock() - start_mg_time ) ) / CLOCKS_PER_SEC << std::endl;

    return linearIsConverged;
  }


  // ********************************************
  void /*LinearImplicitSystemForSSC::*/SolveForSSC( const unsigned& gridn, const bool &kspClean, const int &npre, const int &npost )
  {

    unsigned grid0 = ( gridn <= _gridr ) ? gridn : _gridr;

    for ( unsigned ig = gridn - 1u; ig > grid0 - 1u; ig-- ) {

      if (_factorTest) {
        _LinSolver[ig]->SetRichardsonScaleFactor(_scale / (gridn - 1u));
      }
      else
        _LinSolver[ig]->SetRichardsonScaleFactor(_scale);

      for ( unsigned k = 0; k < npre; k++ ) {
        _LinSolver[ig]->Solve( _VariablesToBeSolvedIndex, kspClean );
      }
      Restrictor( ig );
    }

    _LinSolver[grid0 - 1u]->Solve( _VariablesToBeSolvedIndex, kspClean );

    for ( unsigned ig = grid0; ig < gridn; ig++ ) {
      Prolongator( ig );
      for ( unsigned k = 0; k < npost; k++ ) {
        _LinSolver[ig]->Solve( _VariablesToBeSolvedIndex, kspClean );
      }
    }
  }

  
  
  
  
  
  protected:

      ///@todo just find where these variables are used and take the corresponding functions!
      
      bool _factorTest;
      double _scale;
      bool _sscLevelSmoother;
                        
      unsigned _sscType[4];
                           
      unsigned _gridr; /*Number of uniformly refined levels*/
      
      
  };


  
  
  
  
  
  
  
  
  
  
  
}
#endif


 

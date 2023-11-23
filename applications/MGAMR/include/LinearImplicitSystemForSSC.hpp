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
                                { };
  

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
      
    protected:

      bool _factorTest;
      double _scale;
      bool _sscLevelSmoother;
                        
      unsigned _sscType[4];
                           
  };


  
  
  
  
  
  
  
  
  
  
  
}
#endif


 

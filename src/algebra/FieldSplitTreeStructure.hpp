#ifndef __femus_enums_FieldSplitTreeStructure_hpp__
#define __femus_enums_FieldSplitTreeStructure_hpp__

#include "SolvertypeEnum.hpp"
#include "PrecondtypeEnum.hpp"

#include <vector>
#include <map>

namespace femus {
  class FieldSpliTreeStructure{
  public:
    //single split constructor
    FieldSpliTreeStructure( SolverType solver, PreconditionerType preconditioner, const std::vector < unsigned > fields){
      _solver = solver;
      _preconditioner = preconditioner;
      _numberOfSplits = 1;
      _branch.resize(0); 
      // _branch.resize(1);
      
      _fields.resize(1);
      _fields[0] = fields;
     
      std::map < unsigned, bool > mymap;
      for(unsigned i = 0; i < _numberOfSplits; i++ ){
	for(unsigned j=0; j < _fields[i].size(); i++){
	  mymap[_fields[i][j]] = true;
	}
      } 
     
      _allFields.resize( mymap.size() );
      
      unsigned j = 0;
      for (std::map<unsigned, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it){
	_allFields[j] = it->first;
	j++;
      }
    };
    
    //multiple split constructor
    FieldSpliTreeStructure( SolverType solver, PreconditionerType preconditioner, std::vector < FieldSpliTreeStructure *> branch){
      
      _solver = solver;
      _preconditioner = preconditioner;
      _numberOfSplits = _branch.size();
      
      _branch.resize(_numberOfSplits);
      _fields.resize(_numberOfSplits);
      for(unsigned i = 0; i < _numberOfSplits; i++ ){
	_branch[i] = branch[i];
	for ( unsigned j = 0; _branch[i]->_allFields.size(); j++){ // we do not specify the size of _allFields.size()
	  _fields[i][j] = _branch[i]->_allFields[j];
	}
      }
      
      
      std::map < unsigned, bool > mymap;
      for(unsigned i = 0; i < _numberOfSplits; i++ ){
	for(unsigned j=0; j < _fields[i].size(); i++){
	  mymap[_fields[i][j]] = true;
	}
      } 
      
      _allFields.resize( mymap.size() );
      
      unsigned j = 0;
      for (std::map<unsigned, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it){
	_allFields[j] = it->first;
	j++;
      }   
    }
    
  SolverType GetSolver(){return _solver;} 
  PreconditionerType GetPreconditioner(){return _preconditioner;} 
  unsigned GetNumberOfSplits(){return _numberOfSplits;} 
  
  FieldSpliTreeStructure* GetBranch(const unsigned &i){return _branch[i];}
  std::vector < unsigned > & GetFields(const unsigned &i) {return _fields[i];} 
  // _fields[i][j] is a vector
  // std::vector < unsigned > GetFields(const unsigned &i) {return _fields[i];} 
  std::vector < unsigned > & GetAllFields() {return _allFields;} 
  // Vector <unsigned> GetAllFields() {return _allFields}
  
  private:
      
    SolverType _solver;
    PreconditionerType _preconditioner; 
    unsigned _numberOfSplits;
    std::vector < FieldSpliTreeStructure* > _branch; // I do not understand this part?
    std::vector < std::vector < unsigned > > _fields;
    std::vector < unsigned > _allFields;
  };
  
  
}


#endif

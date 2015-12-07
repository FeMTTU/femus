#ifndef __femus_enums_FieldSplitTreeStructure_hpp__
#define __femus_enums_FieldSplitTreeStructure_hpp__

#include "SolvertypeEnum.hpp"
#include "PrecondtypeEnum.hpp"

#include <cstdlib>
#include <iostream>

#include <vector>
#include <map>

namespace femus {
  class FieldSpliTreeStructure{
  public:
    //single split constructor
    FieldSpliTreeStructure( const SolverType &solver, const PreconditionerType &preconditioner, const std::vector < unsigned > &fields){
      
      _solver = solver;
      _preconditioner = preconditioner;
      _numberOfSplits = 1;
      _childrenBranches.resize(0); 
      
      
      _fieldsInSplit.resize(1);
      _fieldsInSplit[0] = fields;
     
      //BEGIN ALL FIELD COLLECTION
      std::map < unsigned, bool > mymap;
      for(unsigned i = 0; i < _numberOfSplits; i++ ){
	for(unsigned j=0; j < _fieldsInSplit[i].size(); j++){ //I change "i++" to j++
	  mymap[_fieldsInSplit[i][j]] = true;
	}
      } 
     
      _allFields.resize(mymap.size() );
      
      unsigned j = 0;
      for (std::map<unsigned, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it){
	_allFields[j] = it->first;
	j++;
      }
      //END ALL FIELD COLLECTION
    };
    
    //multiple split constructor
    FieldSpliTreeStructure( const SolverType &solver, const PreconditionerType &preconditioner, const std::vector < FieldSpliTreeStructure> &childrenBranches){
      
      _solver = solver;
      _preconditioner = preconditioner;
      _numberOfSplits = childrenBranches.size();
      _fieldsInSplit.resize(_numberOfSplits);
      
      _childrenBranches.resize(_numberOfSplits);
      
      for(unsigned i = 0; i < _numberOfSplits; i++ ){
	_childrenBranches[i] = &childrenBranches[i];
	_fieldsInSplit[i].resize( childrenBranches[i]._allFields.size() );
	for ( unsigned j = 0; j < childrenBranches[i]._allFields.size(); j++){ 
 	  _fieldsInSplit[i][j] = childrenBranches[i]._allFields[j];
	}
      }
      
      //BEGIN ALL FIELD COLLECTION
      std::map < unsigned, bool > mymap;
      for(unsigned i = 0; i < _numberOfSplits; i++ ){
	for(unsigned j=0; j < _fieldsInSplit[i].size(); j++){ //I change "i++" to j++
	  mymap[_fieldsInSplit[i][j]] = true;
	}
      } 
     
      _allFields.resize(mymap.size() );
      
      unsigned j = 0;
      for (std::map<unsigned, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it){
	_allFields[j] = it->first;
	j++;
      }
      //END ALL FIELD COLLECTION
    }
    
  SolverType GetSolver(){return _solver;} 
  PreconditionerType GetPreconditioner(){return _preconditioner;} 
  unsigned GetNumberOfSplits(){return _numberOfSplits;} 
  
  const FieldSpliTreeStructure * GetChildrenBranch(const unsigned &i){
    if(i < _numberOfSplits){
      return _childrenBranches[i];
    }
    else{
      std::cout << "Wrong input (= " << i << ") in function FieldSpliTreeStructure::GetChildrenBranchstd"<<std::cout;
      std::cout << "Number of Splits = "<<_numberOfSplits<<std::endl;
      abort();
    }
    
  }
  const std::vector < unsigned > & GetFieldsInSplit(const unsigned &i) {return _fieldsInSplit[i];} 
  // _fields[i][j] is a vector
  // std::vector < unsigned > GetFields(const unsigned &i) {return _fields[i];} 
  const std::vector < unsigned > & GetAllFields() {return _allFields;} 
  // Vector <unsigned> GetAllFields() {return _allFields}
  
  private:
      
    SolverType _solver;
    PreconditionerType _preconditioner; 
    unsigned _numberOfSplits;
    std::vector < const FieldSpliTreeStructure * > _childrenBranches; // I do not understand this part?
    std::vector < std::vector < unsigned > > _fieldsInSplit;
    std::vector < unsigned > _allFields;
  };
  
  
}


#endif

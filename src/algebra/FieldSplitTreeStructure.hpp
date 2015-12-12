#ifndef __femus_enums_FieldSplitTreeStructure_hpp__
#define __femus_enums_FieldSplitTreeStructure_hpp__

#include "SolvertypeEnum.hpp"
#include "PrecondtypeEnum.hpp"

#include <cstdlib>
#include <iostream>

#include <vector>
#include <map>
#include <string>

namespace femus {
  class FieldSpliTreeStructure{
  public:
    //single split constructor
    FieldSpliTreeStructure( const SolverType &solver, const PreconditionerType &preconditioner, const std::vector < unsigned > &fields, std::string name){
      
      _name=name;
      
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
    FieldSpliTreeStructure( const SolverType &solver, const PreconditionerType &preconditioner, std::vector < FieldSpliTreeStructure> &childrenBranches, std::string name){
      
      _name=name;
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
    
  const SolverType  & GetSolver(){return _solver;} 
  const PreconditionerType & GetPreconditioner(){return _preconditioner;}
  const unsigned & GetNumberOfSplits(){return _numberOfSplits;}
  
  FieldSpliTreeStructure * GetChildranch(const unsigned &i){
    if(i < _numberOfSplits){
      return _childrenBranches[i];
    }
    else{
      std::cout << "Wrong input (= " << i << ") in function FieldSpliTreeStructure::GetChildrenBranchstd"<<std::cout;
      std::cout << "Number of Splits = "<<_numberOfSplits<<std::endl;
      abort();
    }
    
  }
  
  void PrintNestedFields(const unsigned &counter = 0){

    std::string sub=" ";
    for(int i=0; i < counter; i++){
      sub+="sub-";
    }  
    std::cout<< "Fields in the " << _name << sub <<"system:\n"; 
    for(int j = 0; j < _allFields.size(); j++) std::cout << _allFields[j] << " ";
    std::cout<<std::endl;   
    
    if( GetNumberOfSplits() > 1){
      for(unsigned i = 0; i< GetNumberOfSplits(); i++){
	_childrenBranches[i]->PrintNestedFields(counter+1);
      }
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

    std::vector < FieldSpliTreeStructure * > _childrenBranches; 
    std::vector < FieldSpliTreeStructure* > _branch; 
    std::vector < std::vector < unsigned > > _fields;
    std::vector < std::vector < unsigned > > _fieldsInSplit;
    std::vector < unsigned > _allFields;
    std::string _name;
  };
  
  
}


#endif

#ifndef __femus_enums_FieldSplitTreeStructure_hpp__
#define __femus_enums_FieldSplitTreeStructure_hpp__


#include <cstdlib>
#include <iostream>

#include <vector>
#include <map>
#include <string>

#include "SolvertypeEnum.hpp"
#include "PrecondtypeEnum.hpp"
#include "PetscVector.hpp"
#include "PetscPreconditioner.hpp"


namespace femus {
  class FieldSpliTreeStructure {
    public:
      //single split constructor
      FieldSpliTreeStructure (const SolverType &solver, const PreconditionerType &preconditioner, const std::vector < unsigned > &fields, std::string name) {
        _father = NULL;
        _name = name;

        _solver = solver;
        _preconditioner = preconditioner;
        _numberOfSplits = 1;
        _childBranch.resize (0);


        _fieldsSplit.resize (1);
        _fieldsSplit[0] = fields;

        //BEGIN ALL FIELD COLLECTION
        std::map < unsigned, bool > mymap;
        for (unsigned i = 0; i < _numberOfSplits; i++) {
          for (unsigned j = 0; j < _fieldsSplit[i].size(); j++) { //I change "i++" to j++
            mymap[_fieldsSplit[i][j]] = true;
          }
        }

        _fieldsAll.resize (mymap.size());

        unsigned j = 0;
        for (std::map<unsigned, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it) {
          _fieldsAll[j] = it->first;
          j++;
        }
        //END ALL FIELD COLLECTION
      };

      //multiple split constructor
      FieldSpliTreeStructure (const SolverType &solver, const PreconditionerType &preconditioner, std::vector < FieldSpliTreeStructure*> childBranch, std::string name) {

        _father = NULL;

        _name = name;
        _solver = solver;
        _preconditioner = preconditioner;
        _numberOfSplits = childBranch.size();
        _fieldsSplit.resize (_numberOfSplits);

        _childBranch.resize (_numberOfSplits);

        for (unsigned i = 0; i < _numberOfSplits; i++) {
          childBranch[i]->_father = this;
          _childBranch[i] = childBranch[i];
          _fieldsSplit[i].resize (childBranch[i]->_fieldsAll.size());
          for (unsigned j = 0; j < childBranch[i]->_fieldsAll.size(); j++) {
            _fieldsSplit[i][j] = childBranch[i]->_fieldsAll[j];
          }
        }

        //BEGIN ALL FIELD COLLECTION
        std::map < unsigned, bool > mymap;
        for (unsigned i = 0; i < _numberOfSplits; i++) {
          for (unsigned j = 0; j < _fieldsSplit[i].size(); j++) {
            mymap[_fieldsSplit[i][j]] = true;
          }
        }

        _fieldsAll.resize (mymap.size());

        unsigned j = 0;
        for (std::map<unsigned, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it) {
          _fieldsAll[j] = it->first;
          j++;
        }
        //END ALL FIELD COLLECTION
      }

      const SolverType  & GetSolver() {
        return _solver;
      }
      const PreconditionerType & GetPreconditioner() {
        return _preconditioner;
      }
      const unsigned & GetNumberOfSplits() {
        return _numberOfSplits;
      }

      FieldSpliTreeStructure * GetChildBranch (const unsigned &i) {
        if (i < _numberOfSplits) {
          return _childBranch[i];
        }
        else {
          std::cout << "Wrong input (= " << i << ") in function FieldSpliTreeStructure::GetChildrenBranchstd" << std::cout;
          std::cout << "Number of Splits = " << _numberOfSplits << std::endl;
          abort();
        }

      }

      void PrintNestedFields (const unsigned &counter = 0) {

        std::string sub = " ";
        for (int i = 0; i < counter; i++) {
          sub += "sub-";
        }

        std::cout << "Fields in the " << _name << sub << "system:\n";
        for (int j = 0; j < _fieldsAll.size(); j++) std::cout << _fieldsAll[j] << " ";
        std::cout << std::endl;

        if (_father != NULL) {
          std::cout << "My father is " << _father->GetName() << std::endl << std::endl;
        }

        if (GetNumberOfSplits() > 1) {
          for (unsigned i = 0; i < GetNumberOfSplits(); i++) {
            std::cout << "Split " << i << " ";
            std::cout <<"Size= "<< _isSplitIndex[i].size() << " begin= " << _isSplitIndex[i][0] << " end=" << _isSplitIndex[i][_isSplitIndex[i].size() - 1] << std::endl;
            _childBranch[i]->PrintNestedFields (counter + 1);
          }
        }

      }


      void ReorderFields (const std::vector< std::vector < unsigned > > &KKoffset, const unsigned& iproc, const unsigned& nprocs) {

        _isSplitIndex.resize (GetNumberOfSplits());
	_isSplit.resize (GetNumberOfSplits());
        for (unsigned i = 0; i < GetNumberOfSplits(); i++) {

          //on the actual structure

          unsigned size = 0;
          for (unsigned k = 0; k < _fieldsSplit[i].size(); k++) {
            unsigned index = _fieldsSplit[i][k];
            unsigned offset = KKoffset[index][iproc];
            unsigned offsetp1 = KKoffset[index + 1][iproc];
            size += offsetp1 - offset;
          }

          _isSplitIndex[i].resize (size);

          unsigned counter = 0;
          for (int k = 0; k < _fieldsSplit[i].size(); k++) {
            unsigned index = _fieldsSplit[i][k];
            unsigned offset = KKoffset[index][iproc];
            unsigned offsetp1 = KKoffset[index + 1][iproc];
            for (int j = offset; j < offsetp1; j++) {
              _isSplitIndex[i][counter] = j;
              counter++;
            }
          }
          
      	  ISCreateGeneral(MPI_COMM_WORLD, _isSplitIndex[i].size(), &_isSplitIndex[i][0], PETSC_USE_POINTER, &_isSplit[i]);
      

          // on the child branches

          std::vector < std::vector < unsigned > > tempFields = _childBranch[i]->_fieldsSplit;

          for (unsigned j = 0; j < _childBranch[i]->_fieldsAll.size(); j++) {
            unsigned index = _childBranch[i]->_fieldsAll[j];
            _childBranch[i]->_fieldsAll[j] = j;
            for (unsigned k = 0; k < _childBranch[i]->_numberOfSplits; k++) {
              for (unsigned l = 0; l < _childBranch[i]->_fieldsSplit[k].size(); l++) {
                if (tempFields[k][l] == index) {
                  _childBranch[i]->_fieldsSplit[k][l] = j;

                }
              }
            }
          }

          std::vector< std::vector < unsigned > > fieldsInSplitOffset (_fieldsSplit[i].size() + 1);

          for (unsigned j = 0; j < _fieldsSplit[i].size() + 1; j++) {
            fieldsInSplitOffset[j].resize (nprocs);
          }
          fieldsInSplitOffset[0][0] = 0;
          for (unsigned jproc = 0; jproc < nprocs; jproc++) {
            for (unsigned k = 0; k < _fieldsSplit[i].size(); k++) {
              unsigned index = _fieldsSplit[i][k];
              unsigned ownedDofs =  KKoffset[index + 1][jproc] -  KKoffset[index][jproc];
              fieldsInSplitOffset[k + 1][jproc] = fieldsInSplitOffset[k][jproc] + ownedDofs;
            }
            if (jproc < nprocs - 1) {
              fieldsInSplitOffset[0][jproc + 1] = fieldsInSplitOffset[_fieldsSplit[i].size()][jproc];
            }
          }
          
          for (unsigned jproc = 0; jproc < nprocs; jproc++) {
            for (unsigned k = 0; k <= _fieldsSplit[i].size(); k++) {
	      std::cout << fieldsInSplitOffset[k][jproc] << " ";
	    }
	    std::cout<<std::endl;
	  }
          

          if (_childBranch[i]->GetNumberOfSplits() > 1) {
            _childBranch[i]->ReorderFields (fieldsInSplitOffset, iproc, nprocs);
          }
        }
      }

      
      void SetPC(KSP ksp){
	
	PC pc;
	KSPGetPC(ksp, &pc);

	//BEGIN from here
        if( _preconditioner == FIELDSPLIT_PRECOND ){
	  PetscPreconditioner::set_petsc_preconditioner_type(_preconditioner, pc);
	  PCFieldSplitSetType(pc, PC_COMPOSITE_ADDITIVE);
	  for(int i=0; i<_numberOfSplits; i++ ){
	    PCFieldSplitSetIS( pc, NULL, _isSplit[i]);
	  }
	  KSPSetUp(ksp);
	  KSP* subksp;
	  PetscInt nlocal = static_cast < PetscInt > (_numberOfSplits);
	  PCFieldSplitGetSubKSP(pc, &nlocal, &subksp);
	  for(unsigned i =0; i<_numberOfSplits; i++){
	    _childBranch[i]->SetPC(subksp[i]);
	  }
	}
	else{
	  _rtol = 1.e-3;
	  _abstol = 1.e-20;
	  _dtol = 1.e+50;
	  KSPSetType(ksp, (char*) KSPPREONLY);
	  PC pc;
	  KSPGetPC(ksp, &pc);
	  KSPSetTolerances(ksp, _rtol, _abstol, _dtol, 1);
	  KSPSetFromOptions(ksp);
          PetscReal epsilon = 1.e-16;
	  PetscPreconditioner::set_petsc_preconditioner_type(_preconditioner, pc);
	  PCFactorSetZeroPivot(pc, epsilon);
	  PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);	    
	}
      }

      FieldSpliTreeStructure *GetFather() const {
        if (_father != NULL) {
          return _father;
        }
        else {
          std::cout << "Warning this split has no father" << std::endl;
          abort();
        }
      }

      std::string GetName() const {
        return _name;
      }


      const std::vector < unsigned > & GetFieldsInSplit (const unsigned &i) {
        return _fieldsSplit[i];
      }

      const std::vector < unsigned > & GetAllFields() {
        return _fieldsAll;
      }


    private:

      SolverType _solver;
      PreconditionerType _preconditioner;
      unsigned _numberOfSplits;
      FieldSpliTreeStructure *_father;
      std::vector < FieldSpliTreeStructure * > _childBranch;
      std::vector < std::vector < unsigned > > _fieldsSplit;
      std::vector < unsigned > _fieldsAll;
      std::string _name;
      std::vector < std::vector < PetscInt > > _isSplitIndex; 
      std::vector < IS > _isSplit;
      double _rtol;
      double _abstol;
      double _dtol;
  };


}


#endif

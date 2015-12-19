
#include "FieldSplitTree.hpp"


namespace femus {

  FieldSplitTree::FieldSplitTree( const SolverType& solver, const PreconditionerType& preconditioner, const std::vector < unsigned >& fields, std::string name ) {
    _father = NULL;
    _name = name;

    _solver = solver;
    _preconditioner = preconditioner;
    _numberOfSplits = 1;
    _child.resize( 0 );


    _fieldsSplit.resize( 1 );
    _fieldsSplit[0] = fields;

    //BEGIN ALL FIELD COLLECTION
    std::map < unsigned, bool > mymap;
    for( unsigned i = 0; i < _numberOfSplits; i++ ) {
      for( unsigned j = 0; j < _fieldsSplit[i].size(); j++ ) { //I change "i++" to j++
        mymap[_fieldsSplit[i][j]] = true;
      }
    }

    _fieldsAll.resize( mymap.size() );

    unsigned j = 0;
    for( std::map<unsigned, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it ) {
      _fieldsAll[j] = it->first;
      j++;
    }
    //END ALL FIELD COLLECTION
  };

  //multiple split constructor
  FieldSplitTree::FieldSplitTree( const SolverType& solver, const PreconditionerType& preconditioner, std::vector < FieldSplitTree*> childBranch, std::string name ) {

    _father = NULL;

    _name = name;
    _solver = solver;
    _preconditioner = preconditioner;
    _numberOfSplits = childBranch.size();
    _fieldsSplit.resize( _numberOfSplits );

    _child.resize( _numberOfSplits );

    for( unsigned i = 0; i < _numberOfSplits; i++ ) {
      childBranch[i]->_father = this;
      _child[i] = childBranch[i];
      _fieldsSplit[i].resize( childBranch[i]->_fieldsAll.size() );
      for( unsigned j = 0; j < childBranch[i]->_fieldsAll.size(); j++ ) {
        _fieldsSplit[i][j] = childBranch[i]->_fieldsAll[j];
      }
    }

    //BEGIN ALL FIELD COLLECTION
    std::map < unsigned, bool > mymap;
    for( unsigned i = 0; i < _numberOfSplits; i++ ) {
      for( unsigned j = 0; j < _fieldsSplit[i].size(); j++ ) {
        mymap[_fieldsSplit[i][j]] = true;
      }
    }

    _fieldsAll.resize( mymap.size() );

    unsigned j = 0;
    for( std::map<unsigned, bool>::iterator it = mymap.begin(); it != mymap.end(); ++it ) {
      _fieldsAll[j] = it->first;
      j++;
    }
    //END ALL FIELD COLLECTION
  }


  FieldSplitTree::~FieldSplitTree() {
    if( _numberOfSplits > 1 ) {
      for( unsigned i = 0; i < _isSplit.size(); i++ ) {
        for( unsigned j = 0; j < _isSplit[i].size(); j++ ) {
          ISDestroy( &_isSplit[i][j] );
        }
      }
    }
    for( unsigned i = 0; i < _isSplitIndexPt.size(); i++ ) {
      delete [] _isSplitIndexPt[i];
    }
  }

  void FieldSplitTree::PrintFieldSplitTree( const unsigned& counter ) {

    std::string sub = " ";
    for( int i = 0; i < counter; i++ ) {
      sub += "sub-";
    }

    std::cout << "Fields in the " << _name << sub << "system:\n";
    for( int j = 0; j < _fieldsAll.size(); j++ ) std::cout << _fieldsAll[j] << " ";
    std::cout << std::endl;

    if( _father != NULL ) {
      std::cout << "My father is " << _father->GetName() << std::endl << std::endl;
    }

    if( GetNumberOfSplits() > 1 ) {
      for( unsigned i = 0; i < GetNumberOfSplits(); i++ ) {
        _child[i]->PrintFieldSplitTree( counter + 1 );
      }
    }

  }


  void FieldSplitTree::BuildIndexSet( const std::vector< std::vector < unsigned > >& KKoffset, const unsigned& iproc, const unsigned& nprocs, const unsigned& level ) {

    if( _isSplit.size() < level ) _isSplit.resize( level );
    _isSplit[level - 1].resize( GetNumberOfSplits() );

    for( unsigned i = 0; i < GetNumberOfSplits(); i++ ) {

      //on the actual structure
      unsigned size = 0;
      for( unsigned k = 0; k < _fieldsSplit[i].size(); k++ ) {
        unsigned index = _fieldsSplit[i][k];
        unsigned offset = KKoffset[index][iproc];
        unsigned offsetp1 = KKoffset[index + 1][iproc];
        size += offsetp1 - offset;
      }


      PetscInt* isSplitIndex = new PetscInt [size];
      unsigned ptSize = _isSplitIndexPt.size();
      _isSplitIndexPt.resize( ptSize + 1 );
      _isSplitIndexPt[ptSize] = isSplitIndex;

      unsigned counter = 0;
      for( int k = 0; k < _fieldsSplit[i].size(); k++ ) {
        unsigned index = _fieldsSplit[i][k];
        unsigned offset = KKoffset[index][iproc];
        unsigned offsetp1 = KKoffset[index + 1][iproc];
        for( int j = offset; j < offsetp1; j++ ) {
          isSplitIndex[counter] = j;
          counter++;
        }
      }

      ISCreateGeneral( MPI_COMM_WORLD, size, isSplitIndex, PETSC_USE_POINTER, &_isSplit[level - 1][i] );

      // on the child branches

      std::vector < std::vector < unsigned > > tempFields = _child[i]->_fieldsSplit;

      for( unsigned j = 0; j < _child[i]->_fieldsAll.size(); j++ ) {
        unsigned index = _child[i]->_fieldsAll[j];
        _child[i]->_fieldsAll[j] = j;
        for( unsigned k = 0; k < _child[i]->_numberOfSplits; k++ ) {
          for( unsigned l = 0; l < _child[i]->_fieldsSplit[k].size(); l++ ) {
            if( tempFields[k][l] == index ) {
              _child[i]->_fieldsSplit[k][l] = j;

            }
          }
        }
      }

      if( _child[i]->GetNumberOfSplits() > 1 ) {

        std::vector< std::vector < unsigned > > fieldsInSplitOffset( _fieldsSplit[i].size() + 1 );

        for( unsigned j = 0; j < _fieldsSplit[i].size() + 1; j++ ) {
          fieldsInSplitOffset[j].resize( nprocs );
        }
        fieldsInSplitOffset[0][0] = 0;
        for( unsigned jproc = 0; jproc < nprocs; jproc++ ) {
          for( unsigned k = 0; k < _fieldsSplit[i].size(); k++ ) {
            unsigned index = _fieldsSplit[i][k];
            unsigned ownedDofs =  KKoffset[index + 1][jproc] -  KKoffset[index][jproc];
            fieldsInSplitOffset[k + 1][jproc] = fieldsInSplitOffset[k][jproc] + ownedDofs;
          }
          if( jproc < nprocs - 1 ) {
            fieldsInSplitOffset[0][jproc + 1] = fieldsInSplitOffset[_fieldsSplit[i].size()][jproc];
          }
        }
        _child[i]->BuildIndexSet( fieldsInSplitOffset, iproc, nprocs, level );
      }
    }
  }


  void FieldSplitTree::SetPC( KSP& ksp, const unsigned& level ) {

    PC pc;
    KSPGetPC( ksp, &pc );

    //BEGIN from here
    if( _preconditioner == FIELDSPLIT_PRECOND ) {
      PetscPreconditioner::set_petsc_preconditioner_type( _preconditioner, pc );
      PCFieldSplitSetType( pc, PC_COMPOSITE_ADDITIVE );
      for( int i = 0; i < _numberOfSplits; i++ ) {
        PCFieldSplitSetIS( pc, NULL, _isSplit[level - 1][i] );
      }
      KSPSetUp( ksp );
      KSP* subksp;
      PetscInt nlocal = static_cast < PetscInt >( _numberOfSplits );
      PCFieldSplitGetSubKSP( pc, &nlocal, &subksp );
      for( unsigned i = 0; i < _numberOfSplits; i++ ) {
        _child[i]->SetPC( subksp[i], level );
      }
    }
    else {
      _rtol = 1.e-3;
      _abstol = 1.e-20;
      _dtol = 1.e+50;
      _maxits = 1;
      KSPSetType( ksp, ( char* ) KSPPREONLY );
      PC pc;
      KSPGetPC( ksp, &pc );
      KSPSetTolerances( ksp, _rtol, _abstol, _dtol, _maxits );
      KSPSetFromOptions( ksp );
      PetscReal epsilon = 1.e-16;
      PetscPreconditioner::set_petsc_preconditioner_type( _preconditioner, pc );
      PCFactorSetZeroPivot( pc, epsilon );
      PCFactorSetShiftType( pc, MAT_SHIFT_NONZERO );
    }
  }

  FieldSplitTree* FieldSplitTree::GetFather() const {
    if( _father != NULL ) {
      return _father;
    }
    else {
      std::cout << "Warning this split has no father" << std::endl;
      abort();
    }
  }

  FieldSplitTree* FieldSplitTree::GetChild( const unsigned& i ) {
    if( i < _numberOfSplits ) {
      return _child[i];
    }
    else {
      std::cout << "Wrong input (= " << i << ") in function FieldSplitTree::GetChildrenBranchstd" << std::cout;
      std::cout << "Number of Splits = " << _numberOfSplits << std::endl;
      abort();
    }

  }



}

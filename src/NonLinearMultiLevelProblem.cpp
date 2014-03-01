#include "NonLinearMultiLevelProblem.hpp"
#include "ElemType.hpp"
#include "Elem.hpp"
#include "NumericVector.hpp"
#include "SparseMatrix.hpp"
#include "LinearSolver.hpp"
#include "FEMTTUConfig.h"
#include "Parameter.hpp"
#include "hdf5.h"

//C++ include
#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>   
#include <string.h>
using std::cout;
using std::endl;
using std::min;
using std::string;

bool (* mesh::_SetRefinementFlag)(const double &x, const double &y, const double &z, 
				  const int &ElemGroupNumber,const int &level) = NULL;

//---------------------------------------------------------------------------------------------------
NonLinearMultiLevelProblem::~NonLinearMultiLevelProblem() {

  for (unsigned i=0; i<_gridn; i++) {
    delete _solution[i];
    delete _msh[i];
  }

  for (unsigned i=0; i<3; i++)
    for (unsigned j=0; j<3; j++)
      if (1!=i || 2!=j) delete type_elem[i][j];

  for (unsigned i=3; i<5; i++)
    for (unsigned j=0; j<3; j++)
      if (4!=i || 2!=j) delete type_elem[i][j];

  delete type_elem[5][0];
  delete type_elem[5][1];

  delete type_elem[0][3];
  delete type_elem[0][4];

  delete type_elem[3][3];
  delete type_elem[3][4];


  for (unsigned i=0; i<SolName.size(); i++) delete [] SolName[i];
  for (unsigned i=0; i<SolName.size(); i++) delete [] BdcType[i];


};

//---------------------------------------------------------------------------------------------------
NonLinearMultiLevelProblem::NonLinearMultiLevelProblem(const unsigned short &igridn,const unsigned short &igridr,
						       const char mesh_file[], const char GaussOrder[],
						       const double Lref, bool (* SetRefinementFlag)(const double &x, const double &y, const double &z, 
												     const int &ElemGroupNumber,const int &level)):_gridn(igridn), _gridr(igridr) {
		       
  MPI_Comm_rank(MPI_COMM_WORLD, &_iproc);
  MPI_Comm_size(MPI_COMM_WORLD, &_nprocs);

  type_elem[0][0]=new const elem_type("hex","linear",GaussOrder);
  type_elem[0][1]=new const elem_type("hex","quadratic",GaussOrder);
  type_elem[0][2]=new const elem_type("hex","biquadratic",GaussOrder);
  type_elem[0][3]=new const elem_type("hex","constant",GaussOrder);
  type_elem[0][4]=new const elem_type("hex","disc_linear",GaussOrder);


  type_elem[1][0]=new const elem_type("tet","linear",GaussOrder);
  type_elem[1][1]=new const elem_type("tet","biquadratic",GaussOrder);
  type_elem[1][2]=type_elem[1][1];

  type_elem[2][0]=new const elem_type("wedge","linear",GaussOrder);
  type_elem[2][1]=new const elem_type("wedge","quadratic",GaussOrder);
  type_elem[2][2]=new const elem_type("wedge","biquadratic",GaussOrder);

  type_elem[3][0]=new const elem_type("quad","linear",GaussOrder);
  type_elem[3][1]=new const elem_type("quad","quadratic",GaussOrder);
  type_elem[3][2]=new const elem_type("quad","biquadratic",GaussOrder);
  type_elem[3][3]=new const elem_type("quad","constant",GaussOrder);
  type_elem[3][4]=new const elem_type("quad","disc_linear",GaussOrder);


  type_elem[4][0]=new const elem_type("tri","linear",GaussOrder);
  type_elem[4][1]=new const elem_type("tri","biquadratic",GaussOrder);
  type_elem[4][2]=type_elem[4][1];

  type_elem[5][0]=new const elem_type("line","linear",GaussOrder);
  type_elem[5][1]=new const elem_type("line","biquadratic",GaussOrder);
  type_elem[5][2]=type_elem[5][1];

  cout << "MESH DATA: " << endl;

  //Steady State simulation
  _test_time = 0;

  //Temporary coordinates vector
  vector <vector <double> > vt;  
  vt.resize(3);

  _msh.resize(_gridn);
  _solution.resize(_gridn);
  
  _msh[0]=new mesh(mesh_file, vt,Lref);
  _solution[0]=new Solution(_msh[0]);
  
  unsigned gridn_temp=_gridn;
  _gridn=1;
  AddSolution("X","biquadratic",1,0);
  AddSolution("Y","biquadratic",1,0);
  AddSolution("Z","biquadratic",1,0);
  sprintf(BdcType[0],"Steady");
  sprintf(BdcType[1],"Steady");
  sprintf(BdcType[2],"Steady");
    
  _solution[0]->ResizeSolutionVector("X");
  _solution[0]->ResizeSolutionVector("Y");
  _solution[0]->ResizeSolutionVector("Z");
    
  _gridn=gridn_temp;
  _solution[0]->SetCoarseCoordinates(vt);
    
  unsigned indX=GetIndex("X");
  unsigned indY=GetIndex("Y");
  unsigned indZ=GetIndex("Z");
    
  for (unsigned i=1; i<_gridr; i++) {
    _solution[i-1u]->SetElementRefiniement(1);
    _msh[i] = new mesh(i,_msh[i-1]->el); 
        
    _solution[i]=new Solution(_msh[i]);
    _solution[i]->AddSolution("X","biquadratic",1,0);
    _solution[i]->AddSolution("Y","biquadratic",1,0);
    _solution[i]->AddSolution("Z","biquadratic",1,0);
    _solution[i]->ResizeSolutionVector("X");
    _solution[i]->ResizeSolutionVector("Y");
    _solution[i]->ResizeSolutionVector("Z");
             
    BuildProlongatorMatrix(i, indX);
    unsigned TypeIndex=SolType[indX];
    
    _solution[i]->_Sol[indX]->matrix_mult(*_solution[i-1]->_Sol[indX],*_solution[i]->_ProjMat[TypeIndex]);
    _solution[i]->_Sol[indY]->matrix_mult(*_solution[i-1]->_Sol[indY],*_solution[i]->_ProjMat[TypeIndex]);
    _solution[i]->_Sol[indZ]->matrix_mult(*_solution[i-1]->_Sol[indZ],*_solution[i]->_ProjMat[TypeIndex]);
    _solution[i]->_Sol[indX]->close();
    _solution[i]->_Sol[indY]->close();
    _solution[i]->_Sol[indZ]->close(); 
  }
  
  for (unsigned i=_gridr; i<_gridn; i++) {
    if(SetRefinementFlag==NULL) {
      cout << "Set Refinement Region flag is not defined! " << endl;
      exit(1);
    }
    else {
      mesh::_SetRefinementFlag = SetRefinementFlag;
      _solution[i-1u]->SetElementRefiniement(2);
    }
    _msh[i] = new mesh(i,_msh[i-1]->el); 
    
    _solution[i]=new Solution(_msh[i]);
    _solution[i]->AddSolution("X","biquadratic",1,0);
    _solution[i]->AddSolution("Y","biquadratic",1,0);
    _solution[i]->AddSolution("Z","biquadratic",1,0);
    _solution[i]->ResizeSolutionVector("X");
    _solution[i]->ResizeSolutionVector("Y");
    _solution[i]->ResizeSolutionVector("Z");
     
    BuildProlongatorMatrix(i, indX);
    unsigned TypeIndex=SolType[indX];
    
    _solution[i]->_Sol[indX]->matrix_mult(*_solution[i-1]->_Sol[indX],*_solution[i]->_ProjMat[TypeIndex]);
    _solution[i]->_Sol[indY]->matrix_mult(*_solution[i-1]->_Sol[indY],*_solution[i]->_ProjMat[TypeIndex]);
    _solution[i]->_Sol[indZ]->matrix_mult(*_solution[i-1]->_Sol[indZ],*_solution[i]->_ProjMat[TypeIndex]);
    _solution[i]->_Sol[indX]->close();
    _solution[i]->_Sol[indY]->close();
    _solution[i]->_Sol[indZ]->close();
        
  }
  _solution[_gridn-1u]->SetElementRefiniement(0);	
  elr_old.resize(_msh[_gridr-1u]->GetElementNumber());
    
  unsigned refindex = _msh[0]->GetRefIndex();
  // cout << "******" << refindex << endl;
  elem_type::_refindex=refindex;

  cout << endl;
  _ntime_steps=0;
  _dt=1.;
  _time_step0 = 1;
  _time=0.;
  _time_step=0;
  _moving_mesh=0;
  _Schur=0;
  _NSchurVar=1;
  _is_nonlinear = false;
  _non_linear_toll = 1.e-03;
  _non_linear_algorithm = 0;
  _init_func_set=false;
  _bdc_func_set=false;
  _VankaIsSet = true;

  BuildProlongatorMatrices();
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::MarkStructureNode() {
  for (unsigned i=0; i<_gridn; i++) _msh[i]->AllocateAndMarkStructureNode();
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::Add_Fluid(Fluid *fluid) {
  _fluid = fluid;
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::Add_Solid(Solid *solid) {
  _solid = solid;
}

//---------------------------------------------------------------------------------------------------
unsigned NonLinearMultiLevelProblem::GetNumberOfGrid() {
  return _gridn;
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::SetMatrixProperties(const char pdename[], const char property[]) {
  unsigned ipde=GetPdeIndex(pdename);
  if (!strcmp(property,"Symmetric")) {
    const bool mprop = true;
    for (unsigned i=0; i<_gridn; i++) _LinSolver[ipde][i]->SetMatrixProperties(mprop);
  } else {
    cout<<"Error! This option is not admitted \"All\""<<endl;
    exit(1);
  }
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::AddStabilization(const char pdename[], const bool stab, const double compressibility) {
  unsigned ipde=GetPdeIndex(pdename);
  for (unsigned i=0; i<_gridn; i++) _LinSolver[ipde][i]->AddStabilization(stab, compressibility);
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::SetVankaSchurOptions(bool Schur, short unsigned NSchurVar) {
  
  if(Schur==1 && NSchurVar ==0){
    cout<<"Error incompatible options in NonLinearMultiLevelProblem::SetVankaSchurOptions "<<endl;
    exit(0);
  }
  _Schur=Schur;
  _NSchurVar=NSchurVar;
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::SetTolerances(const char pdename[],const double rtol, const double atol,
					       const double divtol, const unsigned maxits) {
  unsigned ipde=GetPdeIndex(pdename);					       
  for (unsigned i=1; i<_gridn; i++) {
    _LinSolver[ipde][i]->set_tolerances(rtol,atol,divtol,maxits,0);
  }
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::SetSchurTolerances(const char pdename[], const double rtol, const double atol,
						    const double divtol, const unsigned maxits) {
  unsigned ipde=GetPdeIndex(pdename);
  for (unsigned i=1; i<_gridn; i++) {
    _LinSolver[ipde][i]->set_tolerances(rtol,atol,divtol,maxits,1);
  }
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::SetDirichletBCsHandling(const char pdename[],const char DirichletMode[]) {
  unsigned ipde=GetPdeIndex(pdename);     
  unsigned int DirichletBCsHandlingMode;
  
  if (!strcmp(DirichletMode,"Penalty")) {
    DirichletBCsHandlingMode = 0;   
  }
  else if (!strcmp(DirichletMode,"Elimination")) {
    DirichletBCsHandlingMode = 1;   
  } 
  else {
    cout << "Error! The Dirichlet BCs Handling method " << DirichletMode  << " is not implemented"<<endl;
    exit(1);
  }
  
  for (unsigned i=0; i<_gridn; i++) {
    _LinSolver[ipde][i]->set_dirichletBCsHandling(DirichletBCsHandlingMode);
  }
}

//---------------------------------------------------------------------------------------------------
unsigned NonLinearMultiLevelProblem::GetTmOrder(const unsigned i) {
  return SolTmorder[i];
};

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::SetSmoother(const char smoothername[]) {
  if (!strcmp(smoothername,"Vanka")) {
    _VankaIsSet = true;
  } else if (!strcmp(smoothername,"Gmres")) {
    _VankaIsSet = false;
  } else {
    cout<<"Error! The smoother " << smoothername  << " is not implemented"<<endl;
    exit(1);
  }
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::SetDimVankaBlock(const char pdename[],unsigned const dim_vanka_block) {

  unsigned ipde=GetPdeIndex(pdename);
  const unsigned dim = _msh[0]->GetDimension();
  const unsigned base = pow(2,dim);
  unsigned num_vanka_block = pow(base,dim_vanka_block);

  for (unsigned i=1; i<_gridn; i++) {
    unsigned num_vanka_block2 = min(num_vanka_block,_msh[i]->GetElementNumber());
    _LinSolver[ipde][i]->set_num_elem_vanka_block(num_vanka_block2);
  }
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::SetDimVankaBlock(const char pdename[], const char dim_vanka_block[]) {
  unsigned ipde=GetPdeIndex(pdename);
  if (!strcmp(dim_vanka_block,"All")) {
    for (unsigned i=1; i<_gridn; i++) {
      unsigned num_vanka_block2=_msh[i]->IS_Mts2Gmt_elem_offset[_msh[i]->_iproc+1]-_msh[i]->IS_Mts2Gmt_elem_offset[_msh[i]->_iproc]; 
      _LinSolver[ipde][i]->set_num_elem_vanka_block(num_vanka_block2);
    }
  } 
  else {
    cout<<"Error! option \""<< dim_vanka_block<<"\" is not admitted in function"
        <<" NonLinearMultiLevelProblem::SetDimVankaBlock(const char [], const char [])"<<endl;
    exit(1);
  }
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::SetSolverFineGrids(const char pdename[], const char solvertype[]) {
  unsigned ipde=GetPdeIndex(pdename);
  if (!strcmp(solvertype,"GMRES")) {
    for (unsigned i=1; i<_gridn; i++) {
      _LinSolver[ipde][i]->set_solver_type(GMRES);
    }
  } else {
    cout<<"Error! The solver " <<  solvertype << " is not implemented"<<endl;
    exit(1);
  }
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::SetPreconditionerFineGrids(const char pdename[],const char preconditioner_type[]) {
  unsigned ipde=GetPdeIndex(pdename);
  if (!strcmp(preconditioner_type,"LU")) {
    for (unsigned i=1; i<_gridn; i++) {
      _LinSolver[ipde][i]->set_preconditioner_type(MLU_PRECOND);
    }
  } else if (!strcmp(preconditioner_type,"ILU")) {
    for (unsigned i=1; i<_gridn; i++) {
      _LinSolver[ipde][i]->set_preconditioner_type(ILU_PRECOND);
    }
  } else if (!strcmp(preconditioner_type,"JACOBI")) {
    for (unsigned i=1; i<_gridn; i++) {
      _LinSolver[ipde][i]->set_preconditioner_type(JACOBI_PRECOND);
    }
  } else if (!strcmp(preconditioner_type,"NO_PRECONDITIONING")) {
    for (unsigned i=1; i<_gridn; i++) {
      _LinSolver[ipde][i]->set_preconditioner_type(IDENTITY_PRECOND);
    }
  } else {
    cout<<"Error! The " <<  preconditioner_type << " preconditioner is not implemented"<<endl;
    exit(1);
  }
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::SetNonLinearAlgorithm(const bool isnonlinear,const char nnlinalg[], 
						       const double nl_toll) {

  _is_nonlinear = isnonlinear;
  _non_linear_toll = nl_toll;

  if (isnonlinear) {

    if (!strcmp(nnlinalg,"Quasi-Newton")) {
      _non_linear_algorithm = 1;
    } else if (!strcmp(nnlinalg,"Newton")) {
      _non_linear_algorithm = 2;
    } else if (!strcmp(nnlinalg,"Linear")) {
      _non_linear_algorithm = 0;
    } else {
      cout<<"Error! The "<< nnlinalg  << " algorithm is not implemented"<<endl;
      exit(1);
    }
  } else {
    _non_linear_algorithm = 0;
  }

}

//---------------------------------------------------------------------------------------------------
unsigned NonLinearMultiLevelProblem::GetNonLinearAlgorithm() {
  return _non_linear_algorithm;
}

//---------------------------------------------------------------------------------------------------
bool NonLinearMultiLevelProblem::GetNonLinearCase() {
  return _is_nonlinear;
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::SetMovingMesh(std::vector<std::string>& movvars_in) {
  _moving_mesh = 1;
  _moving_vars = movvars_in;
}

// 
// TODO Create a function that compute the l2norm of a quantity
//
//---------------------------------------------------------------------------------------------------
double NonLinearMultiLevelProblem::ComputeL2norm() {

  //   double _Td=2;
  //   int node2[27];
  //   double vx[3][27];
  //   double phi2[27],gradphi2[27][3],Weight2;
  // 
  //   unsigned order_ind2 = _LinSolver[0][0]->SolType[GetIndex("T")];
  //   unsigned end_ind2   = _LinSolver[0][0]->END_IND[order_ind2];

  double Error=0;

  return Error;
}

// 
// TODO Create a function that compute the Force on the boundary
//
//---------------------------------------------------------------------------------------------------
int NonLinearMultiLevelProblem::ComputeBdStress(int bd, double Cforce[3]) {

  int node2[27];
  double vx[3][27];
  double phi2[27],gradphi2[27][3],Weight2;
  double normal[3];
  //const double rho = _fluid->get_density();
  //const double mu = _fluid->get_viscosity();
  //const double uref = _fluid->_parameter.Get_reference_velocity();

  double BdIntx=0;
  double BdInty=0;
  double BdIntz=0;
  PetscScalar *MYSOL[1];
  unsigned kP=GetIndex("P");
  //unsigned kU=GetIndex("U");
  PetscErrorCode ierr;
  const short unsigned NV1[6][2]= {{9,9},{6,6},{9,6},{3,3},{3,3},{1,1}};
  const unsigned FELT[6][2]= {{3,3},{4,4},{3,4},{5,5},{5,5}};

  unsigned indX=GetIndex("X");
  unsigned indY=GetIndex("Y");
  unsigned indZ=GetIndex("Z");


  for (unsigned ig=0; ig<_gridn; ig++) {

   
    //     PetscVector* petsc_vec_solP = static_cast<PetscVector*>(_solution[ig]->_Sol[kP]);
    //     ierr = VecGetArray(petsc_vec_solP->vec(),&MYSOL[0]);
    //     CHKERRQ(ierr);

    //tested up to now for quad9 in 2D
    unsigned order_ind = SolType[GetIndex("U")];

    for (unsigned iel=0; iel<_msh[ig]->GetElementNumber(); iel++) {
      if (_msh[ig]->el->GetRefinedElementIndex(iel)==0 || ig==_gridn-1u) {

        short unsigned kelt=_msh[ig]->el->GetElementType(iel);
        unsigned nfaces = _msh[ig]->el->GetElementFaceNumber(iel);

        for (unsigned jface=0; jface<nfaces; jface++) {
          if (_msh[ig]->el->GetFaceElementIndex(iel,jface)<0) {

            int dom = -( (_msh[ig]->el->GetFaceElementIndex(iel,jface))+1);
            if (dom==bd) {

              unsigned nve=NV1[kelt][jface<(_msh[ig]->el->GetElementFaceNumber(iel,0))];
	      //               //linear pressure
              unsigned felt = FELT[kelt][jface< (_msh[ig]->el->GetElementFaceNumber(iel,0))];

	      for(unsigned i=0;i<nve;i++) {
                unsigned inode=_msh[ig]->el->GetFaceVertexIndex(iel,jface,i)-1u;
                unsigned inode_Metis=_msh[ig]->GetMetisDof(inode,2);
		
		vx[0][i]=(*_solution[ig]->_Sol[indX])(inode_Metis);  
                vx[1][i]=(*_solution[ig]->_Sol[indY])(inode_Metis);
                vx[2][i]=(*_solution[ig]->_Sol[indZ])(inode_Metis);
		
              }
	      
              for (unsigned igs=0; igs < type_elem[felt][order_ind]->GetGaussPointNumber(); igs++) {
                (type_elem[felt][order_ind]->*type_elem[felt][order_ind]->Jacobian_sur_ptr)(vx,igs,Weight2,phi2,gradphi2,normal);
                for (unsigned i=0; i<nve; i++) {
		  BdIntx += phi2[i]*MYSOL[0][iel]*normal[0]*Weight2;
		  BdInty += phi2[i]*MYSOL[0][iel]*normal[1]*Weight2;
		  BdIntz += phi2[i]*MYSOL[0][iel]*normal[2]*Weight2;
                }
              } //end loop on gauss point
            }
          }
        }
      }
    }

    //ierr = VecRestoreArray(petsc_vec_solP->vec(),&MYSOL[0]);
    CHKERRQ(ierr);

  }

  Cforce[0] = BdIntx;
  Cforce[1] = BdInty;
  Cforce[2] = BdIntz;

  return ierr;
}

// 
// This function computes the integral on the boundary of a Pde like Navier-stokes or Heat equation
// It has been tested on flat boundary and HEX27 and Quad9
//
//---------------------------------------------------------------------------------------------------
int NonLinearMultiLevelProblem::ComputeBdIntegral(const char pdename[],const char var_name[], const unsigned & kel, const unsigned & jface, unsigned level, unsigned dir) {
  
  
  unsigned ipde=GetPdeIndex(pdename);
  //   PetscVector* RESp=static_cast<PetscVector*> (_LinSolver[ipde][level]->_RES);  //TODO
  //   Vec RES=RESp->vec(); //TODO
  
  int ierr;
  double tau;
  double vx[3][27];
  double phi[27],gradphi[27][3],Weight;
  double normal[3];
  int node[27];
     
  unsigned indexvar = GetSolPdeIndex(pdename,var_name);
  short unsigned kelt = _msh[level]->el->GetElementType(kel);
  unsigned order_ind = SolType[GetIndex(var_name)];
  unsigned indX=GetIndex("X");
  unsigned indY=GetIndex("Y");
  unsigned indZ=GetIndex("Z");
		
  bool test = _SetBoundaryConditionFunction(0.,0.,0.,var_name,tau,-(_msh[level]->el->GetFaceElementIndex(kel,jface)+1),_time);
  if(!test) {
		
    const short unsigned NV1[6][2]={{9,9},{6,6},{9,6},{3,3},{3,3},{1,1}};
    unsigned nve=NV1[kelt][jface<_msh[level]->el->GetElementFaceNumber(kel,0)];
    const unsigned FELT[6][2]={{3,3},{4,4},{3,4},{5,5},{5,5}};
    unsigned felt = FELT[kelt][jface<_msh[level]->el->GetElementFaceNumber(kel,0)];

    //		cout << "felt : "  << felt << endl;
    // 		cout << "node: " << nve << endl;
		
    for(unsigned i=0;i<nve;i++) {
      unsigned inode=_msh[level]->el->GetFaceVertexIndex(kel,jface,i)-1u;
      node[i] = inode + _LinSolver[ipde][level]->KKIndex[indexvar];
      unsigned inode_Metis=_msh[level]->GetMetisDof(inode,2);
      
      vx[0][i]=(*_solution[level]->_Sol[indX])(inode_Metis);  
      vx[1][i]=(*_solution[level]->_Sol[indY])(inode_Metis);
      vx[2][i]=(*_solution[level]->_Sol[indZ])(inode_Metis);
      
      
    }

    for(unsigned igs=0;igs < type_elem[felt][order_ind]->GetGaussPointNumber(); igs++) {
      (type_elem[felt][order_ind]->*type_elem[felt][order_ind]->Jacobian_sur_ptr)(vx,igs,Weight,phi,gradphi,normal);
      //   		  cout << Weight << endl;
      //  		  cout << "normal_x: " << normal[0] << endl;
      //  		  cout << "normal_y: " << normal[1] << endl;
      // 		  cout << "normal_z: " << normal[2] << endl;
      //  		  cout << type_elem[felt][order_ind]->GetGaussPointNumber() << endl;
      //  		  cout << "  " << endl;
      for(unsigned i=0;i<nve;i++) {
	//     	          cout << phi[i] << "  " << nve << endl;
	_SetBoundaryConditionFunction(vx[0][i],vx[1][i],vx[2][i],var_name,tau,-(_msh[level]->el->GetFaceElementIndex(kel,jface)+1),_time);
 	    
	//Manca la moltiplicazione x la normale che e'ancora da fare
	PetscScalar value = -phi[i]*tau*normal[dir]*Weight*_dt;
		    
	// Non voglio chiamare Vecsetvalue ma aggiungere il valore direttamente a F
	// per fare questo mi serve la relazione tra i(node locale di surface) e il nodo locale di volume
	_LinSolver[ipde][level]->_RES->add(node[i],value);
      }
    }
  }
  return ierr;
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::CreatePdeStructure() {
  _LinSolver.resize(_PdeIndex.size());
  for(unsigned ipde=0;ipde<_PdeIndex.size();ipde++){
    _LinSolver[ipde].resize(_gridn);
    for(unsigned i=0;i<_gridn;i++){
      _LinSolver[ipde][i]=LinearSolver::build(i,_msh[i]).release();
      //_LinSolver[ipde][i]->SetBdcPointer(&_solution[i]->_Bdc);
    }
    
    for (unsigned i=0; i<_gridn; i++) {
      _LinSolver[ipde][i]->InitPde(_SolPdeIndex[ipde],SolType,SolName,&_solution[i]->_Bdc,_gridr,_gridn);
    }
    for (unsigned ig=1; ig<_gridn; ig++) {
      BuildProlongatorMatrix(ig,_PdeName[ipde]);
    }
  }
  
  return;
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::DeletePdeStructure() {
  for(unsigned ipde=0;ipde<_PdeIndex.size();ipde++)
    for (unsigned ig=0; ig<_gridn; ig++) {
      _LinSolver[ipde][ig]->DeletePde();
      delete _LinSolver[ipde][ig];
    }
  
  for (unsigned i=0; i<_PdeName.size(); i++){ 
    delete [] _PdeName[i];
  }
  
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::AttachAssembleFunction (  int (*function)(NonLinearMultiLevelProblem &mg, unsigned level, const unsigned &gridn, const unsigned &ipde, const bool &assembe_matrix) ) {
  _assemble_function = function;
  return;
}

//---------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::AttachSetBoundaryConditionFunction ( bool (* SetBoundaryConditionFunction) (const double &x, const double &y, const double &z,const char name[], 
													     double &value, const int FaceName, const double time) ) {
 
  _bdc_func_set = true;
  _SetBoundaryConditionFunction = SetBoundaryConditionFunction;
  return;
}

// *******************************************************
void NonLinearMultiLevelProblem::AttachInitVariableFunction ( double (* InitVariableFunction)(const double &x, const double &y, 
											      const double &z,const char name[]) ) {
  _init_func_set = true;
  _InitVariableFunction = InitVariableFunction;
  return;
}

//--------------------------------------------------------------------------------------------------
void NonLinearMultiLevelProblem::Solve(const char pdename[], unsigned const &Vcycle_number, unsigned const &npre, 
				       unsigned const &npost, const char multigrid_type[], const bool &test_linear) {
  clock_t start_mg_time = clock();
  
  unsigned nonlinear_Vcycle_numeber    = (test_linear == false) ? Vcycle_number:1;
  unsigned linear_Vcycle_numeber       = (test_linear == true ) ? Vcycle_number:1;
  bool full_cycle = (!strcmp(multigrid_type,"V-Cycle")) ? 0 : 1;
  
  unsigned ipde = GetPdeIndex(pdename);
    
  std::pair<int, double> solver_info;
     
  for ( unsigned igridn=full_cycle + (!full_cycle)*_gridr; igridn <= _gridn; igridn++) {
    cout << endl << "    ************* Level Max: " << igridn << " *************\n" << endl;
    for ( unsigned nonlinear_cycle = 0; nonlinear_cycle < nonlinear_Vcycle_numeber; nonlinear_cycle++ ) { //non linear cycle
      clock_t start_cycle_time = clock();
      cout << endl << "    ************** Cycle: " << nonlinear_cycle + 1 << " ****************\n" << endl;
     
      // ============== Fine level Assembly ==============
      clock_t start_time = clock();
      _LinSolver[ipde][igridn-1u]->SetResZero();
      _LinSolver[ipde][igridn-1u]->SetEpsZero();
      bool assemble_matrix = true; //Be carefull!!!! this is needed in the _assemble_function
      _assemble_function(*this,igridn-1u, igridn-1u,ipde, assemble_matrix);
      cout << "Grid: " << igridn-1 << "\t        ASSEMBLY TIME:\t"<<static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << endl;
 
      for(int linear_cycle = 0; linear_cycle < linear_Vcycle_numeber; linear_cycle++){ //linear cycle

	for (unsigned ig = igridn-1u; ig > 0; ig--) {
	  
	  // ============== Presmoothing ============== 
	  for (unsigned k = 0; k < npre; k++) {
	    solver_info = (_VankaIsSet) ? _LinSolver[ipde][ig]->solve(VankaIndex, _NSchurVar, _Schur) : _LinSolver[ipde][ig]->solve();
	  }
	  // ============== Non-Standard Multigrid Restriction ==============
	  start_time = clock();
	  Restrictor(ipde, ig, igridn, nonlinear_cycle, linear_cycle, full_cycle);
	  cout << "Grid: " << ig << "-->" << ig-1 << "  RESTRICTION TIME:\t"<<static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << endl;
	}
      
	// ============== Coarse Direct Solver ==============
	solver_info = ( _VankaIsSet ) ? _LinSolver[ipde][0]->solve(VankaIndex, _NSchurVar, _Schur) : _LinSolver[ipde][0]->solve();
	
            
	for (unsigned ig = 1; ig < igridn; ig++) {
	  
	  // ============== Standard Prolongation ==============
	  start_time=clock();
	  Prolongator(ipde, ig);
	  cout << "Grid: " << ig-1 << "-->" << ig << " PROLUNGATION TIME:\t" << static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << endl;

	  // ============== PostSmoothing ==============    
	  for (unsigned k = 0; k < npost; k++) {
	    solver_info = ( _VankaIsSet ) ? _LinSolver[ipde][ig]->solve(VankaIndex, _NSchurVar, _Schur) : _LinSolver[ipde][ig]->solve();
	  }
	}
	// ============== Update Solution ( _gridr-1 <= ig <= igridn-2 ) ==============
	for (unsigned ig = _gridr-1; ig < igridn-1; ig++) {
	  _solution[ig]->SumEpsToSol(_SolPdeIndex[ipde], _LinSolver[ipde][ig]->_EPS, _LinSolver[ipde][ig]->_RES, _LinSolver[ipde][ig]->KKoffset );	
	}
	// ============== Test for linear Convergence ==============
	if(test_linear && igridn>2) GetConvergence(pdename, igridn-2);
      }
      // ============== Update Solution ( ig = igridn )==============
      _solution[igridn-1]->SumEpsToSol(_SolPdeIndex[ipde], _LinSolver[ipde][igridn-1]->_EPS, _LinSolver[ipde][igridn-1]->_RES, _LinSolver[ipde][igridn-1]->KKoffset );
      // ============== Test for non-linear Convergence ==============
      bool conv = GetConvergence(pdename, igridn-1);
      if (conv == true) nonlinear_cycle = nonlinear_Vcycle_numeber + 1;
       
      cout << endl;
      cout << "COMPUTATION RESIDUAL: \t"<<static_cast<double>((clock()-start_time))/CLOCKS_PER_SEC << endl;

      cout << "CYCLE TIME:           \t"<<static_cast<double>((clock()-start_cycle_time))/CLOCKS_PER_SEC << endl;
    }
    // ==============  Solution Prolongation ==============
    if (igridn < _gridn) {
      ProlongatorSol(pdename, igridn);
    }
  }

  for (int ig = _gridr-1; ig < _gridn-1; ig++) {
    _LinSolver[ipde][ig]->_CC_flag = 0;
  }
    
  cout << "SOLVER TIME:   \t\t\t"<<static_cast<double>((clock()-start_mg_time))/CLOCKS_PER_SEC << endl;

}



//---------------------------------------------------------------------------------------------------------------
bool NonLinearMultiLevelProblem::GetConvergence(const char pdename[], const unsigned gridn) {

  unsigned ipde=GetPdeIndex(pdename);
  bool conv=true;
  double ResMax;
  double L2normEps;
  
  //for debugging purpose
  for (unsigned k=0; k<_SolPdeIndex[ipde].size(); k++) {
    unsigned indexSol=_SolPdeIndex[ipde][k];
    
    L2normEps    = _solution[gridn]->_Eps[indexSol]->l2_norm();
    ResMax       = _solution[gridn]->_Res[indexSol]->linfty_norm();

    cout << "level=" << _msh[gridn]->GetGridNumber() << "\tLinftynormRes" << SolName[indexSol] << "=" << ResMax    <<endl;
    cout << "level=" << _msh[gridn]->GetGridNumber() << "\tL2normEps"     << SolName[indexSol] << "=" << L2normEps <<endl;
    
    if (L2normEps < _non_linear_toll && conv==true) {
      conv=true;
    } 
    else {
      conv=false;
    }
  }
  return conv;
}

//*****************************************************************************************
int NonLinearMultiLevelProblem::FreeMultigrid() {
  PetscErrorCode ierr;
  for (int igridn=0; igridn<_gridn; igridn++) {
    for (int itype=0; itype<3; itype++) {
      for (int jtype=0; jtype<3; jtype++) {
 	delete ProlQitoQj_[itype][jtype][igridn];
      }
    }
    _solution[igridn]->FreeSolutionVectors();
  }
  return 1;
}

//*******************************************************************************************
void NonLinearMultiLevelProblem::AddSolution(const char name[], const char order[],
					     unsigned tmorder, const bool &Pde_type) {
  
  unsigned n=SolType.size();
  SolType.resize(n+1u);
  SolName.resize(n+1u);
  BdcType.resize(n+1u);
  SolTmorder.resize(n+1u);
  TestIfPressure.resize(n+1u);

  if (tmorder==0) {
    if (_test_time==0) {
      tmorder=1;
    } else {
      tmorder=2;
    }
  }

  TestIfPressure[n]=0;
  if (!strcmp(order,"linear")) {
    SolType[n]=0;
  } else if (!strcmp(order,"quadratic")) {
    SolType[n]=1;
  } else if (!strcmp(order,"biquadratic")) {
    SolType[n]=2;
  } else if (!strcmp(order,"constant")) {
    SolType[n]=3;
  } else if (!strcmp(order,"disc_linear")) {
    SolType[n]=4;
  } else {
    cout<<"Error! Invalid Finite Element entry for variable " << name << " in AddSolution function"<<endl;
    exit(0);
  }
  SolName[n]  = new char [8];
  BdcType[n]  = new char [20];
  strcpy(SolName[n],name);
  SolTmorder[n]=tmorder;
 
  cout << "Add variable " << std::setw(3) << SolName[n] << " discretized with FE type "
       << std::setw(12) << order << " and time discretzation order " << tmorder-1 << endl;

  for (unsigned ig=0; ig<_gridn; ig++) {
    _solution[ig]->AddSolution(name,order,tmorder,Pde_type);
  }
}

// *******************************************************
void NonLinearMultiLevelProblem::AssociatePropertyToSolution(const char solution_name[], const char solution_property[]){
  unsigned index=GetIndex(solution_name);
  if( !strcmp(solution_property,"pressure") || !strcmp(solution_property,"Pressure") ) TestIfPressure[index]=1;
  else if( !strcmp(solution_property,"default") || !strcmp(solution_property,"Default") ) TestIfPressure[index]=0;
  else {
    cout<<"Error invalid property in function NonLinearMultiLevelProblem::AssociatePropertyToSolution"<<endl;
    exit(0);
  }
}

// *******************************************************
void NonLinearMultiLevelProblem::Initialize(const char name[]) {
  
  unsigned i_start;
  unsigned i_end;
  if (!strcmp(name,"All")) {
    i_start=3;
    i_end=SolType.size();
  } else {
    i_start=GetIndex(name);
    i_end=i_start+1u;
  }

  unsigned indX=GetIndex("X");
  unsigned indY=GetIndex("Y");
  unsigned indZ=GetIndex("Z");
  double xx=0;
  double yy=0;
  double zz=0;
  
  
  double value;
  for (unsigned i=i_start; i<i_end; i++) {
    //CheckVectorSize(i);
    unsigned sol_type = SolType[i];
    for (unsigned ig=0; ig<_gridn; ig++) {
      unsigned num_el = _msh[ig]->GetElementNumber();
      _solution[ig]->ResizeSolutionVector(SolName[i]);
      if (ig>0) BuildProlongatorMatrix(ig,i);     
      //for parallel
      for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
        for (int iel=_msh[ig]->IS_Mts2Gmt_elem_offset[isdom]; 
	     iel < _msh[ig]->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
	  unsigned kel_gmt = _msh[ig]->IS_Mts2Gmt_elem[iel];    
	  unsigned sol_ord = _msh[ig]->GetEndIndex(SolType[i]);
	  unsigned nloc_dof= _msh[ig]->el->GetElementDofNumber(kel_gmt,sol_ord);
	  if(sol_type<3) {
            for(int j=0; j<nloc_dof; j++) {
	      unsigned inode=(sol_type<3)?(_msh[ig]->el->GetElementVertexIndex(kel_gmt,j)-1u):(kel_gmt+j*num_el);
	      int inode_Metis=_msh[ig]->GetMetisDof(inode,sol_type);
	      unsigned icoord_Metis=_msh[ig]->GetMetisDof(inode,2);
	      xx=(*_solution[ig]->_Sol[indX])(icoord_Metis);  
	      yy=(*_solution[ig]->_Sol[indY])(icoord_Metis);
	      zz=(*_solution[ig]->_Sol[indZ])(icoord_Metis);
	      
	      value = (sol_type<3)?_InitVariableFunction(xx,yy,zz,SolName[i]):0;
	      _solution[ig]->_Sol[i]->set(inode_Metis,value);
	      if (SolTmorder[i]==2) {
		_solution[ig]->_SolOld[i]->set(inode_Metis,value);
	      }
	    }
	  }
	}
      }
      _solution[ig]->_Sol[i]->close();
      if (SolTmorder[i]==2) {
	_solution[ig]->_SolOld[i]->close();
      }
    }
  }
}
  
  
void NonLinearMultiLevelProblem::AddPde(const char pdename[]){
  unsigned n=_PdeName.size();
  _PdeName.resize(n+1u);
  _PdeName[n]  = new char [30];
  strcpy(_PdeName[n],pdename);
}
  
// *******************************************************************  
unsigned NonLinearMultiLevelProblem::GetPdeIndex(const char pdename[]) const{
  unsigned index=0;
  while (strcmp(_PdeName[index],pdename)) {
    index++;
    if (index==_PdeName.size()) {
      cout<<"error! invalid Pde name ("<< pdename<<") in function GetPdeIndex(...)"<<endl;
      exit(0);
    }
  }
  return index; 
}
// *******************************************************
unsigned NonLinearMultiLevelProblem::GetIndex(const char name[]) const {
  unsigned index=0;
  while (strcmp(SolName[index],name)) {
    index++;
    if (index==SolType.size()) {
      cout<<"error! invalid solution name "<< name <<"entry GetIndex(...)"<<endl;
      exit(0);
    }
  }
  return index;
}

// *******************************************************
unsigned NonLinearMultiLevelProblem::GetSolType(const char name[]) {
  unsigned index=0;
  while (strcmp(SolName[index],name)) {
    index++;
    if (index==SolType.size()) {
      cout<<"error! invalid name entry GetSolType(...)"<<endl;
      exit(0);
    }
  }
  return SolType[index];
}

// *******************************************************
void NonLinearMultiLevelProblem::ClearSolPdeIndex() {
  for(unsigned i=0;i<_SolPdeIndex.size();i++){
    _SolPdeIndex[i].clear();
  }
  _SolPdeIndex.clear();
  _PdeIndex.clear(); 
}

// *******************************************************
void NonLinearMultiLevelProblem::AddSolutionToSolPdeIndex( const char pdename[], const char solname[]){
  int ipde=-1;
  for(unsigned i=0;i<_PdeIndex.size();i++){
    if(!strcmp(_PdeName[_PdeIndex[i]],pdename)) {
      ipde=i;
      break;
    }
  }
  if(ipde==-1){
    ipde=_PdeIndex.size();
    _PdeIndex.resize(ipde+1u);
    _PdeIndex[ipde]=GetPdeIndex(pdename);
    _SolPdeIndex.resize(ipde+1u);
  }
  unsigned jsol=0;
  for(unsigned j=0;j<_SolPdeIndex[ipde].size();j++){
    if(strcmp(SolName[_SolPdeIndex[ipde][j]],solname)) jsol++;
  }
  if(jsol==_SolPdeIndex[ipde].size()){
    _SolPdeIndex[ipde].resize(jsol+1u);
    _SolPdeIndex[ipde][jsol]=GetIndex(solname);
  }
}

// *******************************************************
unsigned NonLinearMultiLevelProblem::GetSolPdeIndex(const char pdename[], const char solname[]) {
  unsigned ipde=GetPdeIndex(pdename);  
  unsigned index=0;
  while (strcmp(SolName[_SolPdeIndex[ipde][index]],solname)) {
    index++;
    if (index==_SolPdeIndex[ipde].size()) {
      cout<<"error! invalid name entry NonLinearMultiLevelProblem::GetSolPdeIndex(const char pdename[], const char solname[])"<<endl;
      exit(0);
    }
  }
  return index;
}

// *******************************************************
void NonLinearMultiLevelProblem::ClearVankaIndex() {
  VankaIndex.clear();
};

// *******************************************************
void NonLinearMultiLevelProblem::AddToVankaIndex(const char pdename[], const char solname[]) {
  unsigned ipde=GetPdeIndex(pdename);
  unsigned n=VankaIndex.size();
  VankaIndex.resize(n+1u);
  unsigned varind=GetIndex(solname);

  for (unsigned i=0; i<_SolPdeIndex[ipde].size(); i++) {
    if (_SolPdeIndex[ipde][i]==varind) {
      VankaIndex[n]=i;
      break;
    }
    if (i==_SolPdeIndex[ipde].size()-1u) {
      cout<<"Error! the Vanka variable "<<solname<<" is not included in the Pde variables."<<endl;
      exit(0);
    }
  }
};

// *******************************************************


void NonLinearMultiLevelProblem::Restrictor(const unsigned &ipde, const unsigned &gridf, const unsigned &gridn, 
					    const unsigned &non_linear_iteration, const unsigned &linear_iteration, const bool &full_cycle){
    
  _LinSolver[ipde][gridf-1u]->SetEpsZero();
  _LinSolver[ipde][gridf-1u]->SetResZero();
   
  bool assemble_matrix = (linear_iteration == 0) ? true : false;  //Be carefull!!!! this is needed in the _assemble_function      
  if (gridf>=_gridr) {
    _assemble_function(*this,gridf-1,gridn-1u,ipde, assemble_matrix);
  }
  
  bool matrix_reuse=true;
  if(assemble_matrix){
    if (gridf>=_gridr) {
      if (!_LinSolver[ipde][gridf-1]->_CC_flag) {
	_LinSolver[ipde][gridf-1]->_CC_flag=1;
	_LinSolver[ipde][gridf-1]->_CC->matrix_PtAP(*_LinSolver[ipde][gridf]->_PP,*_LinSolver[ipde][gridf]->_KK,!matrix_reuse);
      } 
      else{
	_LinSolver[ipde][gridf-1]->_CC->matrix_PtAP(*_LinSolver[ipde][gridf]->_PP,*_LinSolver[ipde][gridf]->_KK,matrix_reuse);
      }
      _LinSolver[ipde][gridf-1u]->_KK->matrix_add(1.,*_LinSolver[ipde][gridf-1u]->_CC,"subset_nonzero_pattern");
    } 
    else { //Projection of the Matrix on the lower level
      if (non_linear_iteration==0 && ( full_cycle*(gridf==gridn-1u) || !full_cycle )) {
	_LinSolver[ipde][gridf-1]->_KK->matrix_PtAP(*_LinSolver[ipde][gridf]->_PP,*_LinSolver[ipde][gridf]->_KK,!matrix_reuse);
      }
      else{ 
	_LinSolver[ipde][gridf-1]->_KK->matrix_PtAP(*_LinSolver[ipde][gridf]->_PP,*_LinSolver[ipde][gridf]->_KK,matrix_reuse);
      }	    
    }
  }
      
  _LinSolver[ipde][gridf-1u]->_RESC->matrix_mult_transpose(*_LinSolver[ipde][gridf]->_RES, *_LinSolver[ipde][gridf]->_PP);
  *_LinSolver[ipde][gridf-1u]->_RES += *_LinSolver[ipde][gridf-1u]->_RESC;
}

// *******************************************************
int NonLinearMultiLevelProblem::Prolongator(const unsigned &ipde, const unsigned &gridf) {
  _LinSolver[ipde][gridf]->_EPSC->matrix_mult(*_LinSolver[ipde][gridf-1]->_EPS,*_LinSolver[ipde][gridf]->_PP);
  _LinSolver[ipde][gridf]->UpdateResidual();
  _LinSolver[ipde][gridf]->SumEpsCToEps();
  return 1;
}

// *******************************************************
void NonLinearMultiLevelProblem::ProlongatorSol(const char pdename[], unsigned gridf) {
  
  unsigned ipde = GetPdeIndex(pdename);
  
  for (unsigned k=0; k<_SolPdeIndex[ipde].size(); k++) {
    unsigned SolIndex=_SolPdeIndex[ipde][k];
    unsigned Typeindex=SolType[SolIndex];
    _solution[gridf]->_Sol[SolIndex]->matrix_mult(*_solution[gridf-1]->_Sol[SolIndex],*_solution[gridf]->_ProjMat[Typeindex]);
    _solution[gridf]->_Sol[SolIndex]->close(); 
    
  }
}

//---------------------------------------------------------------------------------------------------
// This routine generates the matrix for the projection of the FE matrix to finer grids 
//---------------------------------------------------------------------------------------------------

int NonLinearMultiLevelProblem::BuildProlongatorMatrix(unsigned gridf, const char pdename[]) {

  unsigned ipde = GetPdeIndex(pdename);
      
  if (gridf<1) {
    cout<<"Error! In function \"BuildProlongatorMatrix\" argument less then 1"<<endl;
    exit(0);
  }
  
  int ierr;
  
  int nf= _LinSolver[ipde][gridf]->KKIndex[_LinSolver[ipde][gridf]->KKIndex.size()-1u];
  int nc= _LinSolver[ipde][gridf-1]->KKIndex[_LinSolver[ipde][gridf-1]->KKIndex.size()-1u];
  int nf_loc = _LinSolver[ipde][gridf]->KKoffset[_LinSolver[ipde][gridf]->KKIndex.size()-1][_iproc]-_LinSolver[ipde][gridf]->KKoffset[0][_iproc];
  int nc_loc = _LinSolver[ipde][gridf-1]->KKoffset[_LinSolver[ipde][gridf-1]->KKIndex.size()-1][_iproc]-_LinSolver[ipde][gridf-1]->KKoffset[0][_iproc];
  _LinSolver[ipde][gridf]->_PP = SparseMatrix::build().release();
  _LinSolver[ipde][gridf]->_PP->init(nf,nc,nf_loc,nc_loc,27,27);
      
  for (unsigned k=0; k<_SolPdeIndex[ipde].size(); k++) {
    unsigned SolIndex=_SolPdeIndex[ipde][k];
       
    // loop on the coarse grid 
    for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
      for (int iel_mts=_msh[gridf-1]->IS_Mts2Gmt_elem_offset[isdom]; 
	   iel_mts < _msh[gridf-1]->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	unsigned iel = _msh[gridf-1]->IS_Mts2Gmt_elem[iel_mts];
	if(_msh[gridf-1]->el->GetRefinedElementIndex(iel)){ //only if the coarse element has been refined
    
	  short unsigned ielt=_msh[gridf-1]->el->GetElementType(iel);
	  type_elem[ielt][SolType[SolIndex]]->prolongation(*_LinSolver[ipde][gridf],*_LinSolver[ipde][gridf-1],iel,
							   _LinSolver[ipde][gridf]->_PP,SolIndex,k);
	
	}
      }
    }
  }
  _LinSolver[ipde][gridf]->_PP->close();
    
  return ierr;
}


//---------------------------------------------------------------------------------------------------
// This routine generates the matrix for the projection of the solution to finer grids  
//---------------------------------------------------------------------------------------------------

void NonLinearMultiLevelProblem::BuildProlongatorMatrix(unsigned gridf, unsigned SolIndex) {
  
  if (gridf<1) {
    cout<<"Error! In function \"BuildProlongatorMatrix\" argument less then 1"<<endl;
    exit(0);
  }
  
  unsigned TypeIndex=SolType[SolIndex];
    
  if(_solution[gridf]->_ProjMatFlag[TypeIndex]==0){
    _solution[gridf]->_ProjMatFlag[TypeIndex]=1;

    int nf     = _msh[gridf]->MetisOffset[SolType[SolIndex]][_nprocs];
    int nc     = _msh[gridf-1]->MetisOffset[SolType[SolIndex]][_nprocs];
    int nf_loc = _msh[gridf]->own_size[SolType[SolIndex]][_iproc];
    int nc_loc = _msh[gridf-1]->own_size[SolType[SolIndex]][_iproc]; 

    _solution[gridf]->_ProjMat[TypeIndex] = SparseMatrix::build().release();
    _solution[gridf]->_ProjMat[TypeIndex]->init(nf,nc,nf_loc,nc_loc,27,27);
 
    // loop on the coarse grid 
    for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
      for (int iel_mts=_msh[gridf-1]->IS_Mts2Gmt_elem_offset[isdom]; 
	   iel_mts < _msh[gridf-1]->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	unsigned iel = _msh[gridf-1]->IS_Mts2Gmt_elem[iel_mts];
	if(_msh[gridf-1]->el->GetRefinedElementIndex(iel)){ //only if the coarse element has been refined
	  short unsigned ielt=_msh[gridf-1]->el->GetElementType(iel);
	  type_elem[ielt][SolType[SolIndex]]->prolongation(*_msh[gridf],*_msh[gridf-1],iel,
							   _solution[gridf]->_ProjMat[TypeIndex]); 
	}
      }
    }
    _solution[gridf]->_ProjMat[TypeIndex]->close();
  }
}

//---------------------------------------------------------------------------------------------------
// This routine generates the matrices for the projection of the solutions from different FE spaces 
//---------------------------------------------------------------------------------------------------

void NonLinearMultiLevelProblem::BuildProlongatorMatrices() {

  ProlQitoQj_[0][0].resize(_gridn);
  ProlQitoQj_[0][1].resize(_gridn);
  ProlQitoQj_[0][2].resize(_gridn);
  
  ProlQitoQj_[1][0].resize(_gridn);
  ProlQitoQj_[1][1].resize(_gridn);
  ProlQitoQj_[1][2].resize(_gridn);
  
  ProlQitoQj_[2][0].resize(_gridn);
  ProlQitoQj_[2][1].resize(_gridn);
  ProlQitoQj_[2][2].resize(_gridn);
  
  for (unsigned igridn=0; igridn<_gridn; igridn++) {
    for(int itype=0;itype<3;itype++){
      int ni = _msh[igridn]->MetisOffset[itype][_nprocs];
      bool *testnode=new bool [ni];
      for (int jtype=0; jtype<3; jtype++) {
        int nj = _msh[igridn]->MetisOffset[jtype][_nprocs];
	memset(testnode,0,ni*sizeof(bool));
	ProlQitoQj_[itype][jtype][igridn] = SparseMatrix::build().release();
	ProlQitoQj_[itype][jtype][igridn]->init(ni,nj,_msh[igridn]->own_size[itype][_iproc],
						_msh[igridn]->own_size[jtype][_iproc],27,27);
		
	for(int isdom=_iproc; isdom<_iproc+1; isdom++) {
	  for (int iel_mts=_msh[igridn]->IS_Mts2Gmt_elem_offset[isdom]; 
	       iel_mts < _msh[igridn]->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
	    unsigned iel = _msh[igridn]->IS_Mts2Gmt_elem[iel_mts];
	    short unsigned ielt=_msh[igridn]->el->GetElementType(iel);
            type_elem[ielt][jtype]->ProlQitoQj(*_msh[igridn],iel,ProlQitoQj_[itype][jtype][igridn],testnode,itype);	  
	  }
	}
	ProlQitoQj_[itype][jtype][igridn]->close();
      }
      delete [] testnode;
    }
  }
}

// *******************************************************************
void NonLinearMultiLevelProblem::GenerateBdc(const char name[], const char bdc_type[]) {

  if(_bdc_func_set==false) {
    cout << "Error: The boundary condition user-function is not set! Please call the AttachSetBoundaryConditionFunction routine" 
         << endl;
    exit(1); 
  }
  
  const short unsigned NV1[6][2]= {{9,9},{6,6},{9,6},{3,3},{3,3},{1,1}};
  
  unsigned indX=GetIndex("X");
  unsigned indY=GetIndex("Y");
  unsigned indZ=GetIndex("Z");
  double xx;
  double yy;
  double zz;

  unsigned i_start;
  unsigned i_end;
  if (!strcmp(name,"All")) {
    i_start=0;
    i_end=SolType.size();
    for (unsigned k=i_start; k<i_end; k++) {
      if(_solution[0]->_ResEpsBdcFlag[k]){
	sprintf(BdcType[k],"Steady");
	cout << "Set " << std::setw(15) << BdcType[k] << " Boundary_condition"
	     << " for variable " << std::setw(3) << SolName[k] << endl;
      }
      else {
	sprintf(BdcType[k],"Not-available");
      }
	
    }
  } 
  else {
    i_start=GetIndex(name);
    i_end=i_start+1u;
    if(_solution[0]->_ResEpsBdcFlag[i_start]){
      if (!strcmp(bdc_type,"Steady")) {
	strcpy(BdcType[i_start],bdc_type);
      } else if (!strcmp(bdc_type,"Time_dependent")) {
	strcpy(BdcType[i_start],bdc_type);
      } else {
	cout << "Error! Invalid boundary condition specified for " << SolName[i_start]
	     << " in GenerateBdc function" << endl;
	exit(1);
      }
      cout << "Set " << std::setw(14) <<BdcType[i_start] << " Boundary_condition"
	   << " for variable " << std::setw(3) << SolName[i_start] << endl;
    }
    else {
      sprintf(BdcType[i_start],"Not-available");
    }
  }
  // 2 Default Neumann
  // 1 DD Dirichlet
  // 0 Dirichlet
  for (unsigned igridn=0; igridn<_gridn; igridn++) {
    for (unsigned i=i_start; i<i_end; i++) {
      if(_solution[igridn]->_ResEpsBdcFlag[i]){
	for (unsigned j=_msh[igridn]->MetisOffset[SolType[i]][_iproc]; j<_msh[igridn]->MetisOffset[SolType[i]][_iproc+1]; j++) {
	  _solution[igridn]->_Bdc[i]->set(j,2.);
	}
	if (SolType[i]<3) {  
	  for(int isdom=_iproc; isdom<_iproc+1; isdom++) {   // 1 DD Dirichlet
	    for (int iel=_msh[igridn]->IS_Mts2Gmt_elem_offset[isdom]; 
		 iel < _msh[igridn]->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
	      unsigned iel_gmt = _msh[igridn]->IS_Mts2Gmt_elem[iel];
	      for (unsigned jface=0; jface<_msh[igridn]->el->GetElementFaceNumber(iel_gmt); jface++) {
		if (_msh[igridn]->el->GetFaceElementIndex(iel_gmt,jface)==0) { //Domain Decomposition Dirichlet
		  short unsigned ielt=_msh[igridn]->el->GetElementType(iel_gmt);
		  unsigned nv1=(!TestIfPressure[i])?
		    NV1[ielt][jface<_msh[igridn]->el->GetElementFaceNumber(iel_gmt,0)]: //non-pressure type
		    _msh[igridn]->el->GetElementDofNumber(iel,_msh[igridn]->GetEndIndex(SolType[i]) ); //pressure type
		  for (unsigned iv=0; iv<nv1; iv++) {
		    unsigned inode=(!TestIfPressure[i])? 
		      _msh[igridn]->el->GetFaceVertexIndex(iel_gmt,jface,iv)-1u: //non-pressure type
		      _msh[igridn]->el->GetElementVertexIndex(iel_gmt,iv)-1u;    //pressure type
		    unsigned inode_Metis=_msh[igridn]->GetMetisDof(inode,SolType[i]);
		    _solution[igridn]->_Bdc[i]->set(inode_Metis,1.);
		  }
		}
	      }
	    }
	  }
	  for(int isdom=_iproc; isdom<_iproc+1; isdom++) {  // 0 Dirichlet
	    for (int iel=_msh[igridn]->IS_Mts2Gmt_elem_offset[isdom]; 
		 iel < _msh[igridn]->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
	      unsigned iel_gmt = _msh[igridn]->IS_Mts2Gmt_elem[iel];
	      for (unsigned jface=0; jface<_msh[igridn]->el->GetElementFaceNumber(iel_gmt); jface++) {
		if (_msh[igridn]->el->GetFaceElementIndex(iel_gmt,jface)<0) { //Dirichlet
		  short unsigned ielt=_msh[igridn]->el->GetElementType(iel_gmt);
		  unsigned nv1=NV1[ielt][jface<_msh[igridn]->el->GetElementFaceNumber(iel_gmt,0)];
		  for (unsigned iv=0; iv<nv1; iv++) {
		    unsigned inode=_msh[igridn]->el->GetFaceVertexIndex(iel_gmt,jface,iv)-1u;
		    unsigned inode_coord_Metis=_msh[igridn]->GetMetisDof(inode,2);
		    double value;
		    xx=(*_solution[igridn]->_Sol[indX])(inode_coord_Metis);  
		    yy=(*_solution[igridn]->_Sol[indY])(inode_coord_Metis);
		    zz=(*_solution[igridn]->_Sol[indZ])(inode_coord_Metis);
		    bool test=_SetBoundaryConditionFunction(xx,yy,zz,SolName[i],value,-(_msh[igridn]->el->GetFaceElementIndex(iel_gmt,jface)+1),_time);
		    if (test) {
		      unsigned inode_Metis=_msh[igridn]->GetMetisDof(inode,SolType[i]);
		      _solution[igridn]->_Bdc[i]->set(inode_Metis,0.);
		      _solution[igridn]->_Sol[i]->set(inode_Metis,value);
		    }
		  }
		}
	      }
	    }
	  }
	}
	else if(TestIfPressure[i]){ // 1 DD Dirichlet for pressure types only
	  for(int isdom=_iproc; isdom<_iproc+1; isdom++) {   
	    unsigned nel=_msh[igridn]->GetElementNumber();
	    for (int iel=_msh[igridn]->IS_Mts2Gmt_elem_offset[isdom]; 
		 iel < _msh[igridn]->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
	      unsigned iel_gmt = _msh[igridn]->IS_Mts2Gmt_elem[iel];
	      for (unsigned jface=0; jface<_msh[igridn]->el->GetElementFaceNumber(iel_gmt); jface++) {
		if (_msh[igridn]->el->GetFaceElementIndex(iel_gmt,jface)==0) { //Domain Decomposition Dirichlet
		  short unsigned ielt=_msh[igridn]->el->GetElementType(iel_gmt);
		  unsigned nv1=_msh[igridn]->el->GetElementDofNumber(iel_gmt,_msh[igridn]->GetEndIndex(SolType[i]));
	  	  for (unsigned iv=0; iv<nv1; iv++) {
		    unsigned inode=(iel_gmt+iv*nel);
		    unsigned inode_Metis=_msh[igridn]->GetMetisDof(inode,SolType[i]);
		    _solution[igridn]->_Bdc[i]->set(inode_Metis,1.);
		  }
		}
	      }
	    }
	  }
	}	
	_solution[igridn]->_Sol[i]->close();
	_solution[igridn]->_Bdc[i]->close();
      }
    }
  }
}



// *******************************************************************
int  NonLinearMultiLevelProblem::printsol_gmv_binary(const char type[],unsigned igridn, bool debug) const {
 
  if (igridn==0) igridn=_gridn;
  
  unsigned igridr=(_gridr <= igridn)?_gridr:igridn;

  // ********** linear -> index==0 *** quadratic -> index==1 **********
  unsigned index=(strcmp(type,"linear"))?1:0;

  char *filename = new char[60];
  sprintf(filename,"./output/mesh.level%d.%d.%s.gmv",igridn,_time_step,type);

  std::ofstream fout;
  
  if(_iproc!=0) {
    fout.rdbuf();   //redirect to dev_null
  }
  else {
    fout.open(filename);
    if (!fout) {
      cout << "Output mesh file "<<filename<<" cannot be opened.\n";
      exit(0);
    }
  }
  
  PetscErrorCode ierr;

  unsigned nvt=0;
  unsigned nvt_max=0;
  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
    unsigned nvt_ig=_msh[ig]->MetisOffset[index][_nprocs];
    nvt_max=(nvt_max>nvt_ig)?nvt_max:nvt_ig;
    nvt+=nvt_ig;
  }
  
  double *var_nd=new double [nvt_max+1]; //TO FIX Valgrind complaints! In reality it should be only nvt
  vector <NumericVector*> Mysol(igridn);
  for(unsigned ig=igridr-1u; ig<_gridn; ig++) {
    Mysol[ig] = NumericVector::build().release();
    Mysol[ig]->init(_msh[ig]->MetisOffset[index][_nprocs],_msh[ig]->own_size[index][_iproc],true,AUTOMATIC);
  }
     
  // ********** Header **********
  char *det= new char[10];
  sprintf(det,"%s","gmvinput");
  fout.write((char *)det,sizeof(char)*8);
  sprintf(det,"%s","ieeei4r8");
  fout.write((char *)det,sizeof(char)*8);

  // ********** Start printing node coordinates  **********
  sprintf(det,"%s","nodes");
  fout.write((char *)det,sizeof(char)*8);
  fout.write((char *)&nvt,sizeof(unsigned));
    
  unsigned indXYZ[3];
  indXYZ[0]=GetIndex("X");
  indXYZ[1]=GetIndex("Y");
  indXYZ[2]=GetIndex("Z");
   
  for (int i=0; i<3; i++) {
    for (unsigned ig=igridr-1u; ig<igridn; ig++) {
      Mysol[ig]->matrix_mult(*_solution[ig]->_Sol[indXYZ[i]],*ProlQitoQj_[index][SolType[indXYZ[i]]][ig]);
      vector <double> v_local;
      Mysol[ig]->localize_to_one(v_local,0);
      unsigned nvt_ig=_msh[ig]->MetisOffset[index][_nprocs];      
      if(_iproc==0){ 
	for (unsigned ii=0; ii<nvt_ig; ii++) 
	  var_nd[ii]= v_local[ii];
      }
      if (_moving_mesh) {
	unsigned indDXDYDZ=GetIndex(_moving_vars[i].c_str());
	Mysol[ig]->matrix_mult(*_solution[ig]->_Sol[indDXDYDZ],*ProlQitoQj_[index][SolType[indDXDYDZ]][ig]);
	Mysol[ig]->localize_to_one(v_local,0);
	unsigned nvt_ig=_msh[ig]->MetisOffset[index][_nprocs];      
	if(_iproc==0){ 
	  for (unsigned ii=0; ii<nvt_ig; ii++) 
	    var_nd[ii]+= v_local[ii];
	}
      }
      fout.write((char *)&var_nd[0],nvt_ig*sizeof(double)); 
    }
  }
  // ********** End printing node coordinates  **********

  // ********** Start printing cell connectivity  **********
  const int eltp[2][6]= {{8,4,6,4,3,2},{20,10,15,8,6,3}};
  sprintf(det,"%s","cells");
  fout.write((char *)det,sizeof(char)*8);

  unsigned nel=0;
  for (unsigned ig=igridr-1u; ig<igridn-1u; ig++)
    nel+=( _msh[ig]->GetElementNumber() - _msh[ig]->el->GetRefinedElementNumber());
  nel+=_msh[igridn-1u]->GetElementNumber();
  fout.write((char *)&nel,sizeof(unsigned));

  unsigned topology[27];
  unsigned offset=1;
  
  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
    for (unsigned ii=0; ii<_msh[ig]->GetElementNumber(); ii++) {
      if ( ig==igridn-1u || 0==_msh[ig]->el->GetRefinedElementIndex(ii)) {
        short unsigned ielt=_msh[ig]->el->GetElementType(ii);
        if (ielt==0) sprintf(det,"phex%d",eltp[index][0]);
        else if (ielt==1) sprintf(det,"ptet%d",eltp[index][1]);
        else if (ielt==2) sprintf(det,"pprism%d",eltp[index][2]);
        else if (ielt==3) {
          if (eltp[index][3]==8) sprintf(det,"%dquad",eltp[index][3]);
          else sprintf(det,"quad");
        } else if (ielt==4) {
          if (eltp[index][4]==6) sprintf(det,"%dtri",eltp[index][4]);
          else sprintf(det,"tri");
        } else if (ielt==5) {
          if (eltp[index][5]==3) sprintf(det,"%dline",eltp[index][5]);
          else sprintf(det,"line");
        }
        fout.write((char *)det,sizeof(char)*8);
        fout.write((char *)&NVE[ielt][index],sizeof(unsigned));
	for(unsigned j=0;j<NVE[ielt][index];j++){
	  
	  unsigned jnode=_msh[ig]->el->GetElementVertexIndex(ii,j)-1u;
	  unsigned jnode_Metis = _msh[ig]->GetMetisDof(jnode,index);
	  	  
	  topology[j]=jnode_Metis+offset;
	}
	fout.write((char *)topology,sizeof(unsigned)*NVE[ielt][index]);
      }
    }
    offset+=_msh[ig]->MetisOffset[index][_nprocs];
  }
  // ********** End printing cell connectivity  **********
  
  double *var_el=new double [nel+1]; //TO FIX Valgrind complaints! In reality it should be only nel
  
  // ********** Start printing Variables **********
  const unsigned zero=0u;
  const unsigned one=1u;
  sprintf(det,"%s","variable");
  fout.write((char *)det,sizeof(char)*8);

  // ********** Start printing Regions **********
  strcpy(det,"Regions");
  fout.write((char *)det,sizeof(char)*8);
  fout.write((char *)&zero,sizeof(unsigned));

  int icount=0;
  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
    for (unsigned ii=0; ii<_msh[ig]->GetElementNumber(); ii++) {
      if ( ig==igridn-1u || 0==_msh[ig]->el->GetRefinedElementIndex(ii)) {
	var_el[icount]=_msh[ig]->el->GetElementGroup(ii);
        icount++;
      }
    }
  }
  fout.write((char *)&var_el[0],nel*sizeof(double));
  
  if(_nprocs>=1){
    strcpy(det,"Reg_proc");
    fout.write((char *)det,sizeof(char)*8);
    fout.write((char *)&zero,sizeof(unsigned));

    int icount=0;
    for (unsigned ig=igridr-1u; ig<igridn; ig++) {
      for (unsigned ii=0; ii<_msh[ig]->GetElementNumber(); ii++) {
	if ( ig==igridn-1u || 0==_msh[ig]->el->GetRefinedElementIndex(ii)) {
	  var_el[icount]=_msh[ig]->epart[ii];
	  icount++;
	}
      }
    }
    fout.write((char *)&var_el[0],nel*sizeof(double));
  }
  
  // ********** End printing Regions **********
  
  // ********** Start printing Solution **********
  for (unsigned i=0; i<SolType.size(); i++) {
    for(int name=0;name<4;name++){
      if (name==0){
	sprintf(det,"%s %s","Sol",SolName[i]);
      }
      else if (name==1){
	sprintf(det,"%s %s","Bdc",SolName[i]);
      }
      else if (name==2){
	sprintf(det,"%s %s","Res",SolName[i]);
      }
      else{
	sprintf(det,"%s %s","Eps",SolName[i]);
      }
      if(name==0 || (debug && _solution[igridn-1u]->_ResEpsBdcFlag[i])){
	if (SolType[i]<3) {  // **********  on the nodes **********
	  fout.write((char *)det,sizeof(char)*8);
	  fout.write((char *)&one,sizeof(unsigned));
	  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
	    if (name==0){
	      Mysol[ig]->matrix_mult(*_solution[ig]->_Sol[i],*ProlQitoQj_[index][SolType[i]][ig]);  
	    }
	    else if (name==1){
	      Mysol[ig]->matrix_mult(*_solution[ig]->_Bdc[i],*ProlQitoQj_[index][SolType[i]][ig]);
	    }
	    else if (name==2){
	      Mysol[ig]->matrix_mult(*_solution[ig]->_Res[i],*ProlQitoQj_[index][SolType[i]][ig]);
	    }
	    else{
	      Mysol[ig]->matrix_mult(*_solution[ig]->_Eps[i],*ProlQitoQj_[index][SolType[i]][ig]);
	    }
	    std::vector<double> v_local;
	    Mysol[ig]->localize_to_one(v_local,0);
	    fout.write((char *)&v_local[0],v_local.size()*sizeof(double));
	  }
	}
	else { // ********** on the elements **********
	  fout.write((char *)det,sizeof(char)*8);
	  fout.write((char *)&zero,sizeof(unsigned));
	  int icount=0;
	  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
	    std::vector<double> v_local;
	    if (name==0){
	      _solution[ig]->_Sol[i]->localize_to_one(v_local,0); 
	    }
	    if (name==1){
	      _solution[ig]->_Bdc[i]->localize_to_one(v_local,0); 
	    }
	    else if (name==2){
	      _solution[ig]->_Res[i]->localize_to_one(v_local,0);
	    }
	    else{
	      _solution[ig]->_Eps[i]->localize_to_one(v_local,0);
	    }
	    for (unsigned ii=0; ii<_msh[ig]->GetElementNumber(); ii++) {
	      if ( ig==igridn-1u || 0==_msh[ig]->el->GetRefinedElementIndex(ii)) {
		unsigned iel_Metis = _msh[ig]->GetMetisDof(ii,SolType[i]);
		var_el[icount]=v_local[iel_Metis];
		icount++;
	      }
	    }
	  }
	  fout.write((char *)&var_el[0],nel*sizeof(double));
	}	
      }
    }
  }
  // ********** End printing Solution **********
  sprintf(det,"%s","endvars");
  fout.write((char *)det,sizeof(char)*8);
  
  // ********** End printing Variables **********
  sprintf(det,"%s","endgmv");
  fout.write((char *)det,sizeof(char)*8);
  fout.close();
  // ********** End printing file **********
  
  // Free memory
  delete [] var_el;
  delete [] var_nd;
  for (unsigned ig=igridr-1u; ig<igridn; ig++) {
    delete Mysol[ig];
  }
  delete [] det;
  delete [] filename;
  return 1;
}



// *************************************************************************
// This function prints the solution in vtu binary (64-encoded) format
// *************************************************************************

void  NonLinearMultiLevelProblem::printsol_vtu_inline(const char type[], std::vector<std::string>& vars) const {
  
  bool test_all=!(vars[0].compare("All"));
  
  int icount;
  unsigned index=0;
  unsigned index_nd=0;
  if (!strcmp(type,"linear")) {   //linear
    index=0;
    index_nd=0;
  } else if (!strcmp(type,"quadratic")) { //quadratic
    index=1;
    index_nd=1;
  } else if (!strcmp(type,"biquadratic")) { //biquadratic
    index=3;
    index_nd=2;
  }

  const int eltp[4][6]= {{12,10,13,9,5,3},{25,24,26,23,22,21},{},{29,1,1,28,22,1}};
  
  char *filename= new char[60];
  sprintf(filename,"./output/mesh.level%d.%d.%s.vtu",_gridn,_time_step,type);
  std::ofstream fout;
  
  if(_iproc!=0) {
    fout.rdbuf();   //redirect to dev_null
  }
  else {
    fout.open(filename);
    if (!fout) {
      cout << "Output mesh file "<<filename<<" cannot be opened.\n";
      exit(0);
    }
  }
  // haed ************************************************
  fout<<"<?xml version=\"1.0\"?>" << endl;
  fout<<"<VTKFile type = \"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
  fout << " <UnstructuredGrid>" << endl;
  //-----------------------------------------------------------------------------------------------------------

  //----------------------------------------------------------------------------------------------------------
  
  unsigned nvt=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    unsigned nvt_ig=_msh[ig]->MetisOffset[index_nd][_nprocs];
    nvt+=nvt_ig;
  }  

  unsigned nel=0;
  unsigned counter=0;
  for (unsigned ig=_gridr-1u; ig<_gridn-1u; ig++) {
    nel+=( _msh[ig]->GetElementNumber() - _msh[ig]->el->GetRefinedElementNumber());
    counter+=(_msh[ig]->el->GetElementNumber("Hex")-_msh[ig]->el->GetRefinedElementNumber("Hex"))*NVE[0][index];
    counter+=(_msh[ig]->el->GetElementNumber("Tet")-_msh[ig]->el->GetRefinedElementNumber("Tet"))*NVE[1][index];
    counter+=(_msh[ig]->el->GetElementNumber("Wedge")-_msh[ig]->el->GetRefinedElementNumber("Wedge"))*NVE[2][index];
    counter+=(_msh[ig]->el->GetElementNumber("Quad")-_msh[ig]->el->GetRefinedElementNumber("Quad"))*NVE[3][index];
    counter+=(_msh[ig]->el->GetElementNumber("Triangle")-_msh[ig]->el->GetRefinedElementNumber("Triangle"))*NVE[4][index];
    counter+=(_msh[ig]->el->GetElementNumber("Line")-_msh[ig]->el->GetRefinedElementNumber("Line"))*NVE[5][index];
  }
  nel+=_msh[_gridn-1u]->GetElementNumber();
  counter+=_msh[_gridn-1u]->el->GetElementNumber("Hex")*NVE[0][index];
  counter+=_msh[_gridn-1u]->el->GetElementNumber("Tet")*NVE[1][index];
  counter+=_msh[_gridn-1u]->el->GetElementNumber("Wedge")*NVE[2][index];
  counter+=_msh[_gridn-1u]->el->GetElementNumber("Quad")*NVE[3][index];
  counter+=_msh[_gridn-1u]->el->GetElementNumber("Triangle")*NVE[4][index];
  counter+=_msh[_gridn-1u]->el->GetElementNumber("Line")*NVE[5][index]; 
  
  const unsigned dim_array_coord [] = { nvt*3*sizeof(float) };  
  const unsigned dim_array_conn[]   = { counter*sizeof(int) };
  const unsigned dim_array_off []   = { nel*sizeof(int) };
  const unsigned dim_array_type []  = { nel*sizeof(short unsigned) };
  const unsigned dim_array_reg []   = { nel*sizeof(short unsigned) };
  const unsigned dim_array_elvar [] = { nel*sizeof(float) };
  const unsigned dim_array_ndvar [] = { nvt*sizeof(float) };

  // initialize common buffer_void memory 
  unsigned buffer_size=(dim_array_coord[0]>dim_array_conn[0])?dim_array_coord[0]:dim_array_conn[0];
  void *buffer_void=new char [buffer_size];
  char *buffer_char=static_cast <char *>(buffer_void);
  
  size_t cch;
  cch = b64::b64_encode(&buffer_char[0], buffer_size , NULL, 0);  
  vector <char> enc;
  enc.resize(cch);
  char *pt_char;
  
  fout << "  <Piece NumberOfPoints= \"" << nvt << "\" NumberOfCells= \"" << nel << "\" >" << endl;
  
  //-----------------------------------------------------------------------------------------------
  // print coordinates *********************************************Solu*******************************************
  fout << "   <Points>" << endl;
  fout << "    <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">" << endl;
  
  unsigned indXYZ[3];
  indXYZ[0]=GetIndex("X");
  indXYZ[1]=GetIndex("Y");
  indXYZ[2]=GetIndex("Z");
  
  vector <NumericVector*> mysol(_gridn);
  for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    mysol[ig] = NumericVector::build().release();
    mysol[ig]->init(_msh[ig]->MetisOffset[index_nd][_nprocs],_msh[ig]->own_size[index_nd][_iproc],
		    true,AUTOMATIC);
  }
  
  // point pointer to common mamory area buffer of void type;
  float *var_coord= static_cast<float*>(buffer_void);

  unsigned offset_nvt3=0;
  for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    std::vector<double> v_local;
    unsigned nvt_ig=_msh[ig]->MetisOffset[index_nd][_nprocs];
    for(int kk=0;kk<3;kk++) {
      mysol[ig]->matrix_mult(*_solution[ig]->_Sol[indXYZ[kk]],*ProlQitoQj_[index_nd][SolType[indXYZ[kk]]][ig]);
      mysol[ig]->localize_to_one(v_local,0);
      if(_iproc==0) { 
	for (unsigned i=0; i<nvt_ig; i++) {
	  var_coord[offset_nvt3+i*3+kk] = v_local[i];
	}
      } //if iproc
    } //loop over dimension
    offset_nvt3+=3*nvt_ig;
  }
  
  if (_moving_mesh) {
    
    unsigned offset_nvt3=0;
    unsigned indDXDYDZ[3];
    indDXDYDZ[0]=GetIndex(_moving_vars[0].c_str());
    indDXDYDZ[1]=GetIndex(_moving_vars[1].c_str());
    if(_msh[0]->GetDimension() == 3) {
      indDXDYDZ[2]=GetIndex(_moving_vars[2].c_str());
    }
      
    for(unsigned ig=_gridr-1u; ig<_gridn; ig++){
      std::vector<double> v_local;
      unsigned nvt_ig=_msh[ig]->MetisOffset[index_nd][_nprocs];
      for(int kk=0;kk<_msh[0]->GetDimension();kk++) {
	mysol[ig]->matrix_mult(*_solution[ig]->_Sol[indDXDYDZ[kk]],*ProlQitoQj_[index_nd][SolType[indDXDYDZ[kk]]][ig]);
        mysol[ig]->localize_to_one(v_local,0);
	if(_iproc==0) { 
	  for (unsigned i=0; i<nvt_ig; i++) {
	    var_coord[offset_nvt3+i*3+kk] += v_local[i];
	  }
	} //if iproc
      } //loop over dimension
      offset_nvt3+=3*nvt_ig;
    }
  }
  
  if(_iproc==0) {
    //print coordinates dimension
    cch = b64::b64_encode(&dim_array_coord[0], sizeof(dim_array_coord), NULL, 0);  
    b64::b64_encode(&dim_array_coord[0], sizeof(dim_array_coord), &enc[0], cch);
    pt_char=&enc[0];
    for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
    //print coordinates array
    cch = b64::b64_encode(&var_coord[0], dim_array_coord[0] , NULL, 0);  
    b64::b64_encode(&var_coord[0], dim_array_coord[0], &enc[0], cch);
    pt_char=&enc[0];
    for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    fout << endl;
  }
  fout << "    </DataArray>" << endl;
  fout << "   </Points>" << endl;
  //-----------------------------------------------------------------------------------------------
  
  //-----------------------------------------------------------------------------------------------
  // Printing of element connectivity - offset - format type  *
  fout << "   <Cells>" << endl;
  
  //-----------------------------------------------------------------------------------------------
  //print connectivity
  fout << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">" << endl;
  
  // point pointer to common mamory area buffer of void type;
  int *var_conn = static_cast <int*> (buffer_void);
  icount = 0;
  unsigned offset_nvt=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned iel=0; iel<_msh[ig]->GetElementNumber(); iel++) {
      if (_msh[ig]->el->GetRefinedElementIndex(iel)==0 || ig==_gridn-1u) {
        for (unsigned j=0; j<_msh[ig]->el->GetElementDofNumber(iel,index); j++) {
	  unsigned jnode=_msh[ig]->el->GetElementVertexIndex(iel,j)-1u;
	  unsigned jnode_Metis = _msh[ig]->GetMetisDof(jnode,index_nd);
	  var_conn[icount] = offset_nvt+jnode_Metis;
	  icount++;
	}
      }
    }
    offset_nvt+=_msh[ig]->MetisOffset[index_nd][_nprocs];
  }
  
  //print connectivity dimension
  cch = b64::b64_encode(&dim_array_conn[0], sizeof(dim_array_conn), NULL, 0);  
  b64::b64_encode(&dim_array_conn[0], sizeof(dim_array_conn), &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  //print connectivity array
  cch = b64::b64_encode(&var_conn[0], dim_array_conn[0] , NULL, 0);  
  b64::b64_encode(&var_conn[0], dim_array_conn[0], &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  fout << endl;
  
  fout << "    </DataArray>" << endl;
  //------------------------------------------------------------------------------------------------
  
  //-------------------------------------------------------------------------------------------------
  //printing offset
  fout << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">" << endl;
 
  // point pointer to common mamory area buffer of void type;
  int *var_off=static_cast <int*>(buffer_void);
  icount = 0;
  int offset_el=0;
  //print offset array
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned iel=0; iel<_msh[ig]->GetElementNumber(); iel++) {
      if (_msh[ig]->el->GetRefinedElementIndex(iel)==0 || ig==_gridn-1u) {
	unsigned iel_Metis = _msh[ig]->GetMetisDof(iel,3);
        offset_el += _msh[ig]->el->GetElementDofNumber(iel_Metis,index);
        var_off[icount] = offset_el;
	icount++;
      }
    }
  }
  
  //print offset dimension
  cch = b64::b64_encode(&dim_array_off[0], sizeof(dim_array_off), NULL, 0);
  b64::b64_encode(&dim_array_off[0], sizeof(dim_array_off), &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  //print offset array
  cch = b64::b64_encode(&var_off[0], dim_array_off[0] , NULL, 0);  
  b64::b64_encode(&var_off[0], dim_array_off[0], &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  fout << endl;
  
  fout << "    </DataArray>" << endl;
  //--------------------------------------------------------------------------------------------------
  
  //--------------------------------------------------------------------------------------------------
  //Element format type : 23:Serendipity(8-nodes)  28:Quad9-Biquadratic
  fout << "    <DataArray type=\"UInt16\" Name=\"types\" format=\"binary\">" << endl;
   
  // point pointer to common mamory area buffer of void type;
  unsigned short *var_type = static_cast <unsigned short*> (buffer_void);
  icount=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned ii=0; ii<_msh[ig]->GetElementNumber(); ii++) {
      if (_msh[ig]->el->GetRefinedElementIndex(ii)==0 || ig==_gridn-1u) {
	unsigned iel_Metis = _msh[ig]->GetMetisDof(ii,3);
        short unsigned ielt=_msh[ig]->el->GetElementType(iel_Metis);
	var_type[icount] = (short unsigned)(eltp[index][ielt]);
	icount++;
      }
    }
  }
  
  //print element format dimension
  cch = b64::b64_encode(&dim_array_type[0], sizeof(dim_array_type), NULL, 0);  
  b64::b64_encode(&dim_array_type[0], sizeof(dim_array_type), &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  //print element format array
  cch = b64::b64_encode(&var_type[0], dim_array_type[0] , NULL, 0);  
  b64::b64_encode(&var_type[0], dim_array_type[0], &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  fout << endl;
  fout << "    </DataArray>" << endl;
  //----------------------------------------------------------------------------------------------------
  
  fout << "   </Cells>" << endl;
  //--------------------------------------------------------------------------------------------------

  
  
  // /Print Cell Data ****************************************************************************
  fout << "   <CellData Scalars=\"scalars\">" << endl;

  //--------------------------------------------------------------------------------------------
  // Print Regions
  fout << "    <DataArray type=\"UInt16\" Name=\"Regions\" format=\"binary\">" << endl;
  
  // point pointer to common mamory area buffer of void type;
  unsigned short* var_reg=static_cast <unsigned short*> (buffer_void);
  icount=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned ii=0; ii<_msh[ig]->GetElementNumber(); ii++) {
      if ( ig==_gridn-1u || 0==_msh[ig]->el->GetRefinedElementIndex(ii)) {
	unsigned iel_Metis = _msh[ig]->GetMetisDof(ii,3);
	var_reg[icount]= _msh[ig]->el->GetElementGroup(ii);
	icount++;
      }
    }
  }
  
  //print regions dimension
  cch = b64::b64_encode(&dim_array_reg[0], sizeof(dim_array_reg), NULL, 0);  
  b64::b64_encode(&dim_array_reg[0], sizeof(dim_array_reg), &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  
  //print regions array
  cch = b64::b64_encode(&var_reg[0], dim_array_reg[0] , NULL, 0);  
  b64::b64_encode(&var_reg[0], dim_array_reg[0], &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  fout << endl;
  fout << "    </DataArray>" << endl;
  //-----------------------------------------------------------------------------------------------------   
  // Print Metis Partitioning
  fout << "    <DataArray type=\"UInt16\" Name=\"Domain_partition\" format=\"binary\">" << endl;
  
  // point pointer to common mamory area buffer of void type;
  unsigned short* var_proc=static_cast <unsigned short*> (buffer_void);
  icount=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned ii=0; ii<_msh[ig]->GetElementNumber(); ii++) {
      if ( ig==_gridn-1u || 0==_msh[ig]->el->GetRefinedElementIndex(ii)) {
	unsigned iel_Metis = _msh[ig]->GetMetisDof(ii,3);
	var_proc[icount]=(unsigned short)(_msh[ig]->epart[ii]);
	icount++;
      }
    }
  }
  
  //print regions dimension
  cch = b64::b64_encode(&dim_array_reg[0], sizeof(dim_array_reg), NULL, 0);  
  b64::b64_encode(&dim_array_reg[0], sizeof(dim_array_reg), &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  
  //print regions array
  cch = b64::b64_encode(&var_proc[0], dim_array_reg[0] , NULL, 0);  
  b64::b64_encode(&var_proc[0], dim_array_reg[0], &enc[0], cch);
  pt_char=&enc[0];
  for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
  fout << endl;
  fout << "    </DataArray>" << endl;
  

  //Print Solution (on element) ***************************************************************
  for (unsigned i=0; i<(1-test_all)*vars.size()+test_all*SolName.size(); i++) {
    unsigned indx=(test_all==0)?GetIndex(vars[i].c_str()):i;
    if (3 <= SolType[indx]) {
      fout << "    <DataArray type=\"Float32\" Name=\"" << SolName[indx] <<"\" format=\"binary\">" << endl;
      // point pointer to common memory area buffer of void type;
      float *var_el = static_cast< float*> (buffer_void);
      icount=0;
      for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
	vector<double> sol_local;
	_solution[ig]->_Sol[indx]->localize_to_one(sol_local,0);
	for (unsigned ii=0; ii<_msh[ig]->GetElementNumber(); ii++) {
	  if (ig==_gridn-1u || 0==_msh[ig]->el->GetRefinedElementIndex(ii)) {
	    unsigned iel_Metis = _msh[ig]->GetMetisDof(ii,SolType[indx]);
	    var_el[icount]=sol_local[iel_Metis];
	    icount++;
	  }
	}
      }
      
      if(_iproc==0) {
        //print solution on element dimension
        cch = b64::b64_encode(&dim_array_elvar[0], sizeof(dim_array_elvar), NULL, 0);  
        b64::b64_encode(&dim_array_elvar[0], sizeof(dim_array_elvar), &enc[0], cch);
        pt_char=&enc[0];
        for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
    
        //print solution on element array
        cch = b64::b64_encode(&var_el[0], dim_array_elvar[0] , NULL, 0);  
        b64::b64_encode(&var_el[0], dim_array_elvar[0], &enc[0], cch);
        pt_char=&enc[0];
        for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
        fout << endl;
        fout << "    </DataArray>" << endl;
      }
      //----------------------------------------------------------------------------------------------------
    }
  }
  fout << "   </CellData>" << endl;
  //   //------------------------------------------------------------------------------------------------
  // 
  //   //------------------------------------------------------------------------------------------------
  // / Print Solution (on nodes) ********************************************************************
  fout<< " <PointData Scalars=\"scalars\"> " << endl;
  //Loop on variables
   
  // point pointer to common mamory area buffer of void type;
  float* var_nd = static_cast<float*>(buffer_void);
  for (unsigned i=0; i<(1-test_all)*vars.size()+test_all*SolName.size(); i++) {
    unsigned indx=(test_all==0)?GetIndex(vars[i].c_str()):i;
    if (SolType[indx]<3) {
      fout << " <DataArray type=\"Float32\" Name=\"" << SolName[indx] <<"\" format=\"binary\">" << endl;
      //print solutions on nodes dimension
      cch = b64::b64_encode(&dim_array_ndvar[0], sizeof(dim_array_ndvar), NULL, 0);  
      b64::b64_encode(&dim_array_ndvar[0], sizeof(dim_array_ndvar), &enc[0], cch);
      pt_char=&enc[0];
      for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
      
      unsigned offset_nvt=0;
      for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
	mysol[ig]->matrix_mult(*_solution[ig]->_Sol[indx],*ProlQitoQj_[index_nd][SolType[indx]][ig]);
	vector<double> sol_local;
	mysol[ig]->localize_to_one(sol_local,0);
	unsigned nvt_ig=_msh[ig]->MetisOffset[index_nd][_nprocs];
	for (unsigned ii=0; ii<nvt_ig; ii++) {
	  var_nd[ii+offset_nvt] = sol_local[ii];
	}
	offset_nvt+=nvt_ig;
      }
      
      if(_iproc==0) {
        cch = b64::b64_encode(&var_nd[0], dim_array_ndvar [0], NULL, 0);  
        b64::b64_encode(&var_nd[0], dim_array_ndvar [0], &enc[0], cch);
        pt_char=&enc[0];
        for( unsigned i =0; i<cch;i++,pt_char++) fout << *pt_char; 
        fout << endl;

        fout << "    </DataArray>" << endl;
      }
    } //endif
  } // end for sol
  fout << "   </PointData>" << endl;
  
  //------------------------------------------------------------------------------------------------
  

  fout << "  </Piece>" << endl;
  fout << " </UnstructuredGrid>" << endl;
  fout << "</VTKFile>" << endl;
  fout.close();

  
  
  //-----------------------------------------------------------------------------------------------------
  //free memory
  for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    delete mysol[ig];
  }
  delete [] filename;
  delete [] var_nd;  
  
  //--------------------------------------------------------------------------------------------------------
}


// *************************************************************************
// This function prints the solution in Xdmf - hdf5 format
// *************************************************************************

void  NonLinearMultiLevelProblem::printsol_xdmf_hdf5(const char type[], std::vector<std::string>& vars) const {
  
  
  if(_iproc!=0) return;
  
  bool test_all=!(vars[0].compare("All"));
    
  unsigned index=0;
  unsigned index_nd=0;
  if(!strcmp(type,"linear")) {    //linear
    index=0;
    index_nd=0;
  }
  else if(!strcmp(type,"quadratic")) {  //quadratic
    index=1;
    index_nd=1;
  }
  else if(!strcmp(type,"biquadratic")) { //biquadratic
    index=3;
    index_nd=2;
  }

  const string type_el[4][6] = {{"Hexahedron","Tetrahedron","Wedge","Quadrilateral","Triangle","Edge"},
                                {"Hexahedron_20","Tetrahedron_10","Not_implemented","Quadrilateral_8","Triangle_6","Edge_3"},
			        {"Not_implemented","Not_implemented","Not_implemented","Not_implemented",
				 "Not_implemented","Not_implemented"},
                                {"Not_implemented","Not_implemented","Not_implemented","Quadrilateral_9",
				 "Not_implemented","Not_implemented"}};
			 
  
  //I assume that the mesh is not mixed
  string type_elem;
  type_elem = type_el[index][_msh[_gridn-1u]->el->GetElementType(0)];
  
  if (type_elem.compare("Not_implemented") == 0) exit(1);
  
  unsigned nvt=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    unsigned nvt_ig=_msh[ig]->GetDofNumber(index_nd);
    nvt+=nvt_ig;
  } 
  
  // Printing connectivity
  unsigned nel=0;
  for(unsigned ig=0;ig<_gridn-1u;ig++) {
    nel+=( _msh[ig]->GetElementNumber() - _msh[ig]->el->GetRefinedElementNumber());
  }
  nel+=_msh[_gridn-1u]->GetElementNumber();
  
  unsigned icount;
  unsigned el_dof_number  = _msh[_gridn-1u]->el->GetElementDofNumber(0,index);
  int *var_int             = new int [nel*el_dof_number];
  float *var_el_f         = new float [nel];
  float *var_nd_f         = new float [nvt];

  char *filename= new char[60];
  std::ofstream fout;
  
  //--------------------------------------------------------------------------------------------------
  // Print The Xdmf wrapper
  sprintf(filename,"./output/mesh.level%d.%d.%s.xmf",_gridn,_time_step,type);
  //std::ofstream fout;
  fout.open(filename);
  if (!fout) {
    cout << "Output mesh file "<<filename<<" cannot be opened.\n";
    exit(0);
  }
  
  // Print The HDF5 file
  sprintf(filename,"./mesh.level%d.%d.%s.h5",_gridn,_time_step,type);
  // haed ************************************************
  fout<<"<?xml version=\"1.0\" ?>" << endl;
  fout<<"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd []\">"<<endl;
  fout<<"<Xdmf>"<<endl;
  fout<<"<Domain>"<<endl;
  fout<<"<Grid Name=\"Mesh\">"<<endl;
  fout<<"<Time Value =\""<< _time_step<< "\" />"<<endl;
  fout<<"<Topology TopologyType=\""<< type_elem <<"\" NumberOfElements=\""<< nel <<"\">"<<endl;
  //Connectivity
  fout<<"<DataStructure DataType=\"Int\" Dimensions=\""<< nel*el_dof_number <<"\"" << "  Format=\"HDF\">" << endl;
  fout << filename << ":CONNECTIVITY" << endl;
  fout <<"</DataStructure>" << endl;
  fout << "</Topology>" << endl;
  fout << "<Geometry Type=\"X_Y_Z\">" << endl;
  //Node_X
  fout<<"<DataStructure DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nvt << "  1\"" << "  Format=\"HDF\">" << endl;
  fout << filename << ":NODES_X1" << endl;
  fout <<"</DataStructure>" << endl;
  //Node_Y
  fout<<"<DataStructure DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nvt << "  1\"" << "  Format=\"HDF\">" << endl;
  fout << filename << ":NODES_X2" << endl;
  fout <<"</DataStructure>" << endl;
  //Node_Z
  fout<<"<DataStructure DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nvt << "  1\"" << "  Format=\"HDF\">" << endl;
  fout << filename << ":NODES_X3" << endl;
  fout <<"</DataStructure>" << endl;
  fout <<"</Geometry>" << endl;
  //Regions
  fout << "<Attribute Name=\""<< "Regions"<<"\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl;
  fout << "<DataItem DataType=\"Int\" Dimensions=\""<< nel << "\"" << "  Format=\"HDF\">" << endl;
  fout << filename << ":REGIONS" << endl;
  fout << "</DataItem>" << endl;
  fout << "</Attribute>" << endl;
  // Solution Variables
  for (unsigned i=0; i<vars.size(); i++) {
    unsigned indx=GetIndex(vars[i].c_str());  
    //Printing biquadratic solution on the nodes
    if(SolType[indx]<3) {  
      fout << "<Attribute Name=\""<< SolName[indx]<<"\" AttributeType=\"Scalar\" Center=\"Node\">" << endl;
      fout << "<DataItem DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nvt << "  1\"" << "  Format=\"HDF\">" << endl;
      fout << filename << ":" << SolName[indx] << endl;
      fout << "</DataItem>" << endl;
      fout << "</Attribute>" << endl;
    }
    else if (SolType[indx]>=3) {  //Printing picewise constant solution on the element
      fout << "<Attribute Name=\""<< SolName[indx]<<"\" AttributeType=\"Scalar\" Center=\"Cell\">" << endl;
      fout << "<DataItem DataType=\"Float\" Precision=\"4\" Dimensions=\""<< nel << "\"  Format=\"HDF\">" << endl;
      fout << filename << ":" << SolName[indx] << endl;
      fout << "</DataItem>" << endl;
      fout << "</Attribute>" << endl;
    }
  }

  fout <<"</Grid>" << endl;
  fout <<"</Domain>" << endl;
  fout <<"</Xdmf>" << endl;
  fout.close();
  //----------------------------------------------------------------------------------------------------------
  
  //----------------------------------------------------------------------------------------------------------
  hid_t file_id;
  sprintf(filename,"./output/mesh.level%d.%d.%s.h5",_gridn,_time_step,type);
  file_id = H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  hsize_t dimsf[2];
  herr_t status;
  hid_t dataspace;
  hid_t dataset;
  
  //-----------------------------------------------------------------------------------------------------------
  // Printing nodes coordinates 
  
  PetscScalar *MYSOL[1]; //TODO
  unsigned varind[3];
  varind[0]=GetIndex("X");
  varind[1]=GetIndex("Y");
  varind[2]=GetIndex("Z");

  
  for (int i=0; i<3; i++) {
    unsigned offset_nvt=0;
    for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
      NumericVector* mysol;
      mysol = NumericVector::build().release();
      mysol->init(_msh[ig]->GetDofNumber(index_nd),_msh[ig]->GetDofNumber(index_nd),true,SERIAL);
      mysol->matrix_mult(*_solution[ig]->_Sol[varind[i]],*ProlQitoQj_[index_nd][SolType[varind[i]]][ig]);
      unsigned nvt_ig=_msh[ig]->GetDofNumber(index_nd);
      for (unsigned ii=0; ii<nvt_ig; ii++) var_nd_f[ii+offset_nvt] = (*mysol)(ii);
      if (_moving_mesh) {
	unsigned varind_DXDYDZ=GetIndex(_moving_vars[i].c_str());
	mysol->matrix_mult(*_solution[ig]->_Sol[varind_DXDYDZ],*ProlQitoQj_[index_nd][SolType[varind_DXDYDZ]][ig]);
	for (unsigned ii=0; ii<nvt_ig; ii++) var_nd_f[ii+offset_nvt] += (*mysol)(ii);
      }
      offset_nvt+=nvt_ig;
      delete mysol;
    }
    
    dimsf[0] = nvt ;  dimsf[1] = 1;
    std::ostringstream Name; Name << "/NODES_X" << i+1;
    dataspace = H5Screate_simple(2,dimsf, NULL);
    dataset   = H5Dcreate(file_id,Name.str().c_str(),H5T_NATIVE_FLOAT,
			  dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_nd_f[0]);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    
  }
  
  //-------------------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------------------------------
  //connectivity
  icount = 0;
  unsigned offset_conn=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned iel=0; iel<_msh[ig]->GetElementNumber(); iel++) {
      if (_msh[ig]->el->GetRefinedElementIndex(iel)==0 || ig==_gridn-1u) {
        for (unsigned j=0; j<_msh[ig]->el->GetElementDofNumber(iel,index); j++) {
	  unsigned jnode=_msh[ig]->el->GetElementVertexIndex(iel,j)-1u;
	  unsigned jnode_Metis = _msh[ig]->GetMetisDof(jnode,index_nd);
	  var_int[icount] = offset_conn + jnode_Metis;
	  icount++;
	}
      }
    }
    offset_conn += _msh[ig]->GetDofNumber(index_nd);
  }
  
  dimsf[0] = nel*el_dof_number ;  dimsf[1] = 1;
  dataspace = H5Screate_simple(2,dimsf, NULL);
  dataset   = H5Dcreate(file_id,"/CONNECTIVITY",H5T_NATIVE_INT,
			dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status   = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_int[0]);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  //------------------------------------------------------------------------------------------------------
  
  
  //-------------------------------------------------------------------------------------------------------
  // print regions
  icount=0;
  for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
    for (unsigned ii=0; ii<_msh[ig]->GetElementNumber(); ii++) {
      if (ig==_gridn-1u || 0==_msh[ig]->el->GetRefinedElementIndex(ii)) {
	unsigned iel_Metis = _msh[ig]->GetMetisDof(ii,3);
	var_int[icount] = _msh[ig]->el->GetElementGroup(iel_Metis);
	icount++;
      }
    }
  } 
   
  dimsf[0] = nel;  dimsf[1] = 1;
  dataspace = H5Screate_simple(2,dimsf, NULL);
  dataset   = H5Dcreate(file_id,"/REGIONS",H5T_NATIVE_INT,
			dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status   = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_int[0]);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  
  //-------------------------------------------------------------------------------------------------------
  // printing element variables
  for (unsigned i=0; i<(1-test_all)*vars.size()+test_all*SolName.size(); i++) {
    unsigned indx=(test_all==0)?GetIndex(vars[i].c_str()):i;
    if (SolType[indx]>=3) {
      icount=0;
      for (unsigned ig=_gridr-1u; ig<_gridn; ig++) {
	for (unsigned ii=0; ii<_msh[ig]->GetElementNumber(); ii++) {
	  if (ig==_gridn-1u || 0==_msh[ig]->el->GetRefinedElementIndex(ii)) {
	    unsigned iel_Metis = _msh[ig]->GetMetisDof(ii,SolType[indx]);
	    var_el_f[icount]=(*_solution[ig]->_Sol[indx])(iel_Metis);
	    icount++;
	  }
	}
      } 
     
      dimsf[0] = nel;  dimsf[1] = 1;
      dataspace = H5Screate_simple(2,dimsf, NULL);
      dataset   = H5Dcreate(file_id,SolName[indx],H5T_NATIVE_FLOAT,
			    dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status   = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_el_f[0]);
      H5Sclose(dataspace);
      H5Dclose(dataset);
     
    }
  }
  
  //-------------------------------------------------------------------------------------------------------
  // printing nodes variables
  for (unsigned i=0; i<(1-test_all)*vars.size()+test_all*SolName.size(); i++) {
    unsigned indx=(test_all==0)?GetIndex(vars[i].c_str()):i;
    if (SolType[indx] < 3) {
      unsigned offset_nvt=0;
      for(unsigned ig=_gridr-1u; ig<_gridn; ig++) {
        NumericVector* mysol;
	mysol = NumericVector::build().release();
        mysol->init(_msh[ig]->GetDofNumber(index_nd),_msh[ig]->GetDofNumber(index_nd),true,SERIAL);
	mysol->matrix_mult(*_solution[ig]->_Sol[indx],*ProlQitoQj_[index_nd][SolType[indx]][ig]);
	unsigned nvt_ig=_msh[ig]->GetDofNumber(index_nd);
	for (unsigned ii=0; ii<nvt_ig; ii++) var_nd_f[ii+offset_nvt] = (*mysol)(ii);
	offset_nvt+=nvt_ig;
	delete mysol;
      }
     
      dimsf[0] = nvt;  dimsf[1] = 1;
      dataspace = H5Screate_simple(2,dimsf, NULL);
      dataset   = H5Dcreate(file_id,SolName[indx],H5T_NATIVE_FLOAT,
			    dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status   = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT,&var_nd_f[0]);
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
  }
  //-------------------------------------------------------------------------------------------------------
    
  // Close the file -------------
  H5Fclose(file_id);
 
  //free memory
  delete [] filename;
  delete [] var_int;
  delete [] var_el_f;
  delete [] var_nd_f;
  
}

// hid_t NonLinearMultiLevelProblem::print_Dhdf5(hid_t file,const std::string & name, hsize_t dimsf[],double data[]) {
//   hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
//   hid_t dataset   = H5Dcreate(file,name.c_str(),H5T_NATIVE_DOUBLE,
// 			      dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//   hid_t  status   = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
//   H5Sclose(dataspace);
//   H5Dclose(dataset);
//   return status;
// }




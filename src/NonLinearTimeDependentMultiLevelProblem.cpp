#include "NonLinearTimeDependentMultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "ElemType.hpp"
#include "LinearSolver.hpp"
#include "PetscVector.hpp"
#include <iostream>
#include <fstream>
using std::cout;
using std::endl;

// *******************************************************************
// ****************** MGTimeLoop CLASS ********************************
// *******************************************************************

///Constructor
NonLinearTimeDependentMultiLevelProblem::NonLinearTimeDependentMultiLevelProblem(const unsigned short &igridn,const unsigned short &igridr, const char mesh_file[],
                                     const char GaussOrder[], const double Lref, bool (* SetRefinementFlag)(const double &x, const double &y, const double &z, 
                     const int &ElemGroupNumber,const int &level)) :
  NonLinearMultiLevelProblem(igridn, igridr, mesh_file, GaussOrder, Lref,SetRefinementFlag) {

  _print_step    = 100000;
  _save_step     = 100000;
  _ats_flag      = 0;
  _test_time     = 1;
  _ntime_steps   =1;
  _dt            =1.;
}


///Destructor
NonLinearTimeDependentMultiLevelProblem::~NonLinearTimeDependentMultiLevelProblem() {}

//------------------------------------------------------------------------------------------------------

void NonLinearTimeDependentMultiLevelProblem::SetTimeStep(const double dt) {
  _dt = dt;
}

//------------------------------------------------------------------------------------------------------

void NonLinearTimeDependentMultiLevelProblem::SetPrintTimeStep(const unsigned print_step)  {
  _print_step = print_step;
}

//------------------------------------------------------------------------------------------------------

void NonLinearTimeDependentMultiLevelProblem::SetSaveTimeStep(const unsigned save_step)  {
  _save_step = save_step;
}

//------------------------------------------------------------------------------------------------------

double NonLinearTimeDependentMultiLevelProblem::GetTimeStep() {
  return _dt;
}

//------------------------------------------------------------------------------------------------------

double NonLinearTimeDependentMultiLevelProblem::GetTime() {
  return _time;
}

//------------------------------------------------------------------------------------------------------

void NonLinearTimeDependentMultiLevelProblem::SetInitTimeStep(const unsigned time_step0) {
  _time_step0 = time_step0;
}

//------------------------------------------------------------------------------------------------------

void NonLinearTimeDependentMultiLevelProblem::SetNumTimeSteps(const double ntimesteps) {
  _ntime_steps = ntimesteps;
}

//------------------------------------------------------------------------------------------------------

unsigned NonLinearTimeDependentMultiLevelProblem::GetInitTimeStep() const {
  return _time_step0;
}

//------------------------------------------------------------------------------------------------------

unsigned NonLinearTimeDependentMultiLevelProblem::GetNumTimeSteps() const {
  return _ntime_steps;
}

//------------------------------------------------------------------------------------------------------

unsigned NonLinearTimeDependentMultiLevelProblem::GetPrintTimeStep() const {
  return _print_step;
}

//------------------------------------------------------------------------------------------------------

unsigned NonLinearTimeDependentMultiLevelProblem::GetSaveTimeStep() const {
  return _save_step;
}

//------------------------------------------------------------------------------------------------------

void NonLinearTimeDependentMultiLevelProblem::AttachSetTimeStepFunction (double (* set_time_step_function)(const double time)) {
  _set_time_step_function = set_time_step_function;
  _ats_flag = 1;
}

//------------------------------------------------------------------------------------------------------

int NonLinearTimeDependentMultiLevelProblem::SaveData() const {
  char *filename = new char[80];
  PetscVector* petsc_vec_sol;

  PetscErrorCode ierr;
  PetscViewer bin_viewer;
  for (unsigned ig=0; ig<gridn; ig++) {
    for (unsigned i=0; i<SolType.size(); i++) {
      sprintf(filename,"./save/save.time_%d.level_%d.%s.bin",_time_step,ig,SolName[i]);
      ierr = PetscViewerBinaryOpen(MPI_COMM_WORLD,filename,FILE_MODE_WRITE,&bin_viewer);
      CHKERRQ(ierr);
      //petsc_vec_sol  = static_cast<PetscVector*>(Lin_Solver_[ig]->Sol_[i]);
      petsc_vec_sol  = static_cast<PetscVector*>(_solution[ig]->_Sol[i]);
      ierr = VecView(petsc_vec_sol->vec(),bin_viewer);
      CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&bin_viewer);
    }
  }
  delete [] filename;
  return ierr;
}

//------------------------------------------------------------------------------------------------------

int NonLinearTimeDependentMultiLevelProblem::InitializeFromRestart(unsigned restart_time_step) {

  // Set the restart time step
  SetInitTimeStep(restart_time_step+1);
  
  //Set the restart time
  if(_ats_flag==0) {
    _time = restart_time_step*_dt;
  } 
  else {
   _time = 0;
   for(unsigned i=0; i<restart_time_step; i++) {
     _dt = _set_time_step_function(_time);
     _time += _dt;
   }
  }
   
  //Set the restart time step
  _time_step = restart_time_step+1;
     
  char *filename = new char[80];
  PetscVector* petsc_vec_sol;

  PetscErrorCode ierr;
  PetscViewer bin_viewer;
  for(unsigned ig=0;ig<gridn;ig++){
    for(unsigned i=0;i<SolType.size();i++){
      sprintf(filename,"./save/save.time_%d.level_%d.%s.bin",restart_time_step,ig,SolName[i]);
      ierr = PetscViewerBinaryOpen(MPI_COMM_WORLD,filename,FILE_MODE_READ,&bin_viewer); CHKERRQ(ierr);
      //petsc_vec_sol  = static_cast<PetscVector*>(Lin_Solver_[ig]->Sol_[i]);
      petsc_vec_sol  = static_cast<PetscVector*>(_solution[ig]->_Sol[i]);
      ierr = VecLoad(petsc_vec_sol->vec(),bin_viewer); CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&bin_viewer);
    }
    _solution[ig]->UpdateSolution();
  }
  
  delete [] filename;
  return ierr;
}

//--------------------------------------------------------------------------------------------------------

void NonLinearTimeDependentMultiLevelProblem::_UpdateSolution() {
 
    for (int ig=0; ig<gridn; ig++) {
      _solution[ig]->UpdateSolution();
    }
  
  return;
}


//---------------------------------------------------------------------------------------------------------

void NonLinearTimeDependentMultiLevelProblem::_NewmarkAccUpdate() {

  //Static conversion from the Base class (MultiGrid) to the Derived class (NonLinearTimeDependentMultiLevelProblem)
//   NonLinearTimeDependentMultiLevelProblem& mg = static_cast<NonLinearTimeDependentMultiLevelProblem&>(mg2);

//   const double dt    = GetTimeStep();
//   cout << dt << endl;
  const double gamma = 0.5;
  const double a5    = -1.*(1. - gamma)/gamma;
  const double a1    = 1./(gamma*_dt);
  const double a2    = -1./(gamma*_dt);
  
  unsigned dim = _msh[0]->GetDimension();

  unsigned axyz[3];
  unsigned vxyz[3];
  const char accname[3][3] = {"AX","AY","AZ"};
  const char velname[3][2] = {"U","V","W"};
  
  for(unsigned i=0; i<dim; i++) {
    axyz[i] = GetIndex(&accname[i][0]);
    vxyz[i] = GetIndex(&velname[i][0]);
  }
  
  for (int ig=0;ig<gridn;ig++) {
//     unsigned ax=GetIndex("AX");
//     unsigned ay=GetIndex("AY");
//     unsigned az=GetIndex("AZ");
//     unsigned vx=GetIndex("U");
//     unsigned vy=GetIndex("V");
//     unsigned vz=GetIndex("W");
    for(unsigned i=0; i<dim; i++) {
//       Lin_Solver_[ig]->Sol_[axyz[i]]->scale(a5);
//       Lin_Solver_[ig]->Sol_[axyz[i]]->add(a1,*(Lin_Solver_[ig]->Sol_[vxyz[i]]));
//       Lin_Solver_[ig]->Sol_[axyz[i]]->add(a2,*(Lin_Solver_[ig]->Sol_old_[vxyz[i]]));
      
      _solution[ig]->_Sol[axyz[i]]->scale(a5);
      _solution[ig]->_Sol[axyz[i]]->add(a1,*(_solution[ig]->_Sol[vxyz[i]]));
      _solution[ig]->_Sol[axyz[i]]->add(a2,*(_solution[ig]->_SolOld[vxyz[i]]));
      
      
    }


//     Lin_Solver_[ig]->Sol_[ax]->scale(a5);
//     Lin_Solver_[ig]->Sol_[ax]->add(a1,*(Lin_Solver_[ig]->Sol_[vx]));
//     Lin_Solver_[ig]->Sol_[ax]->add(a2,*(Lin_Solver_[ig]->Sol_old_[vx]));
//     
//     Lin_Solver_[ig]->Sol_[ay]->scale(a5);
//     Lin_Solver_[ig]->Sol_[ay]->add(a1,*(Lin_Solver_[ig]->Sol_[vy]));
//     Lin_Solver_[ig]->Sol_[ay]->add(a2,*(Lin_Solver_[ig]->Sol_old_[vy]));
//     
//     Lin_Solver_[ig]->Sol_[az]->scale(a5);
//     Lin_Solver_[ig]->Sol_[az]->add(a1,*(Lin_Solver_[ig]->Sol_[vz]));
//     Lin_Solver_[ig]->Sol_[az]->add(a2,*(Lin_Solver_[ig]->Sol_old_[vz]));
  }

}

//---------------------------------------------------------------------------------------------------------

///Time Dependent Full Multigrid
int NonLinearTimeDependentMultiLevelProblem::FullMultiGrid(unsigned const &ncycle,  unsigned const &npre, 
							   unsigned const &npost, const char mg_type[]) {
  
  
  clock_t start_time, end_time, start_cycle_time, end_cycle_time, start_mg_time, end_mg_time;
  int conv;
  
  bool flagmc = 1;
  if (!strcmp(mg_type,"V-Cycle")) {
    flagmc = 0;
  }
  
  cout << "flagmc: " << flagmc << endl;
  if (_ats_flag) {
    //update time step
    _dt = _set_time_step_function(_time);
  }

  //update time
  _time += _dt;

  //update time step
  _time_step++;
    
  //update boundary condition
  UpdateBdc();

  start_mg_time = clock();

  cout << endl;
  cout << " ************** time step " << _time_step << " **************** " << endl;
  cout << " ************** time      " << _time << " **************** " << endl;
  
  for (unsigned igridn=1*flagmc + (!flagmc)*gridn; igridn<=gridn; igridn++) {
    cout << endl;
    cout << "    ************* Level Max: " << igridn << " *************";
    cout << endl;
    for (unsigned icycle=0; icycle<ncycle; icycle++) {
      start_cycle_time = clock();
      cout << endl;
      cout << "    ************** Cycle: " << icycle + 1 << " **************** " << endl;
      cout << endl;

      start_time=clock();
      Lin_Solver_[igridn-1u]->SetResZero();
      Lin_Solver_[igridn-1u]->SetEpsZero();
      end_time=clock();
      cout<<"Grid: "<<igridn-1<<"      INITIALIZATION TIME:      "
	  <<static_cast<double>((end_time-start_time))/CLOCKS_PER_SEC<<endl;

      start_time=clock();
      //assemble residual and matrix on the finer grid at igridn level
      _assemble_function(*this,igridn-1u,igridn-1u);

      for (unsigned ig=igridn-1u; ig>0; ig--) {
	start_time=clock();

        Lin_Solver_[ig-1u]->SetResZero();
        Lin_Solver_[ig-1u]->SetEpsZero();
	
        if (ig>=gridr) {
          //assemble residual only on the part of the coarse grid that is not refined
          //Domain Decomposition matrix restriction =========================
          _assemble_function(*this,ig-1,igridn-1u);
	   
	  if (!Lin_Solver_[ig-1]->CC_flag) {
            MatPtAP(Lin_Solver_[ig]->KK, Lin_Solver_[ig]->PP,  MAT_INITIAL_MATRIX ,1.0,&Lin_Solver_[ig-1]->CC);
	    Lin_Solver_[ig-1]->CC_flag=1;
          } else MatPtAP(Lin_Solver_[ig]->KK, Lin_Solver_[ig]->PP,  MAT_REUSE_MATRIX ,1.0,&Lin_Solver_[ig-1]->CC);
	  MatAXPY(Lin_Solver_[ig-1u]->KK,1,Lin_Solver_[ig-1u]->CC, SUBSET_NONZERO_PATTERN);
        } 
	else {
          if (icycle==0 && ( flagmc*(ig==igridn-1u) || !flagmc )) {
            MatDestroy(&Lin_Solver_[ig-1]->KK);
            MatPtAP(Lin_Solver_[ig]->KK, Lin_Solver_[ig]->PP, MAT_INITIAL_MATRIX ,1.0,&Lin_Solver_[ig-1]->KK);
          }
          //Projection of the Matrix on the lower level
          else
            MatPtAP(Lin_Solver_[ig]->KK, Lin_Solver_[ig]->PP, MAT_REUSE_MATRIX ,1.0,&Lin_Solver_[ig-1]->KK);
        }
        
        end_time=clock();
        cout<<"Grid: "<<ig-1<<"      ASSEMBLY + RESIDUAL TIME: "
	    <<static_cast<double>((end_time-start_time))/CLOCKS_PER_SEC<<endl;
      }
      
      
      end_time=clock();
      cout<<"Grid: "<<igridn-1<<"      ASSEMBLY + RESIDUAL TIME: "
	  <<static_cast<double>((end_time-start_time))/CLOCKS_PER_SEC<<endl;

      /// Presmoothing  
      for (unsigned ig=igridn-1u; ig>0; ig--) {
	for (unsigned k=0; k<npre; k++) {
          if (ig==ig) {
	    if(_VankaIsSet) {
             if (_VankaSchur) Lin_Solver_[ig]->Vanka_Smoother(MGIndex,VankaIndex,_NSchurVar,_Schur);
             else Lin_Solver_[ig]->Vanka_Smoother(MGIndex,VankaIndex);
	    }  else {
              Lin_Solver_[ig]->solve();
	    }
          }
        }

	start_time=clock();

        //standard Multigrid matrix restriction =========================
        Restrictor(ig); //restriction of the residual
        end_time=clock();
        cout<<"Grid: "<<ig<<"-->"<<ig-1<<"  RESTRICTION TIME:         "
	    <<static_cast<double>((end_time-start_time))/CLOCKS_PER_SEC<<endl;

      }
      
        /// Coarse direct solver
        if(_VankaIsSet) {
	  if (_VankaSchur)Lin_Solver_[0]->Vanka_Smoother(MGIndex,VankaIndex,_NSchurVar,_Schur);
          else Lin_Solver_[0]->Vanka_Smoother(MGIndex,VankaIndex);
	} else {
  	  Lin_Solver_[0]->solve();
	}
      

	for (unsigned ig=1; ig<igridn; ig++) {
	  start_time=clock();
	  Prolungator(ig);
	  Lin_Solver_[ig]->UpdateResidual();
	  Lin_Solver_[ig]->SumEpsCToEps();
	  end_time=clock();
	  cout<<"Grid: "<<ig-1<<"-->"<<ig<<"  PROLUNGATION TIME:        "
	      <<static_cast<double>((end_time-start_time))/CLOCKS_PER_SEC<<endl;

	  /// PostSmoothing    
	  for (unsigned k=0; k<npost; k++) {
	    if (ig==ig) {
	      if(_VankaIsSet) {
	        if (_VankaSchur) Lin_Solver_[ig]->Vanka_Smoother(MGIndex,VankaIndex,_NSchurVar,_Schur);
	        else Lin_Solver_[ig]->Vanka_Smoother(MGIndex,VankaIndex);
	      } else {
  	        Lin_Solver_[ig]->solve();
	      }
	    }
	  }
	}

	cout << endl;
	start_time=clock();
	for (unsigned ig=0; ig<igridn; ig++) {
	  _solution[ig]->SumEpsToSol(MGIndex, Lin_Solver_[ig]->EPS, Lin_Solver_[ig]->RES, Lin_Solver_[ig]->KKoffset );
	}

	conv  = GetConvergence(igridn-1);
	if (conv ==true) icycle = ncycle + 1;

	end_time=clock();
	cout << endl;
	cout<<"COMPUTATION RESIDUAL:                  "<<static_cast<double>((end_time-start_time))/CLOCKS_PER_SEC<<endl;

	end_cycle_time=clock();
	cout<<"CYCLE TIME:                            "<<static_cast<double>((end_cycle_time-start_cycle_time))/CLOCKS_PER_SEC<<endl;
    }
    //only for the Full Multicycle
    if (igridn<gridn) {
      ProlungatorSol(igridn);
    }
  }

  for (int ig=gridr-1; ig<gridn-1; ig++) {
    if(Lin_Solver_[ig]->CC_flag){
      MatDestroy(&(Lin_Solver_[ig]->CC));
      Lin_Solver_[ig]->CC_flag=0;
    }
  }
    
  end_mg_time = clock();
  cout<<"STEADYSOLVER TIME:                            "<<static_cast<double>((end_mg_time-start_mg_time))/CLOCKS_PER_SEC<<endl;

  return 1;
  
  
  
  
  
//   clock_t start_time, end_time, start_cycle_time, end_cycle_time,start_time_step, end_time_step;
//   clock_t start_time1, end_time1;
//   int conv;
//   start_time_step=clock();
//   
//     bool flagmc = 1;
//     if (_time_step>1 && !strcmp(mg_type,"V-Cycle")) {
//       flagmc = 0;
//     }
// 
//     if (_ats_flag) {
//       //update time step
//       _dt = _set_time_step_function(_time);
//     }
// 
//     //update time
//     _time += _dt;
// 
//     //update time step
//     _time_step++;
//     
//     //update boundary condition
//     UpdateBdc();
// 
//     cout << endl;
//     cout << " ************** time step " << _time_step << " **************** " << endl;
//     cout << " ************** time      " << _time << " **************** " << endl;
//    
//     // Here I do a V-Multigrid and not a FullMultiGrid
//     for (unsigned igridn=1*flagmc + (!flagmc)*gridn; igridn<=gridn; igridn++) {
//     cout << endl;
//     cout << "    ************* Level Max: " << igridn << " *************";
//     for (unsigned icycle=0; icycle<ncycle; icycle++) {
//       start_cycle_time = clock();
//       cout << endl;
//       cout << "    ************** Cycle: " << icycle + 1 << " **************** " << endl;
//       cout << endl;
// 
//       start_time=clock();
//       Lin_Solver_[igridn-1u]->SetResZero(MGIndex);
//       Lin_Solver_[igridn-1u]->SetEpsZero(MGIndex);
//       end_time=clock();
//       cout<<"Grid: "<<igridn-1<<"      INITIALIZATION TIME:      "
//       <<static_cast<double>((end_time-start_time))/CLOCKS_PER_SEC<<endl;
// 
//       start_time=clock();
//       //assemble residual and matrix on the finer grid at igridn level
//       _assemble_function(*this,igridn-1u,igridn-1u);
// 
// 
//       end_time=clock();
//       cout<<"Grid: "<<igridn-1<<"      ASSEMBLY + RESIDUAL TIME: "
//       <<static_cast<double>((end_time-start_time))/CLOCKS_PER_SEC<<endl;
// 
//       //  PetscViewer viewer;
//       //       PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);
//       //       for(int ig=igridn-1;ig>0;ig--){
//       // 	MatView(D[ig]->KK,viewer);
//       // 	double ff;
//       //         std::cin>>ff;
//       //       }
// 
// 
//       for (unsigned ig=igridn-1u; ig>0; ig--) {
// 
//         for (unsigned k=0; k<npre; k++) {
//           if (ig==ig) {
//             if (_VankaSchur) Lin_Solver_[ig]->Vanka_Smoother(MGIndex,VankaIndex,_NSchurVar,_Schur);
//             else Lin_Solver_[ig]->Vanka_Smoother(MGIndex,VankaIndex);
//           }
//         }
// 
//         start_time=clock();
//         Lin_Solver_[ig-1u]->SetResZero(MGIndex);
//         Lin_Solver_[ig-1u]->SetEpsZero(MGIndex);
// 
//         //standard Multigrid matrix restriction =========================
//         Restrictor(ig); //restriction of the residual
//         end_time=clock();
//         cout<<"Grid: "<<ig<<"-->"<<ig-1<<"  RESTRICTION TIME:         "
//         <<static_cast<double>((end_time-start_time))/CLOCKS_PER_SEC<<endl;
// 
// 
// 
//         start_time=clock();
// 
//         if (ig>=gridr) {
//           //assemble residual only on the part of the coarse grid that is not refined
//           //Domain Decomposition matrix restriction =========================
//           _assemble_function(*this,ig-1,igridn-1u);
// 	   
// 	  if (!Lin_Solver_[ig-1]->CC_flag) {
//             MatPtAP(Lin_Solver_[ig]->KK, Lin_Solver_[ig]->PP,  MAT_INITIAL_MATRIX ,1.0,&Lin_Solver_[ig-1]->CC);
// 	    Lin_Solver_[ig-1]->CC_flag=1;
//           } else MatPtAP(Lin_Solver_[ig]->KK, Lin_Solver_[ig]->PP,  MAT_REUSE_MATRIX ,1.0,&Lin_Solver_[ig-1]->CC);
// 	  MatAXPY(Lin_Solver_[ig-1u]->KK,1,Lin_Solver_[ig-1u]->CC, SUBSET_NONZERO_PATTERN);
//         } else {
//           if (icycle==0 && ig==igridn-1u) {
//             MatDestroy(&Lin_Solver_[ig-1]->KK);
//             MatPtAP(Lin_Solver_[ig]->KK, Lin_Solver_[ig]->PP, MAT_INITIAL_MATRIX ,1.0,&Lin_Solver_[ig-1]->KK);
//           }
//           //Projection of the Matrix on the lower level
//           else
//             MatPtAP(Lin_Solver_[ig]->KK, Lin_Solver_[ig]->PP, MAT_REUSE_MATRIX ,1.0,&Lin_Solver_[ig-1]->KK);
//         }
// 
//         end_time=clock();
//         cout<<"Grid: "<<ig-1<<"      ASSEMBLY + RESIDUAL TIME: "
//         <<static_cast<double>((end_time-start_time))/CLOCKS_PER_SEC<<endl;
// 
//       }
// 
// 
//       //coarse smoother
//       if (_VankaSchur)Lin_Solver_[0]->Vanka_Smoother(MGIndex,VankaIndex,_NSchurVar,_Schur);
//       else Lin_Solver_[0]->Vanka_Smoother(MGIndex,VankaIndex);
// 
// 
//       for (unsigned ig=1; ig<igridn; ig++) {
//         start_time=clock();
//         Prolungator(ig);
//         Lin_Solver_[ig]->UpdateResidual();
//         Lin_Solver_[ig]->SumEpsCToEps(MGIndex);
//         end_time=clock();
//         cout<<"Grid: "<<ig-1<<"-->"<<ig<<"  PROLUNGATION TIME:        "
//         <<static_cast<double>((end_time-start_time))/CLOCKS_PER_SEC<<endl;
// 
//         for (unsigned k=0; k<npost; k++) {
//           if (ig==ig) {
//             if (_VankaSchur) Lin_Solver_[ig]->Vanka_Smoother(MGIndex,VankaIndex,_NSchurVar,_Schur);
//             else Lin_Solver_[ig]->Vanka_Smoother(MGIndex,VankaIndex);
//           }
//         }
//       }
// 
//       cout << endl;
//       start_time=clock();
//       for (unsigned ig=0; ig<igridn; ig++) {
//         Lin_Solver_[ig]->SumEpsToSol(MGIndex);
//       }
// 
//       
//  //     conv = Lin_Solver_[igridn-1]->comp_toll(MGIndex);
//       conv  = GetConvergence(igridn-1);
//       if (conv ==true) icycle = ncycle + 1;
// 
// 
//       end_time=clock();
//       cout << endl;
//       cout<<"COMPUTATION RESIDUAL:                  "<<static_cast<double>((end_time-start_time))/CLOCKS_PER_SEC<<endl;
// 
//       //if(igridn==1) icycle+=100;
//       end_cycle_time=clock();
//       cout<<"CYCLE TIME:                            "<<static_cast<double>((end_cycle_time-start_cycle_time))/CLOCKS_PER_SEC<<endl;
//     }
//     //only for the Full Multicycle
//     if (igridn<gridn) {
//       ProlungatorSol(igridn);
//     }
//   }
// 
//   for (int ig=gridr-1; ig<gridn-1; ig++) {
//     if(Lin_Solver_[ig]->CC_flag){
//       MatDestroy(&(Lin_Solver_[ig]->CC));
//       Lin_Solver_[ig]->CC_flag=0;
//     }
//   }
//     
//   end_time_step=clock();
//   cout<<"TIME_STEP TIME:                        "<<static_cast<double>((end_time_step-start_time_step))/CLOCKS_PER_SEC<<endl;
// 
//   return 1;
}


//--------------------------------------------------------------------------------------------------------

void NonLinearTimeDependentMultiLevelProblem::UpdateBdc() {

  const short unsigned NV1[6][2]= {{9,9},{6,6},{9,6},{3,3},{3,3},{1,1}};
  unsigned dim = _msh[0]->GetDimension() - 1u;
  unsigned indX=GetIndex("X");
  unsigned indY=GetIndex("Y");
  unsigned indZ=GetIndex("Z");
  double xx;
  double yy;
  double zz;

  for (int k=0; k<SolName.size(); k++) {
    if (!strcmp(BdcType[k],"Time_dependent")) {
      unsigned i_start = k;
      unsigned i_end;
      i_end=i_start+1u;
      
      
      // 2 Default Neumann
      // 1 DD Dirichlet
      // 0 Dirichlet
      for (unsigned igridn=0; igridn<gridn; igridn++) {
	if(_solution[igridn]->_ResEpsBdcFlag[k]){
	  for (unsigned j=_msh[igridn]->MetisOffset[SolType[k]][_iproc]; j<_msh[igridn]->MetisOffset[SolType[k]][_iproc+1]; j++) {
	    _solution[igridn]->_Bdc[k]->set(j,2.);
	  }
	  if (SolType[k]<3) {  
	    for(int isdom=_iproc; isdom<_iproc+1; isdom++) {   // 1 DD Dirichlet
	      for (int iel=_msh[igridn]->IS_Mts2Gmt_elem_offset[isdom]; 
		   iel < _msh[igridn]->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
		unsigned kel_gmt = _msh[igridn]->IS_Mts2Gmt_elem[iel];
		for (unsigned jface=0; jface<_msh[igridn]->el->GetElementFaceNumber(kel_gmt); jface++) {
		  if (_msh[igridn]->el->GetFaceElementIndex(kel_gmt,jface)==0) { //Domain Decomposition Dirichlet
		    short unsigned ielt=_msh[igridn]->el->GetElementType(kel_gmt);
		    unsigned nv1=(!TestIfPressure[k])?
				  NV1[ielt][jface<_msh[igridn]->el->GetElementFaceNumber(kel_gmt,0)]:
				  _msh[igridn]->el->GetElementDofNumber(iel,_msh[igridn]->GetEndIndex(SolType[k]));
		    for (unsigned iv=0; iv<nv1; iv++) {
		      unsigned inode=(!TestIfPressure[k])? 
				      _msh[igridn]->el->GetFaceVertexIndex(kel_gmt,jface,iv)-1u:
				      _msh[igridn]->el->GetElementVertexIndex(kel_gmt,iv)-1u;
		      unsigned inode_Metis=_msh[igridn]->GetMetisDof(inode,SolType[k]);
		      _solution[igridn]->_Bdc[k]->set(inode_Metis,1.);
		    }
		  }
		}
	      }
	    }
	    for(int isdom=_iproc; isdom<_iproc+1; isdom++) {  // 0 Dirichlet
	      for (int iel=_msh[igridn]->IS_Mts2Gmt_elem_offset[isdom]; 
		  iel < _msh[igridn]->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
		unsigned kel_gmt = _msh[igridn]->IS_Mts2Gmt_elem[iel];
		for (unsigned jface=0; jface<_msh[igridn]->el->GetElementFaceNumber(kel_gmt); jface++) {
		  if (_msh[igridn]->el->GetFaceElementIndex(kel_gmt,jface)<0) { //Dirichlet
		    short unsigned ielt=_msh[igridn]->el->GetElementType(kel_gmt);
		    unsigned nv1=NV1[ielt][jface<_msh[igridn]->el->GetElementFaceNumber(kel_gmt,0)];
		    for (unsigned iv=0; iv<nv1; iv++) {
		      unsigned inode=_msh[igridn]->el->GetFaceVertexIndex(kel_gmt,jface,iv)-1u;
		      unsigned inode_coord_Metis=_msh[igridn]->GetMetisDof(inode,2);
		      double value;
		      xx=(*_solution[igridn]->_Sol[indX])(inode_coord_Metis);  
		      yy=(*_solution[igridn]->_Sol[indY])(inode_coord_Metis);
		      zz=(*_solution[igridn]->_Sol[indZ])(inode_coord_Metis);
		      bool test=_SetBoundaryConditionFunction(xx,yy,zz,SolName[k],value,-(_msh[igridn]->el->GetFaceElementIndex(kel_gmt,jface)+1),_time);
		      if (test) {
			unsigned inode_Metis=_msh[igridn]->GetMetisDof(inode,SolType[k]);
			_solution[igridn]->_Bdc[k]->set(inode_Metis,0.);
			_solution[igridn]->_Sol[k]->set(inode_Metis,value);
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  else if(TestIfPressure[k]){
	    for(int isdom=_iproc; isdom<_iproc+1; isdom++) {   // 1 DD Dirichlet for pressure variable only
	      unsigned nel=_msh[igridn]->GetElementNumber();
	      for (int iel=_msh[igridn]->IS_Mts2Gmt_elem_offset[isdom]; 
		  iel < _msh[igridn]->IS_Mts2Gmt_elem_offset[isdom+1]; iel++) {
		unsigned kel_gmt = _msh[igridn]->IS_Mts2Gmt_elem[iel];
		for (unsigned jface=0; jface<_msh[igridn]->el->GetElementFaceNumber(kel_gmt); jface++) {
		  if (_msh[igridn]->el->GetFaceElementIndex(kel_gmt,jface)==0) { //Domain Decomposition Dirichlet
		    short unsigned ielt=_msh[igridn]->el->GetElementType(kel_gmt);
		    unsigned nv1=_msh[igridn]->el->GetElementDofNumber(kel_gmt,_msh[igridn]->GetEndIndex(SolType[k]));
		    for (unsigned iv=0; iv<nv1; iv++) {
		      unsigned inode=(kel_gmt+iv*nel);
		      unsigned inode_Metis=_msh[igridn]->GetMetisDof(inode,SolType[k]);
		      _solution[igridn]->_Bdc[k]->set(inode_Metis,1.);
		    }
		  }
		}
	      }
	    }
	  }	
	  _solution[igridn]->_Sol[k]->close();
	  _solution[igridn]->_Bdc[k]->close();
	}
      }
    }
  }
}

//------------------------------------------------------------------------------------------------------------
void NonLinearTimeDependentMultiLevelProblem::printsol_xdmf_archive(const char type[]) const {
  
  char *filename= new char[60];
  // Print The Xdmf transient wrapper
  sprintf(filename,"./output/mesh.level%d.%s.xmf",gridn,type);
  std::ofstream ftr_out;
  ftr_out.open(filename);
  if (!ftr_out) {
    cout << "Transient Output mesh file "<<filename<<" cannot be opened.\n";
    exit(0);
  }
  
    ftr_out << "<?xml version=\"1.0\" ?> \n";
    ftr_out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd []\">"<<endl;
    ftr_out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> " << endl;
    ftr_out << "<Domain> " << endl;
    ftr_out << "<Grid Name=\"Mesh\" GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
    // time loop for grid sequence
    for ( unsigned time_step = _time_step0; time_step < _time_step0 + _ntime_steps; time_step++) {
     if ( !(time_step%_print_step) ) {
       sprintf(filename,"./mesh.level%d.%d.%s.xmf",gridn,time_step,type);
       ftr_out << "<xi:include href=\"" << filename << "\" xpointer=\"xpointer(//Xdmf/Domain/Grid["<< 1 <<"])\">\n";
        ftr_out << "<xi:fallback/>\n";
        ftr_out << "</xi:include>\n";
       }
      }
     ftr_out << "</Grid> \n";
     ftr_out << "</Domain> \n";
     ftr_out << "</Xdmf> \n";
     ftr_out.close();
   ftr_out.close();  
  //----------------------------------------------------------------------------------------------------------
  delete [] filename;
 
}



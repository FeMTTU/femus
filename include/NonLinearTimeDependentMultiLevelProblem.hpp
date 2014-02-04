#ifndef __NonLinearTimeDependentMultiLevelProblem_hpp__
#define __NonLinearTimeDependentMultiLevelProblem_hpp__

#include "NonLinearMultiLevelProblem.hpp"

//******************** Time loop class **************************
//***************************************************************
class NonLinearTimeDependentMultiLevelProblem : public NonLinearMultiLevelProblem {

private:
  unsigned _save_step;
  unsigned _print_step;
  bool _ats_flag;      //adaptive time step flag

  //pointer function to the set time step function
  double (* _set_time_step_function)(const double time);

  void SetInitTimeStep(const unsigned time_step0);
  
  void UpdateBdc();

public:

  ///Constructor
  NonLinearTimeDependentMultiLevelProblem(const unsigned short &igridn, const unsigned short &igridr,
                    const char mesh_file[], const char GaussOrder[], const double Lref=1., bool (* SetRefinementFlag)(const double &x, const double &y, const double &z, 
                     const int &ElemGroupNumber,const int &level)=NULL);

  ///Destructor
  ~NonLinearTimeDependentMultiLevelProblem();
 
  void SetTimeStep(const double dt);

  void SetPrintTimeStep(const unsigned print_step);

  void SetSaveTimeStep(const unsigned save_step);
  
  unsigned GetInitTimeStep() const;

  unsigned GetNumTimeSteps() const;
  
  unsigned GetPrintTimeStep() const;

  unsigned GetSaveTimeStep() const;

  double GetTimeStep();

  double GetTime();
  
  int SaveData() const;
  
  void _NewmarkAccUpdate();

  void _UpdateSolution();
    
  void AttachSetTimeStepFunction (double (* set_time_step_function)(const double time));

  int InitializeFromRestart(unsigned restart_time_step);

  void SetNumTimeSteps(const double ntimesteps);

  int FullMultiGrid(const char pdename[], unsigned const &ncycle,  unsigned const &npost, unsigned const &npre, 
		    const char mg_type[]="F-Cycle");
  
  void printsol_xdmf_archive(const char type[]) const;

};

#endif
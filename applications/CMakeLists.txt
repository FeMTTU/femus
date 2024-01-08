#############################################################################################
### Build Applications
#############################################################################################


### Build options, Applications - BEGIN

# --
OPTION(BUILD_TUTORIAL "Build the tutorial problems" ON)


# --
OPTION(BUILD_POISSON "Build the Poisson problem" ON)


# Ask the user which executable have to be built
OPTION(BUILD_NAVIERSTOKESAPPS "Build the NavierStokes-line applications" ON)

# Ask the user which executable have to be built
OPTION(BUILD_FSIAPPS "Build the FSI applications" ON)

# --
OPTION(BUILD_OPTCONTROLAPPS "Build the optimal control applications" ON)


# --
OPTION(BUILD_AMR "Build the AMR apps" ON)

OPTION(BUILD_MGAMR "Build the MGAMR" ON)

OPTION(BUILD_REFINEMENT "Build the Refinement apps" ON)




# --
OPTION(BUILD_PCFIELDSPLIT "Build the PCFIELDSPLIT apps" ON)


# --
OPTION(BUILD_RK "Build the Runge-Kutta apps" ON)


# --
OPTION(BUILD_NONLOCAL "Build the NONLOCAL apps" ON)

# --
OPTION(BUILD_UQ "Build the UQ apps" ON)



# --
OPTION(BUILD_CONFORMAL "Build the Conformal mesh apps" ON)

OPTION(BUILD_WILLMORE "Build the Willmore problem" ON)





# --
OPTION(BUILD_ISM "Build the Intersection Marker apps" ON)

OPTION(BUILD_MPM_FEM "Build the MPM FEM apps" ON)

OPTION(BUILD_MPM_FSI "Build the MPM FSI apps" ON)




OPTION(BUILD_SW "Build the Shallow Water apps" ON)

OPTION(BUILD_TUMOR "Build the Tumor apps" ON)


OPTION(BUILD_NITSCHE "Build the Nitsche apps" ON)




OPTION(BUILD_ZZZ_SANDBOX "Build the sandbox apps for playing around" ON)

### Build options, Applications - END



### Add subdirectories - BEGIN

IF(BUILD_TUTORIAL)
  ADD_SUBDIRECTORY(000_tutorial/)
ENDIF(BUILD_TUTORIAL)


IF(BUILD_POISSON)
  ADD_SUBDIRECTORY(001_Poisson/)
ENDIF(BUILD_POISSON)


IF(BUILD_NAVIERSTOKESAPPS)
  set(NS_FOLDER 003_NavierStokes/)
  # Steady NavierStokes
  ADD_SUBDIRECTORY(${NS_FOLDER}/SteadyNavierStokesParallel)

  # Time-dependent BE NavierStokes
  ADD_SUBDIRECTORY(${NS_FOLDER}/BackwardEulerNavierStokesParallel/)

  # Time-dependent Generalized-alpha NavierStokes
  ADD_SUBDIRECTORY(${NS_FOLDER}/TaylorGreenVortexFlowGeneralisedalpha/)

  # Steady Stokes
  ADD_SUBDIRECTORY(${NS_FOLDER}/SteadyStokes/)

  ADD_SUBDIRECTORY(${NS_FOLDER}/NavierStokes_AMR/)
  
  set(NSTEMP_FOLDER 004_Boussinesq/)
  # Steady Boussinesq
  ADD_SUBDIRECTORY(${NSTEMP_FOLDER}/NavierStokesTemp_AMR/)
  unset(NS_FOLDER)
ENDIF(BUILD_NAVIERSTOKESAPPS)


IF(BUILD_FSIAPPS)
  set(FSI_FOLDER 005_FSI/)
  ADD_SUBDIRECTORY(${FSI_FOLDER})
  unset(FSI_FOLDER)
  
ENDIF(BUILD_FSIAPPS)


IF(BUILD_OPTCONTROLAPPS)
  ADD_SUBDIRECTORY(010_OptimalControl/)
ENDIF(BUILD_OPTCONTROLAPPS)



IF(BUILD_RK)
  ADD_SUBDIRECTORY(050_RK/)
ENDIF(BUILD_RK)



IF(BUILD_SW)
  ADD_SUBDIRECTORY(090_SW/)
ENDIF(BUILD_SW)

IF(BUILD_TUMOR)
  ADD_SUBDIRECTORY(091_Tumor/)
ENDIF(BUILD_TUMOR)




IF(BUILD_WILLMORE)
  ADD_SUBDIRECTORY(Willmore/)
ENDIF(BUILD_WILLMORE)


IF(BUILD_CONFORMAL)
  ADD_SUBDIRECTORY(Conformal/)
ENDIF(BUILD_CONFORMAL)




IF(BUILD_AMR)
  ADD_SUBDIRECTORY(AMR)
ENDIF(BUILD_AMR)


IF(BUILD_MGAMR)
  ADD_SUBDIRECTORY(MGAMR)
ENDIF(BUILD_MGAMR)


IF(BUILD_REFINEMENT)
  ADD_SUBDIRECTORY(Refinement/)
ENDIF(BUILD_REFINEMENT)





IF(BUILD_PCFIELDSPLIT)
  ADD_SUBDIRECTORY(PCFieldSplit)
ENDIF(BUILD_PCFIELDSPLIT)




IF(BUILD_ISM)
  ADD_SUBDIRECTORY(ISM/)
ENDIF(BUILD_ISM)

IF(BUILD_MPM_FEM)
  ADD_SUBDIRECTORY(MPM_FEM/)
ENDIF(BUILD_MPM_FEM)

IF(BUILD_MPM_FSI)
  ADD_SUBDIRECTORY(MPM_FSI/)
ENDIF(BUILD_MPM_FSI)


IF(BUILD_UQ)
  ADD_SUBDIRECTORY(UQ/)
ENDIF(BUILD_UQ)

IF(BUILD_NONLOCAL)
  ADD_SUBDIRECTORY(NonLocal/)
ENDIF(BUILD_NONLOCAL)



IF(BUILD_NITSCHE)
  ADD_SUBDIRECTORY(Nitsche/)
ENDIF(BUILD_NITSCHE)



IF(BUILD_ZZZ_SANDBOX)
  ADD_SUBDIRECTORY(zzz_sandbox/)
ENDIF(BUILD_ZZZ_SANDBOX)


### Add subdirectories - END

#############################################################################################
### Build Applications
#############################################################################################

IF(BUILD_NAVIERSTOKESAPPS)
  set(NS_FOLDER NavierStokes/)
  # Steady NavierStokes
  ADD_SUBDIRECTORY(${NS_FOLDER}/SteadyNavierStokesParallel)

  # Time-dependent BE NavierStokes
  ADD_SUBDIRECTORY(${NS_FOLDER}/BackwardEulerNavierStokesParallel/)

  # Time-dependent Generalized-alpha NavierStokes
  ADD_SUBDIRECTORY(${NS_FOLDER}/TaylorGreenVortexFlowGeneralisedalpha/)

  # Steady Stokes
  ADD_SUBDIRECTORY(${NS_FOLDER}/SteadyStokes/)

  # Steady Navier Stokes steady
  ADD_SUBDIRECTORY(${NS_FOLDER}/NavierStokesTemp_AMR/)
  unset(NS_FOLDER)
ENDIF(BUILD_NAVIERSTOKESAPPS)


IF(BUILD_FSIAPPS)
  set(FSI_FOLDER FSI/)
  ADD_SUBDIRECTORY(${FSI_FOLDER})
  unset(FSI_FOLDER)
  
ENDIF(BUILD_FSIAPPS)


IF(BUILD_OPTCONTROLAPPS)
  ADD_SUBDIRECTORY(OptimalControl/)
ENDIF(BUILD_OPTCONTROLAPPS)

IF(BUILD_POISSON)
  ADD_SUBDIRECTORY(Poisson/)
ENDIF(BUILD_POISSON)

IF(BUILD_POISSON_AMR)
  ADD_SUBDIRECTORY(Poisson_AMR/)
ENDIF(BUILD_POISSON_AMR)

IF(BUILD_TUTORIAL)
  ADD_SUBDIRECTORY(tutorial/)
ENDIF(BUILD_TUTORIAL)

IF(BUILD_WILLMORE)
  ADD_SUBDIRECTORY(Willmore/)
ENDIF(BUILD_WILLMORE)

IF(BUILD_MGAMR)
  ADD_SUBDIRECTORY(MGAMR)
ENDIF(BUILD_MGAMR)

IF(BUILD_PCFIELDSPLIT)
  ADD_SUBDIRECTORY(PCFieldSplit)
ENDIF(BUILD_PCFIELDSPLIT)

IF(BUILD_AMR)
  ADD_SUBDIRECTORY(AMR)
ENDIF(BUILD_AMR)

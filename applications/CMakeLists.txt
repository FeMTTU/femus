#############################################################################################
### Build Applications 
#############################################################################################

IF(BUILD_NAVIERSTOKESAPPS)
  
  # Steady NavierStokes
  ADD_SUBDIRECTORY(NavierStokesApps/SteadyNavierStokesParallel)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/SteadyNavierStokesParallel/output/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/SteadyNavierStokesParallel/input/)
  FILE(COPY           NavierStokesApps/SteadyNavierStokesParallel/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/SteadyNavierStokesParallel/input/)

  # Time-dependent BE NavierStokes
  ADD_SUBDIRECTORY(NavierStokesApps/BackwardEulerNavierStokesParallel/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/BackwardEulerNavierStokesParallel/output/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/BackwardEulerNavierStokesParallel/input/)
  FILE(COPY           NavierStokesApps/BackwardEulerNavierStokesParallel/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/BackwardEulerNavierStokesParallel/input/)
 
  # Time-dependent Generalized-alpha NavierStokes
  ADD_SUBDIRECTORY(NavierStokesApps/TaylorGreenVortexFlowGeneralisedalpha/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/TaylorGreenVortexFlowGeneralisedalpha/output/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/TaylorGreenVortexFlowGeneralisedalpha/input/)
  FILE(COPY           NavierStokesApps/TaylorGreenVortexFlowGeneralisedalpha/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/TaylorGreenVortexFlowGeneralisedalpha/input/)
 
  # Steady Stokes
  ADD_SUBDIRECTORY(NavierStokesApps/SteadyStokes/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/SteadyStokes/output/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/SteadyStokes/input/)
  FILE(COPY           NavierStokesApps/SteadyStokes/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/SteadyStokes/input/)

  # Steady Navier Stokes steady
  ADD_SUBDIRECTORY(NavierStokesApps/NavierStokesTemp_AMR/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/NavierStokesTemp_AMR/output/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/NavierStokesTemp_AMR/input/)
  FILE(COPY           NavierStokesApps/NavierStokesTemp_AMR/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/NavierStokesApps/NavierStokesTemp_AMR/input/)
 
ENDIF(BUILD_NAVIERSTOKESAPPS)


IF(BUILD_FSIAPPS)

# FSI SteadyState Benchmark Comsol
    ADD_SUBDIRECTORY(FSIApps/FSI_CMS_NWM/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_CMS_NWM/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_CMS_NWM/input/)
    FILE(COPY           FSIApps/FSI_CMS_NWM/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_CMS_NWM/input/)

# FSI SteadyState Benchmark Turek1
    ADD_SUBDIRECTORY(FSIApps/FSI_TUREK_1/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_TUREK_1/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_TUREK_1/input/)
    FILE(COPY           FSIApps/FSI_TUREK_1/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_TUREK_1/input/)

# FSI SteadyState Benchmark Turek3
    ADD_SUBDIRECTORY(FSIApps/FSI_TUREK_3/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_TUREK_3/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_TUREK_3/input/)
    FILE(COPY           FSIApps/FSI_TUREK_3/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/FSIApps/FSI_TUREK_3/input/)

# FSI SteadyState in porous media
    ADD_SUBDIRECTORY(FSIApps/fsi_porous/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/FSIApps/fsi_porous/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/FSIApps/fsi_porous/input/)
    FILE(COPY           FSIApps/fsi_porous/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/FSIApps/fsi_porous/input/)

ENDIF(BUILD_FSIAPPS)


IF(BUILD_GENCASE)
    ADD_SUBDIRECTORY(gencase/)
ENDIF(BUILD_GENCASE)


IF(BUILD_MHDOPT)
       ADD_SUBDIRECTORY(mhdopt/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/mhdopt/input/)
    FILE(COPY mhdopt/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/mhdopt/input/)
ENDIF(BUILD_MHDOPT)


IF(BUILD_TEMPER)
       ADD_SUBDIRECTORY(temper/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/temper/input/)
    FILE(COPY temper/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/temper/input/)
ENDIF(BUILD_TEMPER)


IF(BUILD_POISSON)
       ADD_SUBDIRECTORY(Poisson/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/Poisson/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/Poisson/input/)
    FILE(COPY           Poisson/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/Poisson/input/)
ENDIF(BUILD_POISSON)


IF(BUILD_TEMPERATURE)
       ADD_SUBDIRECTORY(temperature/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/temperature/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/temperature/input/)
    FILE(COPY           temperature/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/temperature/input/)
ENDIF(BUILD_TEMPERATURE)


IF(BUILD_POISSON_AMR)
       ADD_SUBDIRECTORY(Poisson_AMR/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/Poisson_AMR/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/Poisson_AMR/input/)
    FILE(COPY           Poisson_AMR/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/Poisson_AMR/input/)
ENDIF(BUILD_POISSON_AMR)
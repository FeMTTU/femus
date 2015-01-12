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




set(OPTCONTROL_FOLDER OptControlApps/)

IF(BUILD_OPTCONTROLAPPS)

# MHD Optimal Control
    set(OPT_APPLICATION mhdopt)
    ADD_SUBDIRECTORY(${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/input/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/output/)
    FILE(COPY ${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/input/)

# Temperature Optimal Control
    set(OPT_APPLICATION tempopt)
    ADD_SUBDIRECTORY(${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/input/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/output/)
    FILE(COPY ${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/input/)

# Temperature test
    set(OPT_APPLICATION fe_test)
    ADD_SUBDIRECTORY(${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/input/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/output/)
    FILE(COPY ${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/${OPTCONTROL_FOLDER}/${OPT_APPLICATION}/input/)

    unset(OPT_APPLICATION)
    
ENDIF(BUILD_OPTCONTROLAPPS)





IF(BUILD_POISSON)
       ADD_SUBDIRECTORY(Poisson/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/Poisson/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/Poisson/input/)
    FILE(COPY           Poisson/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/Poisson/input/)
ENDIF(BUILD_POISSON)


IF(BUILD_POISSON_AMR)
       ADD_SUBDIRECTORY(Poisson_AMR/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/Poisson_AMR/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/Poisson_AMR/input/)
    FILE(COPY           Poisson_AMR/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/Poisson_AMR/input/)
ENDIF(BUILD_POISSON_AMR)
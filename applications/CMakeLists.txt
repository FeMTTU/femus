#############################################################################################
### Build Applications 
#############################################################################################

IF(BUILD_NAVIERSTOKESAPPS)
  
  # Steady NavierStokes
  ADD_SUBDIRECTORY(NavierStokesApps/SteadyNavierStokesParallel)

  # Time-dependent BE NavierStokes
  ADD_SUBDIRECTORY(NavierStokesApps/BackwardEulerNavierStokesParallel/)
 
  # Time-dependent Generalized-alpha NavierStokes
  ADD_SUBDIRECTORY(NavierStokesApps/TaylorGreenVortexFlowGeneralisedalpha/)
 
  # Steady Stokes
  ADD_SUBDIRECTORY(NavierStokesApps/SteadyStokes/)

  # Steady Navier Stokes steady
  ADD_SUBDIRECTORY(NavierStokesApps/NavierStokesTemp_AMR/)
 
ENDIF(BUILD_NAVIERSTOKESAPPS)


IF(BUILD_FSIAPPS)

# FSI SteadyState Benchmark Comsol
    ADD_SUBDIRECTORY(FSIApps/FSITimeDependentBenchmarks/)

# FSI SteadyState Benchmark Turek1
    ADD_SUBDIRECTORY(FSIApps/FSISteadyStateBenchmarks/)

# FSI SteadyState Benchmark Turek3
    ADD_SUBDIRECTORY(FSIApps/FSI_TUREK_3/)

# FSI SteadyState in porous media
    ADD_SUBDIRECTORY(FSIApps/fsi_porous/)

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
ENDIF(BUILD_POISSON)


IF(BUILD_POISSON_AMR)
       ADD_SUBDIRECTORY(Poisson_AMR/)
ENDIF(BUILD_POISSON_AMR)


IF(BUILD_NUM_PDE_CLASS)

#   ADD_SUBDIRECTORY(NumPdeSpring2015/Eugenio/)
  
  ADD_SUBDIRECTORY(NumPdeSpring2015/Giacomo/)
  
  ADD_SUBDIRECTORY(NumPdeSpring2015/Sara/)
  
  ADD_SUBDIRECTORY(NumPdeSpring2015/Thanuja/)
  
  ADD_SUBDIRECTORY(NumPdeSpring2015/Thanuka/)
  
  ADD_SUBDIRECTORY(NumPdeSpring2015/Bhagya/)
  
  ADD_SUBDIRECTORY(NumPdeSpring2015/Guoyi/)
        
  set(TUTORIAL_FOLDER NumPdeSpring2015/tutorial)  
  ADD_SUBDIRECTORY(${TUTORIAL_FOLDER})
  
#   set(EUGENIO_FOLDER NumPdeSpring2015/Eugenio)  
#   ADD_SUBDIRECTORY(${EUGENIO_FOLDER})  

ENDIF(BUILD_NUM_PDE_CLASS)


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
ENDIF(BUILD_NAVIERSTOKESAPPS)


IF(BUILD_FSIAPPS)

set(FSI_FOLDER FSI/)

# FSI Time Dependent Benchmarks
    ADD_SUBDIRECTORY(${FSI_FOLDER}/FSITimeDependentBenchmarks/)

# FSI SteadyState Benchmarks
    ADD_SUBDIRECTORY(${FSI_FOLDER}/FSISteadyStateBenchmarks/)

# FSI SteadyState Benchmark Turek3
    ADD_SUBDIRECTORY(${FSI_FOLDER}/FSI_TUREK_3/)

# FSI SteadyState in porous media
    ADD_SUBDIRECTORY(${FSI_FOLDER}/fsi_porous/)

ENDIF(BUILD_FSIAPPS)


IF(BUILD_OPTCONTROLAPPS)

set(OPTCONTROL_FOLDER OptimalControl/)

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

# FE test
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
        
  set(TUTORIAL_FOLDER NumPdeSpring2015/tutorial)  
  ADD_SUBDIRECTORY(${TUTORIAL_FOLDER})
  
  set(EUGENIO_FOLDER NumPdeSpring2015/Eugenio)  
  ADD_SUBDIRECTORY(${EUGENIO_FOLDER})  

  set(BHAGYA_FOLDER NumPdeSpring2015/Bhagya)  
  ADD_SUBDIRECTORY(${BHAGYA_FOLDER})
  
  set(SARA_FOLDER NumPdeSpring2015/Sara)  
  ADD_SUBDIRECTORY(${SARA_FOLDER})
  
  set(GIACOMO_FOLDER NumPdeSpring2015/Giacomo)  
  ADD_SUBDIRECTORY(${GIACOMO_FOLDER})
  
  set(THANUJA_FOLDER NumPdeSpring2015/Thanuja)  
  ADD_SUBDIRECTORY(${THANUJA_FOLDER})
  
  set(THANUKA_FOLDER NumPdeSpring2015/Thanuka)  
  ADD_SUBDIRECTORY(${THANUKA_FOLDER})
  
  set(GUOYI_FOLDER NumPdeSpring2015/Guoyi)  
  ADD_SUBDIRECTORY(${GUOYI_FOLDER})
  
ENDIF(BUILD_NUM_PDE_CLASS)


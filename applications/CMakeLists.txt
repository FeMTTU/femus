#############################################################################################
### Build Applications 
#############################################################################################

IF(BUILD_NAVIERSTOKESAPPS)
  
set(NS_FOLDER NavierStokes/)

  # Steady NavierStokes
  ADD_SUBDIRECTORY(${NS_FOLDER}/SteadyNavierStokesParallel)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/SteadyNavierStokesParallel/output/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/SteadyNavierStokesParallel/input/)
  FILE(COPY           ${NS_FOLDER}/SteadyNavierStokesParallel/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/SteadyNavierStokesParallel/input/)

  # Time-dependent BE NavierStokes
  ADD_SUBDIRECTORY(${NS_FOLDER}/BackwardEulerNavierStokesParallel/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/BackwardEulerNavierStokesParallel/output/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/BackwardEulerNavierStokesParallel/input/)
  FILE(COPY           ${NS_FOLDER}/BackwardEulerNavierStokesParallel/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/BackwardEulerNavierStokesParallel/input/)
 
  # Time-dependent Generalized-alpha NavierStokes
  ADD_SUBDIRECTORY(${NS_FOLDER}/TaylorGreenVortexFlowGeneralisedalpha/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/TaylorGreenVortexFlowGeneralisedalpha/output/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/TaylorGreenVortexFlowGeneralisedalpha/input/)
  FILE(COPY           ${NS_FOLDER}/TaylorGreenVortexFlowGeneralisedalpha/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/TaylorGreenVortexFlowGeneralisedalpha/input/)
 
  # Steady Stokes
  ADD_SUBDIRECTORY(${NS_FOLDER}/SteadyStokes/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/SteadyStokes/output/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/SteadyStokes/input/)
  FILE(COPY           ${NS_FOLDER}/SteadyStokes/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/SteadyStokes/input/)

  # Steady Navier Stokes steady
  ADD_SUBDIRECTORY(${NS_FOLDER}/NavierStokesTemp_AMR/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/NavierStokesTemp_AMR/output/)
  FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/NavierStokesTemp_AMR/input/)
  FILE(COPY           ${NS_FOLDER}/NavierStokesTemp_AMR/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/${NS_FOLDER}/NavierStokesTemp_AMR/input/)
 
ENDIF(BUILD_NAVIERSTOKESAPPS)


IF(BUILD_FSIAPPS)

set(FSI_FOLDER FSI/)

# FSI Time Dependent Benchmarks
    ADD_SUBDIRECTORY(${FSI_FOLDER}/FSITimeDependentBenchmarks/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${FSI_FOLDER}/FSITimeDependentBenchmarks/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${FSI_FOLDER}/FSITimeDependentBenchmarks/input/)
    FILE(COPY           ${FSI_FOLDER}/FSITimeDependentBenchmarks/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/${FSI_FOLDER}/FSITimeDependentBenchmarks/input/)

# FSI SteadyState Benchmarks
    ADD_SUBDIRECTORY(${FSI_FOLDER}/FSISteadyStateBenchmarks/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${FSI_FOLDER}/FSISteadyStateBenchmarks/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${FSI_FOLDER}/FSISteadyStateBenchmarks/input/)
    FILE(COPY           ${FSI_FOLDER}/FSISteadyStateBenchmarks/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/${FSI_FOLDER}/FSISteadyStateBenchmarks/input/)

# FSI SteadyState Benchmark Turek3
    ADD_SUBDIRECTORY(${FSI_FOLDER}/FSI_TUREK_3/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${FSI_FOLDER}/FSI_TUREK_3/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${FSI_FOLDER}/FSI_TUREK_3/input/)
    FILE(COPY           ${FSI_FOLDER}/FSI_TUREK_3/input/ DESTINATION 
  ${PROJECT_BINARY_DIR}/applications/${FSI_FOLDER}/FSI_TUREK_3/input/)

# FSI SteadyState in porous media
    ADD_SUBDIRECTORY(${FSI_FOLDER}/fsi_porous/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${FSI_FOLDER}/fsi_porous/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/${FSI_FOLDER}/fsi_porous/input/)
    FILE(COPY           ${FSI_FOLDER}/fsi_porous/input/ DESTINATION  ${PROJECT_BINARY_DIR}/applications/${FSI_FOLDER}/fsi_porous/input/)

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


IF(BUILD_NUM_PDE_CLASS)
    ADD_SUBDIRECTORY(NumPdeSpring2015/Eugenio/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Eugenio/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Eugenio/input/)
    FILE(COPY           NumPdeSpring2015/Eugenio/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Eugenio/input/)
  
  ADD_SUBDIRECTORY(NumPdeSpring2015/Giacomo/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Giacomo/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Giacomo/input/)
    FILE(COPY           NumPdeSpring2015/Giacomo/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Giacomo/input/)
  
  ADD_SUBDIRECTORY(NumPdeSpring2015/Sara/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Sara/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Sara/input/)
    FILE(COPY           NumPdeSpring2015/Sara/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Sara/input/)
  
  ADD_SUBDIRECTORY(NumPdeSpring2015/Thanuja/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Thanuja/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Thanuja/input/)
    FILE(COPY           NumPdeSpring2015/Thanuja/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Thanuja/input/)
  
  ADD_SUBDIRECTORY(NumPdeSpring2015/Thanuka/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Thanuka/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Thanuka/input/)
    FILE(COPY           NumPdeSpring2015/Thanuka/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Thanuka/input/)
  
  ADD_SUBDIRECTORY(NumPdeSpring2015/Bhagya/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Bhagya/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Bhagya/input/)
    FILE(COPY           NumPdeSpring2015/Bhagya/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Bhagya/input/)
  
  ADD_SUBDIRECTORY(NumPdeSpring2015/Guoyi/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Guoyi/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Guoyi/input/)
    FILE(COPY           NumPdeSpring2015/Guoyi/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/Guoyi/input/)
    
  ADD_SUBDIRECTORY(NumPdeSpring2015/tutorial/Ex1)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/tutorial/Ex1/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/tutorial/Ex1/input/)
    FILE(COPY           NumPdeSpring2015/tutorial/Ex1/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/tutorial/Ex1/input/)
    
  ADD_SUBDIRECTORY(NumPdeSpring2015/tutorial/Ex2)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/tutorial/Ex2/output/)
    FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/tutorial/Ex2/input/)
    FILE(COPY           NumPdeSpring2015/tutorial/Ex2/input/ DESTINATION ${PROJECT_BINARY_DIR}/applications/NumPdeSpring2015/tutorial/Ex2/input/)  

ENDIF(BUILD_NUM_PDE_CLASS)


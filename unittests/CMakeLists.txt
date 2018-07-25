#############################################################################################
### Unit tests
#############################################################################################

ADD_SUBDIRECTORY(testNSSteadyDD/)
    
ADD_SUBDIRECTORY(testFSISteady/)

ADD_SUBDIRECTORY(testSalomeIO/)

IF(SLEPC_FOUND)
 ADD_SUBDIRECTORY(testSVD2NormCondNumb/)
ENDIF(SLEPC_FOUND)
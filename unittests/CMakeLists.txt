#############################################################################################
### Unit tests
#############################################################################################

ADD_SUBDIRECTORY(testNSSteadyDD/)
    
ADD_SUBDIRECTORY(testFSISteady/)

ADD_SUBDIRECTORY(testMED_IO/)

IF(SLEPC_FOUND)
 ADD_SUBDIRECTORY(testSVD2NormCondNumb/)
ENDIF(SLEPC_FOUND)

#############################################################################################
### Unit tests
#############################################################################################

ADD_SUBDIRECTORY(testNSSteadyDD/)
#     FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/unittests/testNSSteadyDD/input/)
#     FILE(COPY           testNSSteadyDD/input/ DESTINATION  ${PROJECT_BINARY_DIR}/unittests/testNSSteadyDD/input/)
    
ADD_SUBDIRECTORY(testFSISteady/)
#     FILE(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/unittests/testFSISteady/input/)
#     FILE(COPY           testFSISteady/input/ DESTINATION  ${PROJECT_BINARY_DIR}/unittests/testFSISteady/input/)
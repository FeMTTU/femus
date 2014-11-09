CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(TestNSSteadyDD)

femusMacroBuildUnitTest(testNSSteadyDD testNSSteadyDD)

#ADD_TEST(NAME testNSSteadyDD COMMAND valgrind --xml=yes --xml-file=${CMAKE_BINARY_DIR}/unittests/testNSSteadyDD.xml ./testNSSteadyDD)


 

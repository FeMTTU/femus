#
# Program:   MULTIMOD APPLICATION FRAMEWORK (MAF)
# Module:    $RCSfile: PackagesMacro.cmake,v $
# Language:  CMake 1.2
# Date:      $Date: 2011-06-03 09:29:54 $
# Version:   $Revision: 1.8.10.1 $
#
# Description:
# Macro file for find and unpack a tar.gz package and to apply patches.

############################################################################
##	      macro to find and unpack a tar.gz package and to apply patches	##
##	      FIND_AND_UNPACK_PACKAGE 																				##
##        "INPUT package name" 																						##
##        "PACKAGE_DIR package dir" 																			##
##        "SOURCE_DIR dir where will be expanded by untar" 								##
############################################################################


 MACRO(FIND_AND_UNPACK_PACKAGE PACKAGE_NAME PACKAGE_DIR SOURCE_DIR)
 
 ##package .tgz or tar.gz 
 	FILE(GLOB TARFILE_NAME_FULL "${PACKAGE_DIR}/${PACKAGE_NAME}*gz")
 	IF (TARFILE_NAME_FULL)	
 		GET_FILENAME_COMPONENT (TARFILE ${TARFILE_NAME_FULL} NAME)
 		CONFIGURE_FILE("${TARFILE_NAME_FULL}" "${SOURCE_DIR}/${TARFILE}" COPYONLY IMMEDIATE)
	  IF(UNIX)
	    INCLUDE(${CMAKE_ROOT}/Modules/FindUnixCommands.cmake)
	    IF (TAR-NOTFOUND)
  			MESSAGE(SEND_ERROR "Cannot unpack ${PACKAGE_NAME} archive without tar!")
  		ELSE (TAR-NOTFOUND)
  			MESSAGE(STATUS "${PACKAGE_NAME}: unpacking... ${TARFILE_NAME_FULL} ")
				EXEC_PROGRAM("${TAR} xzf  ${TARFILE} " ${SOURCE_DIR})
			##devo controllare gli errori dell'exec o lascio che il tutto proceda(vedi var OUT e RETURN)??
			ENDIF (TAR-NOTFOUND)
		ENDIF(UNIX)
	   
		IF (WIN32)
	    #searching for decompress programs
	    SET (GNUWIN32_INSTALL_PATH ${MFL_SOURCE_DIR}/Extras)
		
  		#FIND_PROGRAM(GZIP  gzip  ${GNUWIN32_INSTALL_PATH}/bin )
		SET (GZIP ${GNUWIN32_INSTALL_PATH}/bin/gzip.exe)
  		#FIND_PROGRAM(TAR  tar  ${GNUWIN32_INSTALL_PATH}/bin )  
		SET (TAR ${GNUWIN32_INSTALL_PATH}/bin/tar.exe)
		
  		#FIND_PROGRAM(GTAR  gtar  ${GNUWIN32_INSTALL_PATH}/bin )	  		
  		IF(GZIP AND TAR)
	   		#use a batch file to run the unpacking program
	   		SET (UNPACKING_SCRIPT "unpacking.bat")
	   		SET(MY_COMMAND "${GZIP} -df < ${TARFILE} | ${TAR} xf - " )
	 			STRING(REGEX REPLACE "/" "\\\\" UNPACK_COMMAND ${MY_COMMAND})
				FILE(WRITE ${SOURCE_DIR}/${UNPACKING_SCRIPT} "")
				FILE(APPEND  ${SOURCE_DIR}/${UNPACKING_SCRIPT} ${UNPACK_COMMAND})
				MESSAGE(STATUS "${PACKAGE_NAME}: unpacking... ${TARFILE_NAME_FULL} ")
				EXEC_PROGRAM(${SOURCE_DIR}/${UNPACKING_SCRIPT} ${SOURCE_DIR})
				##devo controllare gli errori dell'exec o lascio che il tutto proceda (vedi var OUT e RETURN)??
				
				#use the gtar command (A VERY GOOD SOLUTION!!!??? But the gtar command is not working yet)	
				#EXEC_PROGRAM("${GTAR} -xzf ${TARFILE} " ${SOURCE_DIR})
	   	ELSE(GZIP AND TAR)
	   			MESSAGE(SEND_ERROR "Cannot unpack  ${PACKAGE_NAME} archive: gzip and tar commands not found!")
	    ENDIF(GZIP AND TAR)			
		ENDIF(WIN32)
	ELSE(TARFILE_NAME_FULL)
		MESSAGE(SEND_ERROR "Package ${PACKAGE_NAME} not found")
	ENDIF(TARFILE_NAME_FULL)
	
	MARK_AS_ADVANCED(GZIP TAR)
 ENDMACRO(FIND_AND_UNPACK_PACKAGE)
 
 

 
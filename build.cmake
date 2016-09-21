SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS    ON)


########################################################################
##
## Finding out BUILD_TYPE : default Debug, accepted Debug|Release
##
########################################################################
MESSAGE("===> Detecting BUILD_TYPE [Debug|Release]")
SET(BUILD_TYPE_KNOWN   OFF)

IF("" STREQUAL "${BUILD_TYPE}")
	SET(BUILD_TYPE "Debug")
ENDIF()

IF("Debug" STREQUAL "${BUILD_TYPE}")
	SET(BUILD_TYPE_KNOWN     ON)
ENDIF()

IF("Release" STREQUAL "${BUILD_TYPE}")
	SET(BUILD_TYPE_KNOWN     ON)
ENDIF()

IF(NOT BUILD_TYPE_KNOWN)
	MESSAGE( FATAL_ERROR "Unknown BUILD_TYPE=${BUILD_TYPE}" )
ENDIF()

########################################################################
##
## Finding out GENERATOR : default Unix Makefiles
##
########################################################################
MESSAGE("===> Detecting GENERATOR [CMake -G]")
SET(GENERATOR_KNOWN OFF)
SET(GENERATOR_MULTI OFF) #if one generator handles multiples build types
SET(COMPILERS_TRACE ON)  #if BUILD_PATH keeps trace of COMPILERS

IF( "" STREQUAL "${GENERATOR}" )
	SET(GENERATOR "Unix Makefiles")
ENDIF()

IF( "${GENERATOR}" STREQUAL "Unix Makefiles" )
	SET(GENERATOR_KNOWN ON)
ENDIF()

IF( "${GENERATOR}" STREQUAL "Ninja" )
	SET(GENERATOR_KNOWN ON)
ENDIF()

IF( "${GENERATOR}" MATCHES "Visual Studio.*" )
	SET(GENERATOR_KNOWN ON)
	SET(GENERATOR_MULTI ON)
	SET(COMPILERS "microsoft") #cl.exe in MSVS
	SET(COMPILERS_TRACE OFF)
ENDIF()

IF( "${GENERATOR}" MATCHES "CodeBlocks.*" )
	SET(GENERATOR_KNOWN ON)
ENDIF()

IF( "${GENERATOR}" STREQUAL "Xcode" )
	SET(GENERATOR_KNOWN ON)
	SET(GENERATOR_MULTI ON)
	SET(COMPILERS       clang) #AppleClang in Xcode
	SET(COMPILERS_TRACE OFF)
ENDIF()

IF(NOT GENERATOR_KNOWN)
	MESSAGE( FATAL_ERROR "Unknown GENERATOR=${GENERATOR}")
ENDIF()

########################################################################
##
## Finding out COMPILERS : default gnu
##
########################################################################
MESSAGE("===> Detecting COMPILERS")
SET(COMPILERS_KNOWN OFF)
IF( "" STREQUAL "${COMPILERS}")
	SET(COMPILERS "gnu")
ENDIF()

IF( "${COMPILERS}" STREQUAL "gnu" )
	SET(CC  "gcc${VERSION}")
	SET(CXX "g++${VERSION}")
	SET(COMPILERS_KNOWN ON)
ENDIF()

IF( "${COMPILERS}" STREQUAL "clang" )
	SET(CC  "clang${VERSION}")
	SET(CXX "clang++${VERSION}")
	SET(COMPILERS_KNOWN ON)
ENDIF()

IF( "${COMPILERS}" STREQUAL "intel" )
	SET(CC  "icc")
	SET(CXX "icpc")
	SET(COMPILERS_KNOWN ON)
	SET(VERSION "")
ENDIF()

IF( "${COMPILERS}" STREQUAL "microsoft" )
	SET(CC  "cl.exe")
	SET(CXX "cl.exe")
	SET(COMPILERS_KNOWN ON)
	SET(VERSION "")
ENDIF()

IF(NOT COMPILERS_KNOWN)
	MESSAGE( FATAL_ERROR "Unknown COMPILERS=${COMPILERS}")
ENDIF()	

########################################################################
##
## Finding out SOURCE
##
########################################################################
MESSAGE("===> Detecting SOURCE [default: ./src]")
IF("" STREQUAL "${SOURCE}")
	SET(SOURCE "src")
ENDIF()

########################################################################
##
## prepare the build path
##
########################################################################
IF(COMPILERS_TRACE)
	SET(BUILD_ROOT "forge/${COMPILERS}${VERSION}")
ELSE()
	SET(BUILD_ROOT "forge")
ENDIF()

STRING(REGEX REPLACE " " "" __GENPATH ${GENERATOR})
IF(GENERATOR_MULTI)
	SET(BUILD_PATH "${BUILD_ROOT}/${__GENPATH}")
ELSE()
	SET(BUILD_PATH "${BUILD_ROOT}/${__GENPATH}/${BUILD_TYPE}")
ENDIF()

MESSAGE("=======> GENERATOR  = ${GENERATOR}")
MESSAGE("=======> BUILD_TYPE = ${BUILD_TYPE}")
MESSAGE("=======> COMPILERS  = ${COMPILERS}")
MESSAGE("=======>        CC  = ${CC}")
MESSAGE("=======>        CXX = ${CXX}")
MESSAGE("=======> BUILD_PATH = ${BUILD_PATH}")
MESSAGE("=======> SOURCE     = ${SOURCE}")

MESSAGE("===> Creating ${BUILD_PATH}")
EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E make_directory ${BUILD_PATH})


########################################################################
##
## Generate The Project
##
########################################################################
MESSAGE("===> Configuring with ${CMAKE_COMMAND}")
MESSAGE("")
MESSAGE("-------------------------------------------------------------")
IF(GENERATOR_MULTI)
	EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
							-G${GENERATOR}
							-DCMAKE_C_COMPILER=${CC}
							-DCMAKE_CXX_COMPILER=${CXX}
							${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE}
							WORKING_DIRECTORY ${BUILD_PATH})
ELSE()
	EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND}
							-G${GENERATOR}
							-DCMAKE_C_COMPILER=${CC}
							-DCMAKE_CXX_COMPILER=${CXX}
							-DCMAKE_BUILD_TYPE=${BUILD_TYPE}
							${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE}
							WORKING_DIRECTORY ${BUILD_PATH})
ENDIF()
MESSAGE("-------------------------------------------------------------")	
MESSAGE("")
MESSAGE("===> Ready to be compiled with ${GENERATOR} for ${COMPILERS}${VERSION}")
MESSAGE("")

IF( NOT ("" STREQUAL "${TARGET}") )
	MESSAGE("===> Building ${TARGET}")
	SET(NATIVE_FLAGS "")
	IF( "${GENERATOR}" STREQUAL "Xcode")
		MESSAGE("=======> parallel build for Xcode")
		SET(NATIVE_FLAGS "-parallelizeTargets")
	ENDIF()
	
	IF( "${GENERATOR}" MATCHES ".*Makefile.*" )
		MESSAGE("=======> parallel build for ${GENERATOR}")
		SET(NATIVE_FLAGS "-j8")
	ENDIF()
	
	IF( "${GENERATOR}" MATCHES ".*Visual Studio.*" )
		MESSAGE("=======> parallel build for ${GENERATOR}")
		SET(NATIVE_FLAGS "-m")
	ENDIF()
	
	IF(GENERATOR_MULTI)
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build . --target ${TARGET} --config ${BUILD_TYPE} --use-stderr -- ${NATIVE_FLAGS}
						WORKING_DIRECTORY ${BUILD_PATH})
	ELSE()
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} --build . --target ${TARGET} --use-stderr -- ${NATIVE_FLAGS}
						WORKING_DIRECTORY ${BUILD_PATH})
	ENDIF()
ENDIF()


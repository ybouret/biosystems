SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS ON)
MESSAGE("")
MESSAGE( STATUS "cleaning up...")
SET(to_remove "")

FUNCTION(__ADD_PATTERN pattern)
	FILE(GLOB tmp ${pattern})
	SET(to_remove ${to_remove} ${tmp} PARENT_SCOPE)
ENDFUNCTION(__ADD_PATTERN)

MESSAGE( STATUS "collecting local temporary files")
__ADD_PATTERN(*.dat)
__ADD_PATTERN(src/lithium/*.dat)
__ADD_PATTERN(src/lithium/supmat/*~)
__ADD_PATTERN(bin/*)

IF(APPLE)
MESSAGE( STATUS "collecting MacOS dumps")
FILE(GLOB_RECURSE tmp .DS_Store)
SET(to_remove ${to_remove} ${tmp})
ENDIF()

#specific stuff
MESSAGE( STATUS "collecting gnuplot fit.log" )
FILE(GLOB_RECURSE tmp fit.log)
SET(to_remove ${to_remove} ${tmp})

FOREACH(item IN LISTS to_remove)
	IF(EXISTS "${item}")
		MESSAGE("removing ${item}")
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove -f ${item})
	ENDIF()
ENDFOREACH(item)

MESSAGE("")



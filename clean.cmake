SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS ON)
MESSAGE("")
MESSAGE("===> cleaning up...")
SET(to_remove "")

FUNCTION(__ADD_PATTERN pattern)
	FILE(GLOB tmp ${pattern})
	SET(to_remove ${to_remove} ${tmp} PARENT_SCOPE)
ENDFUNCTION(__ADD_PATTERN)

MESSAGE("=======> collecting local temporary files")
__ADD_PATTERN(*.dat)
__ADD_PATTERN(*.bin)
__ADD_PATTERN(*.png)
__ADD_PATTERN(*.jpg)

IF(APPLE)
MESSAGE("=======> collecting MacOS dumps")
FILE(GLOB_RECURSE tmp .DS_Store)
SET(to_remove ${to_remove} ${tmp})
ENDIF()

#specific stuff

FOREACH(item IN LISTS to_remove)
	MESSAGE("removing ${item}")
	EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove -f ${item})
ENDFOREACH(item)

MESSAGE("")



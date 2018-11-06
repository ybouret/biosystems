SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS ON)
MESSAGE("")
MESSAGE(STATUS "removing forge content")

FILE(GLOB to_purge forge/*)
FOREACH(item IN LISTS to_purge)
	IF(IS_DIRECTORY ${item})
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove_directory ${item})
	ELSE()
		EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove -f ${item})
	ENDIF()
ENDFOREACH(item)

MESSAGE("")


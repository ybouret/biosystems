IF(FALSE)

ADD_EXECUTABLE(transient transient.cpp transient.lua)
TARGET_LINK_YOCTO(transient mk lua++)

ADD_EXECUTABLE(diff_v1 diff_v1.cpp)
TARGET_LINK_YOCTO(diff_v1 mk)

ADD_EXECUTABLE(eval_diff eval_diff.cpp)
TARGET_LINK_YOCTO(eval_diff mk)

ENDIF()

ADD_EXECUTABLE(sim8 sim8.cpp)
TARGET_LINK_YOCTO(sim8 mk)

INSTALL( TARGETS sim8  RUNTIME DESTINATION bin )

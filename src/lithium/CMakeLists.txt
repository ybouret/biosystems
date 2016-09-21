ADD_EXECUTABLE(transient transient.cpp transient.lua)
TARGET_LINK_YOCTO(transient mk lua++)

ADD_EXECUTABLE(diff_v1 diff_v1.cpp)
TARGET_LINK_YOCTO(diff_v1 mk)


INSTALL( TARGETS transient diff_v1 RUNTIME DESTINATION bin )

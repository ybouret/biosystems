ADD_EXECUTABLE(transient transient.cpp transient.lua)
TARGET_LINK_YOCTO(transient mk lua++ lua)

INSTALL( TARGETS transient RUNTIME DESTINATION bin )

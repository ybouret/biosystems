ADD_EXECUTABLE(fit_delta fit_delta.cpp)
TARGET_LINK_YOCTO(fit_delta mk)
INSTALL( TARGETS fit_delta RUNTIME DESTINATION bin )

ADD_EXECUTABLE(fit_delta fit_delta.cpp)
TARGET_LINK_YOCTO(fit_delta mk)
INSTALL( TARGETS fit_delta RUNTIME DESTINATION bin )

ADD_EXECUTABLE(fit_ximax fit_ximax.cpp)
TARGET_LINK_YOCTO(fit_ximax mk)
INSTALL( TARGETS fit_ximax RUNTIME DESTINATION bin )

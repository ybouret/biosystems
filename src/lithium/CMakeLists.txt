ADD_EXECUTABLE(build-proton build-proton.cpp)
TARGET_LINK_Y(build-proton m)
INSTALL( TARGETS build-proton RUNTIME DESTINATION bin )

ADD_EXECUTABLE(intake intake.cpp)
TARGET_LINK_Y(intake m)
INSTALL( TARGETS intake RUNTIME DESTINATION bin )

ADD_EXECUTABLE(bmax bmax.cpp)
TARGET_LINK_Y(bmax m)
INSTALL( TARGETS bmax RUNTIME DESTINATION bin )


#ADD_EXECUTABLE(fit_ximax fit_ximax.cpp)
#TARGET_LINK_YOCTO(fit_ximax mk)
#INSTALL( TARGETS fit_ximax RUNTIME DESTINATION bin )

ADD_EXECUTABLE(build-proton build-proton.cpp)
TARGET_LINK_Y(build-proton m)
INSTALL( TARGETS build-proton RUNTIME DESTINATION bin )

ADD_EXECUTABLE(intake intake.cpp intake.lua)
TARGET_LINK_Y(intake lua++ m)
INSTALL( TARGETS intake RUNTIME DESTINATION bin )

ADD_EXECUTABLE(bmax bmax.cpp fcnli.hpp)
TARGET_LINK_Y(bmax m)
INSTALL( TARGETS bmax RUNTIME DESTINATION bin )

ADD_EXECUTABLE(alpha alpha.cpp)
TARGET_LINK_Y(alpha m)
INSTALL( TARGETS alpha RUNTIME DESTINATION bin )

INCLUDE("${YROOT}/share/y/fltk.cmake")

IF(FLTK_FOUND)
FLUID_UIC(sim)
ADD_EXECUTABLE(flsim flsim.cpp fcnli.hpp ${sim_SOURCES})
TARGET_LINK_Y(flsim fltk)
TARGET_LINK_FLTK(flsim)
INSTALL( TARGETS flsim RUNTIME DESTINATION bin )
ENDIF()

#ADD_EXECUTABLE(fit_ximax fit_ximax.cpp)
#TARGET_LINK_YOCTO(fit_ximax mk)
#INSTALL( TARGETS fit_ximax RUNTIME DESTINATION bin )

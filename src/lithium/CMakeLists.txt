
ADD_EXECUTABLE(linsim linsim.cpp fcnli.hpp linsim.lua)
TARGET_LINK_Y(linsim lua++ m)
INSTALL( TARGETS linsim RUNTIME DESTINATION bin )



ADD_EXECUTABLE(fit-hall fit-hall.cpp)
TARGET_LINK_Y(fit-hall m)
INSTALL( TARGETS fit-hall RUNTIME DESTINATION bin )


ADD_EXECUTABLE(li-tot li-tot.cpp)
TARGET_LINK_Y(li-tot m)
INSTALL( TARGETS li-tot RUNTIME DESTINATION bin )

IF(FALSE)

ADD_EXECUTABLE(build-proton build-proton.cpp)
TARGET_LINK_Y(build-proton m)
INSTALL( TARGETS build-proton RUNTIME DESTINATION bin )

ADD_EXECUTABLE(fit-h fit-h.cpp)
TARGET_LINK_Y(fit-h m lang)
INSTALL( TARGETS fit-h RUNTIME DESTINATION bin )

ADD_EXECUTABLE(alpha alpha.cpp)
TARGET_LINK_Y(alpha m)
INSTALL( TARGETS alpha RUNTIME DESTINATION bin )

ADD_EXECUTABLE(bmax bmax.cpp fcnli.hpp)
TARGET_LINK_Y(bmax m)
INSTALL( TARGETS bmax RUNTIME DESTINATION bin )



ADD_EXECUTABLE(intake intake.cpp intake.lua)
TARGET_LINK_Y(intake lua++ m)
INSTALL( TARGETS intake RUNTIME DESTINATION bin )



ADD_EXECUTABLE(linear linear.cpp fcnli.hpp linear.lua)
TARGET_LINK_Y(linear lua++ m)
INSTALL( TARGETS linear RUNTIME DESTINATION bin )

INCLUDE("${YROOT}/share/y/fltk.cmake")

IF(FLTK_FOUND)
FLUID_UIC(sim)
ADD_EXECUTABLE(flsim flsim.cpp fcnli.hpp ${sim_SOURCES})
TARGET_LINK_Y(flsim fltk)
TARGET_LINK_FLTK(flsim)
INSTALL( TARGETS flsim RUNTIME DESTINATION bin )
ENDIF()

ENDIF(FALSE)


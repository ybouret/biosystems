


ADD_EXECUTABLE(simalpha simalpha.cpp simalpha.lua)
TARGET_LINK_Y(simalpha lua++ m)
INSTALL( TARGETS simalpha RUNTIME DESTINATION bin )


ADD_EXECUTABLE(simpot simpot.cpp)
TARGET_LINK_Y(simpot m)
INSTALL( TARGETS simpot RUNTIME DESTINATION bin )

ADD_EXECUTABLE(fit-leak fit-leak.cpp)
TARGET_LINK_Y(fit-leak lua++ m)
INSTALL( TARGETS fit-leak RUNTIME DESTINATION bin )


ADD_EXECUTABLE(fit-ph fit-ph.cpp)
TARGET_LINK_Y(fit-ph m)
INSTALL( TARGETS fit-ph RUNTIME DESTINATION bin )

ADD_EXECUTABLE(fit-relax fit-relax.cpp)
TARGET_LINK_Y(fit-relax m)
INSTALL( TARGETS fit-relax RUNTIME DESTINATION bin )

ADD_EXECUTABLE(fit-ph-delta fit-ph-delta.cpp)
TARGET_LINK_Y(fit-ph-delta m)
INSTALL( TARGETS fit-ph-delta RUNTIME DESTINATION bin )

IF(FALSE)
ADD_EXECUTABLE(eta eta.cpp li-common.hpp)
TARGET_LINK_Y(eta  m)
INSTALL( TARGETS eta RUNTIME DESTINATION bin )

ADD_EXECUTABLE(simlin simlin.cpp li-common.hpp)
TARGET_LINK_Y(simlin lua++ m)
INSTALL( TARGETS simlin RUNTIME DESTINATION bin )


ADD_EXECUTABLE(fullsim fullsim.cpp fcnli.hpp fullsim.lua)
TARGET_LINK_Y(fullsim lua++ m)
INSTALL( TARGETS fullsim RUNTIME DESTINATION bin )

ADD_EXECUTABLE(fit-hall fit-hall.cpp)
TARGET_LINK_Y(fit-hall m)
INSTALL( TARGETS fit-hall RUNTIME DESTINATION bin )

ADD_EXECUTABLE(build-proton build-proton.cpp)
TARGET_LINK_Y(build-proton m)
INSTALL( TARGETS build-proton RUNTIME DESTINATION bin )

INCLUDE("${YROOT}/share/y/fltk.cmake")

IF(FLTK_FOUND)
FLUID_UIC(sim)
ADD_EXECUTABLE(flsim flsim.cpp fcnli.hpp ${sim_SOURCES})
TARGET_LINK_Y(flsim fltk)
TARGET_LINK_FLTK(flsim)
INSTALL( TARGETS flsim RUNTIME DESTINATION bin )
ENDIF()

ENDIF(FALSE)


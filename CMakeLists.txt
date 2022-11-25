
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(KinFit)

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS} ${CMAKE_CURRENT_SOURCE_DIR})


# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT REQUIRED)
if ( NOT ROOT_FOUND )
message(FATAL_ERROR "ROOT is required but missing !")
endif ( NOT ROOT_FOUND )


Set(INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/include
)

Include_Directories(${INCLUDE_DIRECTORIES})

FILE(GLOB SRCFILES source/KinFitter.cxx source/KFitParticle.cxx source/KFitRootAnalyzer.cxx source/KFitDecayBuilder.cxx source/KFitNeutralCandFinder.cxx source/KFitVertexFinder.cxx) ## source/hvertextools.cxx)
FILE(GLOB HFILES include/KinFitter.h include/KFitParticle.h include/KFitRootAnalyzer.h include/KFitDecayBuilder.h include/KFitNeutralCandFinder.h include/KFitVertexFinder.h) ## include/hvertextools.h)

message (STATUS "${BoldGreen}ROOT${ColourReset} at " ${CMAKE_CURRENT_SOURCE_DIR})

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

ROOT_GENERATE_DICTIONARY(KinFitDict ${HFILES} LINKDEF include/FitterLinkDef.h)
add_library(KinFit SHARED ${SRCFILES} KinFitDict.cxx)
target_link_libraries(KinFit ${ROOT_LIBRARIES})
set_target_properties(KinFit PROPERTIES CXX_STANDARD 17)

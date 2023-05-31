
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(
  KinFit
  VERSION 0.0.1
  LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT REQUIRED Core)
if(NOT ROOT_FOUND)
  message(FATAL_ERROR "ROOT is required but missing !")
endif(NOT ROOT_FOUND)

set(SRCFILES
    source/KinFitter.cxx source/KFitParticle.cxx source/KFitRootAnalyzer.cxx
    source/KFitDecayBuilder.cxx source/KFitDecayCandFinder.cxx
    source/KFitVertexFinder.cxx) # source/hvertextools.cxx)
set(HFILES
    include/KinFitter.h include/KFitParticle.h include/KFitRootAnalyzer.h
    include/KFitDecayBuilder.h include/KFitDecayCandFinder.h
    include/KFitVertexFinder.h) # include/hvertextools.h)

# cmake-format: off
root_generate_dictionary(KinFitDict ${HFILES}
  MODULE ${PROJECT_NAME}
  LINKDEF include/FitterLinkDef.h)
# cmake-format: on
add_library(KinFit SHARED ${SRCFILES} KinFitDict.cxx)
add_library(KinFit::KinFit ALIAS KinFit)
target_link_libraries(KinFit ROOT::Core )
target_include_directories(
  KinFit
  PUBLIC $<INSTALL_INTERFACE:include>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR}
             VERSION ${PROJECT_VERSION}
             PUBLIC_HEADER "${HFILES}")

# Export the package for use from the build-tree (this registers the build-tree
# with a global CMake-registry)
export(PACKAGE ${PROJECT_NAME})

# Install the export set for use with the install-tree
install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  EXPORT ${PROJECT_NAME}Targets
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE KinFit::
  DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  COMPONENT dev)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  Config.cmake.in ${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
  PATH_VARS CMAKE_INSTALL_INCLUDEDIR CMAKE_INSTALL_LIBDIR)

write_basic_package_version_file(
  ${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_CMAKEDIR})

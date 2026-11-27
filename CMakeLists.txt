
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(
  KinFit
  VERSION 0.0.1
  LANGUAGES CXX)

# Without a build type CMake compiles with no optimisation at all (-O0),
# which makes the matrix-heavy fit several times slower.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      RelWithDebInfo
      CACHE STRING "Build type (Debug, Release, RelWithDebInfo, MinSizeRel)"
            FORCE)
endif()

include(GNUInstallDirs)
set(CMAKE_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Locate the ROOT package and define a number of useful targets and variables.
find_package(ROOT REQUIRED Core)
if(NOT ROOT_FOUND)
  message(FATAL_ERROR "ROOT is required but missing !")
endif(NOT ROOT_FOUND)

# AddressSanitizer build for debugging memory errors (e.g. with the QA suite)
option(KINFIT_ENABLE_ASAN "Build with AddressSanitizer" OFF)
if(KINFIT_ENABLE_ASAN)
  add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address)
endif()

# Core fitter (KinFitter, particles, vertexing) and analysis layer
# (KFitAnalyzer, KFitDecayBuilder); headers and sources live side by side.
set(CORE_CLASSES
    KinFitter KinFitParticle KFitVertexFinder KFitVertexFinderConstB
    KFitDecayCandFinder CoordinateConversion)
set(ANALYSIS_CLASSES KFitAnalyzer KFitDecayBuilder)

set(SRCFILES)
set(HFILES)
set(DICT_HEADERS)
foreach(cls ${CORE_CLASSES})
  list(APPEND SRCFILES core/${cls}.cxx)
  list(APPEND HFILES core/${cls}.h)
  list(APPEND DICT_HEADERS ${cls}.h)
endforeach()
foreach(cls ${ANALYSIS_CLASSES})
  list(APPEND SRCFILES analysis/${cls}.cxx)
  list(APPEND HFILES analysis/${cls}.h)
  list(APPEND DICT_HEADERS ${cls}.h)
endforeach()

# The dictionary refers to headers by bare file name (as they are installed),
# so rootcling needs both source directories on its include path.
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/core
                    ${CMAKE_CURRENT_SOURCE_DIR}/analysis)

# cmake-format: off
root_generate_dictionary(KinFitDict ${DICT_HEADERS}
  MODULE ${PROJECT_NAME}
  LINKDEF FitterLinkDef.h)
# cmake-format: on
add_library(KinFit SHARED ${SRCFILES} KinFitDict.cxx)
add_library(KinFit::KinFit ALIAS KinFit)
target_link_libraries(KinFit ROOT::Core ROOT::Matrix ROOT::Physics ROOT::MathCore ROOT::RIO ROOT::Tree)
target_include_directories(
  KinFit
  PUBLIC $<INSTALL_INTERFACE:include>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/core>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/analysis>)

set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES SOVERSION ${PROJECT_VERSION_MAJOR}
             VERSION ${PROJECT_VERSION}
             PUBLIC_HEADER "${HFILES}")

# Standalone examples in docs/examples (off by default)
option(KINFIT_BUILD_EXAMPLES "Build the compiled examples in docs/examples" OFF)
if(KINFIT_BUILD_EXAMPLES)
  foreach(example toy_lambda_fit vertex_constB)
    add_executable(${example} docs/examples/${example}.cxx)
    target_link_libraries(${example} PRIVATE KinFit)
  endforeach()
endif()

# QA test suite (CTest), see QA/CMakeLists.txt
option(KINFIT_BUILD_TESTS "Build the QA test suite" ON)
if(KINFIT_BUILD_TESTS)
  enable_testing()
  add_subdirectory(QA)
endif()

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

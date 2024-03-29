cmake_minimum_required(VERSION 3.12)

project(DarkART VERSION "0.1.0")

# Use C++-14
set(CMAKE_CXX_STANDARD 14)

# Require (at least) it
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)

set(BIN_DIR ${PROJECT_SOURCE_DIR}/bin)
set(GENERATED_DIR ${PROJECT_BINARY_DIR}/generated)
set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(LIB_DIR ${PROJECT_SOURCE_DIR}/lib)
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
set(TESTS_DIR ${PROJECT_SOURCE_DIR}/tests)
set(EXTERNAL_DIR ${PROJECT_SOURCE_DIR}/external)

# External projects
find_library(ARB
  NAMES arb libarb flint-arb
  PATHS /usr/local/lib /usr/lib/ /usr/share/)
find_package(Boost 1.65 REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include(FetchContent)

# libphysica
set(LIBPHYSICA_DIR ${EXTERNAL_DIR}/libphysica)
FetchContent_Declare(libphysica
  GIT_REPOSITORY https://github.com/temken/libphysica.git
  GIT_TAG origin/dev
  SOURCE_DIR "${LIBPHYSICA_DIR}"
  CMAKE_ARGS "-DCMAKE_BUILD_TYPE=Release -DCODE_COVERAGE=OFF")
FetchContent_GetProperties(libphysica)

if(NOT libphysica_POPULATED)
  FetchContent_Populate(libphysica)
  add_subdirectory(${libphysica_SOURCE_DIR})
endif()

# When installing, CMake will clear the RPATH of these targets so they are installed with an empty RPATH.
# This setting avoids the clearing of the libconfig path.
set(LIBCONFIG_DIR ${EXTERNAL_DIR}/libphysica/external/libconfig)

if(IS_DIRECTORY "${LIBCONFIG_DIR}")
  message(STATUS "Libconfig was downloaded by libphysica. Setting RPATH manually for libconfig.")
  list(APPEND CMAKE_INSTALL_RPATH "${LIBCONFIG_DIR}/lib")
endif()

# Googletest
set(GTEST_DIR ${EXTERNAL_DIR}/googletest)
FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
  SOURCE_DIR "${GTEST_DIR}/src"
  BINARY_DIR "${GTEST_DIR}/build"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  TEST_COMMAND "")

# Version header
execute_process( # Git commit hash
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process( # Git commit hash
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process( # Git commit hash
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE CMAKE_GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process( # Git commit hash
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE CMAKE_GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
  ${INCLUDE_DIR}/DarkART/version.hpp.in
  ${GENERATED_DIR}/DarkART/version.hpp
)

# Source and include directories
include_directories(
  ${INCLUDE_DIR}
  ${Boost_INCLUDE_DIR})
add_subdirectory(${SRC_DIR})

# Code Coverage Configuration
if(NOT TARGET coverage_config)
  add_library(coverage_config INTERFACE)
endif()

option(CODE_COVERAGE "Enable coverage reporting" OFF)

if(CODE_COVERAGE)
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(coverage_config INTERFACE
    -O0 # no optimization
    -g # generate debug info
    --coverage # sets all required flags
  )

  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(coverage_config INTERFACE --coverage)
  else()
    target_link_libraries(coverage_config INTERFACE --coverage)
  endif()
endif(CODE_COVERAGE)

option(BUILD_TESTING "Build the testing tree." ON)

# Only build tests if we are the top-level project
# Allows this to be used by super projects with `add_subdirectory`
if(BUILD_TESTING AND(PROJECT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR))
  enable_testing()
  FetchContent_GetProperties(googletest)

  if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(
      ${googletest_SOURCE_DIR}
      ${googletest_BINARY_DIR}
      EXCLUDE_FROM_ALL)
  endif()

  add_subdirectory(${TESTS_DIR})
endif()
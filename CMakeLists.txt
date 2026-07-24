# ##############################################################################
#
# solid_dmft - An example application using triqs
#
# Copyright (C) ...
#
# solid_dmft is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# solid_dmft is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# solid_dmft (in the file COPYING.txt in this directory). If not, see
# <http://www.gnu.org/licenses/>.
#
# ##############################################################################

cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
cmake_policy(VERSION ${CMAKE_VERSION})

# ############
# Define Project
project(solid_dmft VERSION 4.0.0 LANGUAGES CXX)
get_directory_property(IS_SUBPROJECT PARENT_DIRECTORY)

# version displayed in the documentation
set(BRANCH_FOR_DOC ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.x)

# ############
# Load TRIQS
find_package(TRIQS 4.0 REQUIRED)

# Get the git hash & print status
triqs_get_git_hash_of_source_dir(PROJECT_GIT_HASH)
message(STATUS "${PROJECT_NAME} version : ${PROJECT_VERSION}")
message(STATUS "${PROJECT_NAME} Git hash: ${PROJECT_GIT_HASH}")

# Enforce Consistent Versioning
if(NOT ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR} VERSION_EQUAL ${TRIQS_VERSION_MAJOR}.${TRIQS_VERSION_MINOR})
  message(FATAL_ERROR "The ${PROJECT_NAME} version ${PROJECT_VERSION} is not compatible with TRIQS version ${TRIQS_VERSION}.")
endif()

# Default Install directory to TRIQS_ROOT if not given or when provided as relative path.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS "No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
  set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE)
endif()
if(NOT IS_SUBPROJECT)
  message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} --------")
endif()
set(${PROJECT_NAME}_BINARY_DIR ${PROJECT_BINARY_DIR} CACHE STRING "Binary directory of the ${PROJECT_NAME} Project")

# Make additional Find Modules available
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/share/cmake/Modules)

# ############
# CMake Options

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
if(NOT IS_SUBPROJECT)
  message(STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} --------")
endif()

# Python Support
if(NOT TRIQS_WITH_PYTHON_SUPPORT)
  message(FATAL_ERROR "TRIQS was installed without Python support. TRIQS Python support is required.")
endif()

# Documentation
option(Build_Documentation "Build documentation" OFF)

# Testing
option(Build_Tests "Build tests" ON)
if(Build_Tests)
  enable_testing()
  # When built as a dependency of a project configured with the LLVM sanitizers
  # (-DASAN/-DUBSAN), our Python tests load that project's instrumented TRIQS
  # libraries and need the sanitizer runtime preloaded. Locate it here so that
  # SANITIZER_RT_PRELOAD is populated before the test commands are registered.
  # No-op in a standalone build, where ASAN/UBSAN are unset.
  if(ASAN)
    find_package(sanitizer REQUIRED COMPONENTS asan)
  endif()
  if(UBSAN)
    find_package(sanitizer REQUIRED COMPONENTS ubsan)
  endif()
endif()

# Provide GNU Installation directories
include(GNUInstallDirs)

# #############
# Build Project

# Find / Build dependencies
add_subdirectory(deps)

# Tests
if(Build_Tests)
  add_subdirectory(test)
endif()

# Python
add_subdirectory(python/${PROJECT_NAME})

# Docs
if(NOT IS_SUBPROJECT AND Build_Documentation)
  add_subdirectory(doc)
endif()

# Additional configuration files
add_subdirectory(share)

# executable
add_subdirectory(bin)            # Executables

# packaging versions
add_subdirectory(packaging)

# #############
# Debian Package

option(BUILD_DEBIAN_PACKAGE "Build a deb package" OFF)
if(BUILD_DEBIAN_PACKAGE AND NOT IS_SUBPROJECT)
  if(NOT CMAKE_INSTALL_PREFIX STREQUAL "/usr")
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be /usr for packaging")
  endif()
  set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
  set(CPACK_GENERATOR "DEB")
  set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
  set(CPACK_PACKAGE_CONTACT "https://github.com/TRIQS/${PROJECT_NAME}")
  execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CMAKE_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "triqs (>= 4.0), triqs_cthyb (>= 4.0), triqs_dft_tools (>= 4.0), triqs_hartree_fock (>= 4.0), triqs_hubbardi (>= 4.0), triqs_maxent (>= 4.0)")
  set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
  set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)
  include(CPack)
endif()

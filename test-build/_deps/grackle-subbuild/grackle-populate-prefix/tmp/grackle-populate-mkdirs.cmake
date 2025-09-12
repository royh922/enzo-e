# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION}) # this file comes with cmake

# If CMAKE_DISABLE_SOURCE_CHANGES is set to true and the source directory is an
# existing directory in our source tree, calling file(MAKE_DIRECTORY) on it
# would cause a fatal error, even though it would be a no-op.
if(NOT EXISTS "/Users/juitenghsu/enzo-e/test-build/_deps/grackle-src")
  file(MAKE_DIRECTORY "/Users/juitenghsu/enzo-e/test-build/_deps/grackle-src")
endif()
file(MAKE_DIRECTORY
  "/Users/juitenghsu/enzo-e/test-build/_deps/grackle-build"
  "/Users/juitenghsu/enzo-e/test-build/_deps/grackle-subbuild/grackle-populate-prefix"
  "/Users/juitenghsu/enzo-e/test-build/_deps/grackle-subbuild/grackle-populate-prefix/tmp"
  "/Users/juitenghsu/enzo-e/test-build/_deps/grackle-subbuild/grackle-populate-prefix/src/grackle-populate-stamp"
  "/Users/juitenghsu/enzo-e/test-build/_deps/grackle-subbuild/grackle-populate-prefix/src"
  "/Users/juitenghsu/enzo-e/test-build/_deps/grackle-subbuild/grackle-populate-prefix/src/grackle-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/juitenghsu/enzo-e/test-build/_deps/grackle-subbuild/grackle-populate-prefix/src/grackle-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/juitenghsu/enzo-e/test-build/_deps/grackle-subbuild/grackle-populate-prefix/src/grackle-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()

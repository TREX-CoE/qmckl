#===========================================

# Try to find the QMCkl library;
# If found, it will define the following variables (note the plural form):
#  QMCKL_FOUND - System has libqmckl;
#  QMCKL_INCLUDE_DIRS - The QMCKL include directories;
#  QMCKL_LIBRARIES - The libraries needed to use QMCKL;

# If QMCKL has been installed in a non-standard location, one can set an 
# environment variable $QMCKL_DIR in the current shell:
# $ export QMCKL_DIR=<custom_path>
# to indicate the prefix used during the QMCKL installation
# (typically `./configure prefix=<custom_path> ..` or `cmake -DCMAKE_INSTALL_DIR=<custom_path> ..`)

# This file should be located WITHIN your project source tree.
# (e.g. in cmake/FindQMCKL.cmake)
# How to use it in your project CMakeLists.txt:

# This is needed to locate FindQMCKL.cmake file, modify it according to your source tree.
# list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

# find_package(QMCKL)
# if (QMCKL_FOUND)
#   include_directories(${QMCKL_INCLUDE_DIRS})
#   target_link_libraries(your_target ${QMCKL_LIBRARIES})
# endif()

#===========================================

# This file is distirbuted under the BSD 3-Clause License.
# Copyright (c) 2021, TREX Center of Excellence

#===========================================

message("<FindQMCKL.cmake>")

set(QMCKL_SEARCH_PATHS
	~/Library/Frameworks
	/Library/Frameworks
	/usr/local
	/usr
	/sw # Fink
	/opt/local # DarwinPorts
	/opt/csw # Blastwave
	/opt
)

find_path(QMCKL_INCLUDE_DIR 
	  NAMES qmckl.h 
	  HINTS $ENV{QMCKL_DIR}
	  PATH_SUFFIXES include/qmckl include
	  PATHS ${QMCKL_SEARCH_PATHS}
	  )
	        

# No need to specify platform-specific prefix (e.g. libqmckl on Unix) or
# suffix (e.g. .so on Unix or .dylib on MacOS) in NAMES. CMake takes care of that.
find_library(QMCKL_LIBRARY
             NAMES qmckl
	     HINTS $ENV{QMCKL_DIR}
	     PATH_SUFFIXES lib64 lib
	     PATHS ${QMCKL_SEARCH_PATHS}
	     )

message("<FindQMCKL.cmake>")

# Handle the QUIETLY and REQUIRED arguments and set QMCKL_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(QMCKL DEFAULT_MSG QMCKL_LIBRARY QMCKL_INCLUDE_DIR )
MARK_AS_ADVANCED(QMCKL_INCLUDE_DIR QMCKL_LIBRARY)

# Mot setting _INCLUDE_DIR and _LIBRARIES is considered a bug, 
# see https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/How-To-Find-Libraries
set(QMCKL_LIBRARIES ${QMCKL_LIBRARY})
set(QMCKL_INCLUDE_DIRS ${QMCKL_INCLUDE_DIR}) 


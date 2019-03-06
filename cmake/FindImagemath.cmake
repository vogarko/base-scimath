# - Try to find readline, a library for easy editing of command lines.
# Variables used by this module:
#  IMAGEMATH_ROOT_DIR     - Readline root directory
# Variables defined by this module:
#  IMAGEMATH_FOUND - system has IMAGEMATH Common
#  IMAGEMATH_INCLUDE_DIR  - the IMAGEMATH Common/include directory (cached)
#  IMAGEMATH_INCLUDE_DIRS - the IMAGEMATH Common include directories
#                          (identical to IMAGEMATH_INCLUDE_DIR)
#  IMAGEMATH_LIBRARY      - the IMAGEMATH common  library (cached)
#  IMAGEMATH_LIBRARIES    - the IMAGEMATH common library plus the libraries it 
#                          depends on

# Copyright (C) 2019


if(NOT IMAGEMATH_FOUND)

	find_path(IMAGEMATH_INCLUDE_DIR linmos/LinmosAccumulator.h
		HINTS ${IMAGEMATH_ROOT_DIR} PATH_SUFFIXES include/askap/imagemath)
	find_library(IMAGEMATH_LIBRARY askap_imagemath
		HINTS ${IMAGEMATH_ROOT_DIR} PATH_SUFFIXES lib)
	mark_as_advanced(IMAGEMATH_INCLUDE_DIR IMAGEMATH_LIBRARY )
	
	set(IMAGEMATH_INCLUDE_DIRS ${IMAGEMATH_INCLUDE_DIR})
	set(IMAGEMATH_LIBRARIES ${IMAGEMATH_LIBRARY})
        if(CMAKE_VERSION VERSION_LESS "2.8.3")
	   find_package_handle_standard_args(IMAGEMATH DEFAULT_MSG IMAGEMATH_LIBRARY IMAGEMATH_INCLUDE_DIR)
        else ()
	   include(FindPackageHandleStandardArgs)
	   find_package_handle_standard_args(IMAGEMATH DEFAULT_MSG IMAGEMATH_LIBRARY IMAGEMATH_INCLUDE_DIR)
        endif ()


endif(NOT IMAGEMATH_FOUND)

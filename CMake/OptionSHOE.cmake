# Copyright 2012 Sandia Corporation.
# Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the
# U.S. Government. Redistribution and use in source and binary forms, with
# or without modification, are permitted provided that this Notice and any
# statement of authorship are reproduced on all copies.

OPTION(BUILD_SHOE "Compile in support for Sandia's Higher Order Elements" OFF)
IF (BUILD_SHOE)
  OPTION ( BUILD_SHOE2 "Add support for the new higher order API" OFF )
  OPTION ( CXX_HAS_NAMESPACES "Does the C++ compiler support namespaces (including the 'using' directive)?" ON )
  INCLUDE ( ${VTKSNL_SOURCE_DIR}/CMake/FindPythonInterpreter.cmake )
  INCLUDE ( ${VTKSNL_SOURCE_DIR}/CMake/FindCln.cmake )
  INCLUDE ( ${VTKSNL_SOURCE_DIR}/CMake/FindGinac.cmake )
  INCLUDE ( ${VTKSNL_SOURCE_DIR}/CMake/FindGsl.cmake )
  INCLUDE( ${VTKSNL_SOURCE_DIR}/CMake/FindBlas.cmake )
  INCLUDE( ${VTKSNL_SOURCE_DIR}/CMake/FindLapack.cmake )

  IF ( FOUND_BLAS EQUAL 0 )
    MESSAGE( FATAL_ERROR "You must have BLAS. Please make sure that both\nBLAS_LIBRARY and (if you are on Linux)\nBLAS_G2C_LIBRARY are set. BLAS_G2C_LIBRARY\nshould point to libg2c, which is used by libblas." )
  ENDIF ( FOUND_BLAS EQUAL 0 )

  IF ( FOUND_LAPACK EQUAL 0 )
    MESSAGE( FATAL_ERROR "You must have LAPACK. Please make sure that both\nLAPACK_LIBRARY and (if you are on Linux)\nLAPACK_G2C_LIBRARY are set. LAPACK_G2C_LIBRARY\nshould point to libg2c, which is used by libblas." )
  ENDIF ( FOUND_LAPACK EQUAL 0)
ENDIF (BUILD_SHOE)


#
# Copyright 2012 Sandia Corporation.
# Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the
# U.S. Government. Redistribution and use in source and binary forms, with
# or without modification, are permitted provided that this Notice and any
# statement of authorship are reproduced on all copies.
#
SET( ShoeElementsLibSrcs
  ShapeFunctions.cxx
  CriticalPoints.cxx
  Generic.cxx
  Point.cxx
  Curve.cxx
  Triangle.cxx
  Quadrilateral.cxx
  Tetrahedron.cxx
  Hexahedron.cxx
  ${CMAKE_CURRENT_BINARY_DIR}/CellOpsImpl.cxx
  VelhoFunctions.cxx
)

SET_SOURCE_FILES_PROPERTIES( CellOpsImpl.cxx PROPERTIES GENERATED 1 )

ADD_CUSTOM_COMMAND(
  SOURCE ${VTKSNL_SOURCE_DIR}/Common/Elements/CellOpsSpecifications.txt ${VTKSNL_SOURCE_DIR}/Common/Elements/GenCellOps.py
  COMMAND ${PYTHON}
  ARGS ${VTKSNL_SOURCE_DIR}/Common/Elements/GenCellOps.py ${VTKSNL_SOURCE_DIR}/Common/vtkCellOps.h ${VTKSNL_SOURCE_DIR}/Common/Elements/CellOpsSpecifications.txt ${VTKSNL_BINARY_DIR}/Common/Elements/CellOpsImpl.cxx
  OUTPUTS ${VTKSNL_BINARY_DIR}/Common/Elements/CellOpsImpl.cxx
  TARGET PythonInputFiles
  DEPENDS ${VTKSNL_SOURCE_DIR}/Common/Elements/CellOpsSpecifications.txt ${VTKSNL_SOURCE_DIR}/Common/Elements/GenCellOps.py ${VTKSNL_SOURCE_DIR}/Common/vtkCellOps.h # ${VTKSNL_BINARY_DIR}/Common/Elements/CellOpsImpl.cxx
)

ADD_LIBRARY( vtksnlShoeElementsD SHARED ${ShoeElementsLibSrcs} )
TARGET_LINK_LIBRARIES( vtksnlShoeElementsD
  vtksnlCommon
  Pss
  ${GSL_LIBRARY}
  ${GSL_CBLAS_LIBRARY}
  ${LAPACK_LIBRARY}
  ${LAPACK_EXTRA_LIBRARIES}
  vtkGraphics
  vtkFiltering
  vtkCommon
)
ADD_LIBRARY( vtksnlShoeElements MODULE vtksnlShoeElementsInit.cxx )
TARGET_LINK_LIBRARIES( vtksnlShoeElements vtksnlShoeElementsD )

INCLUDE_DIRECTORIES(
  ${GSL_INCLUDE_PATH}
)

IF ( FOUND_GINAC )
  TARGET_LINK_LIBRARIES( vtksnlShoeElementsD
    ${GINAC_LIBRARY}
    ${CLN_LIBRARY}
    ${CLN_EXTRA_LIBRARY}
  )
  INCLUDE_DIRECTORIES(
    ${GINAC_INCLUDE_PATH}
    ${CLN_INCLUDE_PATH}
  )
ENDIF ( FOUND_GINAC )

INCLUDE(${CMAKE_ROOT}/Modules/CheckIncludeFile.cmake)
CHECK_INCLUDE_FILE( malloc.h HAVE_MALLOC_H )

INSTALL_TARGETS ( /lib/vtksnl vtksnlShoeElements )

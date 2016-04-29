# Copyright 2012 Sandia Corporation.
# Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the
# U.S. Government. Redistribution and use in source and binary forms, with
# or without modification, are permitted provided that this Notice and any
# statement of authorship are reproduced on all copies.

ADD_EXECUTABLE( TessellatorQuality
  tessellatorQuality.cxx
  vtkToySubdivision.cxx
  vtkLissajousSubdivision.cxx
  )
IF (BUILD_SHOE2)
  ADD_EXECUTABLE( S3dViz
    s3dViz.cxx
    )
  TARGET_LINK_LIBRARIES( S3dViz
    vtkGenericFiltering
    )
  ADD_EXECUTABLE( benchmarkTessellators
    benchmarkTessellators.cxx
    )
  TARGET_LINK_LIBRARIES( benchmarkTessellators
    vtksnlCommon
    vtkGenericFiltering
    vtkFiltering
    vtkCommon
    )
ENDIF (BUILD_SHOE2)

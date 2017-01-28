#
# Source files
#

SET (LIBRARY_NAME
    mls_vtk
)


SET (SOURCES
    vtkMovingLeastSquare.cpp
)


SET (WRAPPED_LIBRARIES
    vtkCommon
    vtkFiltering
    vtkVolumeRendering
)


SET (LIBRARIES
    ${WRAPPED_LIBRARIES}
)


SET (INCLUDE_DIRS
    ${PROJECT_SOURCE_DIR}/${LIBRARY_NAME}
    ${PROJECT_SOURCE_DIR}/include
)



INCLUDE(../CMakeCommon.txt)



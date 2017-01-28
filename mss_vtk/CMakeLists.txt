#
# Source files
#

SET (LIBRARY_NAME
    mss_vtk
)

SET (SOURCES
    Distance.cpp
    UGridAccelerator.cpp
    vtkMassSpring.cpp
)

SET_SOURCE_FILES_PROPERTIES(
    Distance.cpp
    UGridAccelerator.cpp
  PROPERTIES
    WRAP_EXCLUDE 1
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




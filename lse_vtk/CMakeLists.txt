#
# Source files
#

SET (LIBRARY_NAME
    lse_vtk
)

SET (SOURCES
    vtkLaplacianSurface.cpp
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
)

if(COMMAND cmake_policy)
   cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
ADD_EXECUTABLE(test_LS test_LS.cpp)
TARGET_LINK_LIBRARIES(test_LS vtkRendering lse_vtk)

INCLUDE(../CMakeCommon.txt)



#
# Source files
#

SET (LIBRARY_NAME
    meshless_vtk
)

SET (SOURCES
    vtkMeshLess.cpp
    ${PROJECT_SOURCE_DIR}/lse_vtk/vtkLaplacianSurface.cpp
)

SET_SOURCE_FILES_PROPERTIES(
    ${PROJECT_SOURCE_DIR}/lse_vtk/vtkLaplacianSurface.cpp
 PROPERTIES
 	WRAP_EXCLUDE 1	
)

SET (WRAPPED_LIBRARIES
    vtkCommon
    vtkFiltering
    vtkVolumeRendering
    lse_vtk
)

SET (LIBRARIES
    ${WRAPPED_LIBRARIES}
)

SET (INCLUDE_DIRS
    ${PROJECT_SOURCE_DIR}/${LIBRARY_NAME}
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/lse_vtk
)

if(COMMAND cmake_policy)
   cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
ADD_EXECUTABLE(test_meshless test_meshless.cpp)
TARGET_LINK_LIBRARIES(test_meshless vtkRendering meshless_vtk)

INCLUDE(../CMakeCommon.txt)



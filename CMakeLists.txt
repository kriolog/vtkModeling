CMAKE_MINIMUM_REQUIRED(VERSION 2.0)

# Set the project/library name here.  Classes should be declared as
# "class VTKExtend_EXPORT vtkFoo", where VTKExtend is the name of the
# project set here.
PROJECT(VTKExtend)

SUBDIRS(PREORDER
    lse_vtk
    mls_vtk
    mss_vtk
    meshless_vtk
)

# Make sure the linker can find the matching_vtk library once it is built.
#LINK_DIRECTORIES(${VTKExtend_BINARY_DIR}/matching_vtk)

# Include all sub libraries into the result
SET(sub_libs
    lse_vtk
    mls_vtk
    mss_vtk
    meshless_vtk
)


# Set your list of sources here.  Do not change the name of the
# VTKExtend_SRCS variable.
SET(VTKExtend_SRCS
)


# List the kits from VTK that are needed by this project.
SET(VTKExtend_LIBS
  vtkCommon
)

#-----------------------------------------------------------------------------
# Most users should not need to change anything below this line.

# Need to include class headers and the configuration header.
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

# Setup output paths.
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH
    "Single output directory for building all libraries.")
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH
    "Single output directory for building all executables.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

# Find and load VTK settings.
FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# Give user option of building shared or static libraries.  Default to
# the choice made for VTK.
OPTION(BUILD_SHARED_LIBS "Build with shared libraries."
       ${VTK_BUILD_SHARED_LIBS})

# Set VTKExtend_BUILD_SHARED_LIBS to 1 or 0 for use in creating the
# configuration header.
IF(BUILD_SHARED_LIBS)
  SET(VTKExtend_BUILD_SHARED_LIBS_CONFIG 1)
ELSE(BUILD_SHARED_LIBS)
  SET(VTKExtend_BUILD_SHARED_LIBS_CONFIG 0)
ENDIF(BUILD_SHARED_LIBS)

OPTION(WRAP_PYTHON "WRAP_PYTHON."
		   ${WRAP_PYTHON})

project(LaplacianSurface)

if(ParaView_SOURCE_DIR)
  include_directories(${VTK_INCLUDE_DIRS})
else()
  find_package(ParaView REQUIRED)
  include(${PARAVIEW_USE_FILE})
endif()

add_paraview_plugin(
  ${PROJECT_NAME}SMPlugin "1.0"
  SERVER_MANAGER_XML ${PROJECT_NAME}.xml
  SERVER_MANAGER_SOURCES vtk${PROJECT_NAME}.cpp
)

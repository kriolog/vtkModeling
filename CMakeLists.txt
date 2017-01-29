if(NOT ParaView_SOURCE_DIR)
  cmake_minimum_required(VERSION 2.8)

  # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
endif()

include_directories(common)

add_subdirectory(Filters/LaplacianSurface)

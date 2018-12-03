cmake_minimum_required(VERSION 2.8)

PROJECT(StructuredGrid)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(StructuredGrid MACOSX_BUNDLE StructuredGrid.cpp)

if(VTK_LIBRARIES)
  target_link_libraries(StructuredGrid ${VTK_LIBRARIES})
else()
  target_link_libraries(StructuredGrid vtkHybrid vtkWidgets)
endif()

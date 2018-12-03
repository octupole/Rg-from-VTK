cmake_minimum_required(VERSION 2.8)

PROJECT(StructuredGrid)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(RgfromVTK RgfromVTK.cpp)

if(VTK_LIBRARIES)
  target_link_libraries(RgfromVTK ${VTK_LIBRARIES})
else()
  target_link_libraries(RgfromVTK vtkHybrid vtkWidgets)
endif()

cmake_minimum_required(VERSION 3.15)
project(vof2d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Options
option(VOF2D_BUILD_PY "Build Python extension via pybind11" ON)

include(FetchContent)

# Fetch pybind11
FetchContent_Declare(
  pybind11
  URL https://github.com/pybind/pybind11/archive/refs/tags/v2.12.0.zip
)
FetchContent_MakeAvailable(pybind11)

# Python dev is handled via pybind11; no explicit FindPython3 Development required here.

# Sources
file(GLOB VOF2D_HEADERS CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/vof2d/*.hpp
)
file(GLOB VOF2D_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)

add_library(vof2d_core STATIC ${VOF2D_SOURCES} ${VOF2D_HEADERS})
target_include_directories(vof2d_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

if(VOF2D_BUILD_PY)
  pybind11_add_module(vof2d MODULE python/bindings.cpp)
  target_link_libraries(vof2d PRIVATE vof2d_core)
  target_include_directories(vof2d PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
  # Place extension in a predictable dir
  set_target_properties(vof2d PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
  )
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
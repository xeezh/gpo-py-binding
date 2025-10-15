cmake_minimum_required(VERSION 3.18)
project(uav_py LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON) 
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(PYBIND11_FINDPYTHON ON)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../../ ${CMAKE_CURRENT_BINARY_DIR}/core_build)

add_library(uav_py MODULE module.cpp)

target_include_directories(uav_py PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/../../include
    ${CMAKE_CURRENT_LIST_DIR}/../../core
)
target_link_libraries(uav_py PRIVATE uav_core pybind11::module)
pybind11_extension(uav_py)
pybind11_strip(uav_py)

install(TARGETS uav_py
  LIBRARY DESTINATION uav_py
  RUNTIME DESTINATION uav_py
  ARCHIVE DESTINATION uav_py
)
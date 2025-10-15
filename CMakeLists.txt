cmake_minimum_required(VERSION 3.18)
project(UAV LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(FetchContent)
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.11.3
)
FetchContent_MakeAvailable(nlohmann_json)

file(GLOB CORE_SRC
  src/classes/*.cpp
  src/algorithms/*.cpp
  core/*.cpp
)

add_library(uav_core ${CORE_SRC})
target_compile_features(uav_core PUBLIC cxx_std_17)
target_include_directories(uav_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/core
)
target_link_libraries(uav_core PUBLIC nlohmann_json::nlohmann_json)

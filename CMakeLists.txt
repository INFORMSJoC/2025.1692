cmake_minimum_required(VERSION 3.11)

project(hpg C)
include(FetchContent)

option(CXX "enable C++ compilation" ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -march=native -pthread")

if(CXX)
  enable_language(CXX)
endif()

set(CMAKE_CXX_STANDARD 20)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if(MSVC AND MT)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()

find_package(GUROBI)

if(GUROBI_FOUND)
  message(STATUS "Gurobi found")
  include_directories(${GUROBI_INCLUDE_DIRS}) 

  add_executable(sol_gurobi.exe sol/main.cpp)
  target_link_libraries(sol_gurobi.exe optimized ${GUROBI_CXX_LIBRARY}
                          debug ${GUROBI_CXX_DEBUG_LIBRARY})
  target_link_libraries(sol_gurobi.exe ${GUROBI_LIBRARY})
  target_compile_definitions(sol_gurobi.exe PRIVATE USE_GUROBI)

  add_executable(gkp_gurobi.exe gkp/gkp.cpp)
  target_link_libraries(gkp_gurobi.exe optimized ${GUROBI_CXX_LIBRARY}
                          debug ${GUROBI_CXX_DEBUG_LIBRARY})
  target_link_libraries(gkp_gurobi.exe ${GUROBI_LIBRARY})
  target_compile_definitions(gkp_gurobi.exe PRIVATE USE_GUROBI)
else()
  message(STATUS "Gurobi not found, using highs")

  message(STATUS "Downloading highs...")
  FetchContent_Declare(
    highs
    GIT_REPOSITORY https://github.com/ERGO-Code/HiGHS.git
    GIT_TAG        "v1.7.2" # problem in getObjectiveValue in later versions
  )
  FetchContent_MakeAvailable(highs)

  add_executable(sol_highs.exe sol/main.cpp)
  target_link_libraries(sol_highs.exe PRIVATE highs)
  target_compile_definitions(sol_highs.exe PRIVATE USE_HIGHS)

  add_executable(gkp_highs.exe gkp/gkp.cpp)
  target_link_libraries(gkp_highs.exe PRIVATE highs)
  target_compile_definitions(gkp_highs.exe PRIVATE USE_HIGHS)
endif()

add_executable(enc.exe enc/encoder.cpp)

add_executable(viz.exe viz/main.cpp)
add_custom_command(TARGET viz.exe POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    "${CMAKE_SOURCE_DIR}/viz/visualizer_template.html"
    "$<TARGET_FILE_DIR:viz.exe>/visualizer_template.html"
)

# Uncomment to use soplex

# message(STATUS "Downloading soplex...")
# FetchContent_Declare(
#   soplex
#   GIT_REPOSITORY https://github.com/scipopt/soplex.git
#   GIT_TAG        master
# )
# 
# set(SOPLEX_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
# set(SOPLEX_BUILD_TESTS OFF CACHE BOOL "" FORCE)
# set(SOPLEX_USE_GMP OFF CACHE BOOL "" FORCE)
# FetchContent_MakeAvailable(soplex)
# 
# add_executable(sol_soplex.exe sol/main.cpp)
# target_link_libraries(sol_soplex.exe PRIVATE libsoplex)
# target_include_directories(sol_soplex.exe PRIVATE ${soplex_SOURCE_DIR}/src)
# target_compile_definitions(sol_soplex.exe PRIVATE USE_SOPLEX)
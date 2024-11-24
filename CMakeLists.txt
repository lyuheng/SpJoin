cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_COMPILER mpic++)

set(CMAKE_CXX_STANDARD 20)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -lpthread")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE) # Set to RELEASE later
endif()

project(SpatialJoin LANGUAGES CXX)

include_directories(include)
include_directories(.)

# NUMA related on Polaris
include_directories(/home/lyuheng/tools/numactl/include)
link_directories(/home/lyuheng/tools/numactl/lib)

# NUMA related on Cheaha
# include_directories(/home/lyuan/Desktop/numactl/include)
# link_directories(/home/lyuan/Desktop/numactl/lib)


find_package(OpenMP)

file(GLOB CAPIFiles CONFIGURE_DEPENDS "src/capi/*.h" "src/capi/*.cc")
add_library(capi STATIC ${CAPIFiles})

file(GLOB StorageFiles CONFIGURE_DEPENDS "src/storagemanager/*.h" "src/storagemanager/*.cc")
add_library(storagemanager STATIC ${StorageFiles})

file(GLOB SpatialFiles CONFIGURE_DEPENDS "src/spatialindex/*.h" "src/spatialindex/*.cc")
add_library(spatialindex STATIC ${SpatialFiles})

file(GLOB RTreeFiles CONFIGURE_DEPENDS "src/rtree/*.h" "src/rtree/*.cc")
add_library(rtree STATIC ${RTreeFiles})

file(GLOB ToolFiles CONFIGURE_DEPENDS "src/tools/*.h" "src/tools/*.cc")
add_library(tools STATIC ${ToolFiles})

file(GLOB EngineFiles CONFIGURE_DEPENDS "src/fine_grained/engine.cc" "src/fine_grained/engine.h")
add_library(engine STATIC ${EngineFiles})

add_executable(run run.cpp)
target_link_libraries(run rtree)
target_link_libraries(run capi)
target_link_libraries(run storagemanager)
target_link_libraries(run spatialindex)
target_link_libraries(run tools)
target_link_libraries(run numa)

add_executable(run_fg run_fg.cpp)
target_link_libraries(run_fg rtree)
target_link_libraries(run_fg spatialindex)
target_link_libraries(run_fg capi)
target_link_libraries(run_fg storagemanager)
target_link_libraries(run_fg tools)
target_link_libraries(run_fg engine)
target_link_libraries(run_fg numa)


add_executable(build_tree build_tree.cpp)
target_link_libraries(build_tree rtree)
target_link_libraries(build_tree capi)
target_link_libraries(build_tree storagemanager)
target_link_libraries(build_tree spatialindex)
target_link_libraries(build_tree tools)

add_executable(test_tree test_tree.cpp)
target_link_libraries(test_tree rtree)
target_link_libraries(test_tree capi)
target_link_libraries(test_tree storagemanager)
target_link_libraries(test_tree spatialindex)
target_link_libraries(test_tree tools)

add_executable(parse_data_tiger parse_data_tiger.cpp)
add_executable(parse_data_osm parse_data_osm.cpp)
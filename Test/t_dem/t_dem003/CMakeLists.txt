cmake_minimum_required(VERSION 3.3)

SET(cpp_file_name "t_dem_settling")

project(${cpp_file_name})

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX HL)
    if(HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIR})
        set(_hdf5_libs ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES})
    else()
        message(FATAL_ERROR "HDF5 not found!")    
    endif()

find_package(CGAL)
set(CMAKE_BUILD_TYPE "Release")

include_directories(~/ComFluSoM/Library/DEM)
include_directories(~/ComFluSoM/Library/Geometries)

add_executable(${cpp_file_name} t_dem_settling.cpp)
target_link_libraries (${cpp_file_name} Eigen3::Eigen OpenMP::OpenMP_CXX ${_hdf5_libs} CGAL::CGAL)

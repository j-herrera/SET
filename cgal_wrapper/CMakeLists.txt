cmake_minimum_required (VERSION 3.0)
project (SET VERSION 1.0.0 DESCRIPTION "Wrapper for C++ side of SET")

set(CMAKE_BUILD_TYPE Release) 

# Find CGAL
find_package(CGAL REQUIRED COMPONENTS Core)
include( ${CGAL_USE_FILE} )

# Find Python
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# Find Boost
find_package(Boost REQUIRED COMPONENTS system thread python3)
include_directories(${Boost_INCLUDE_DIRS})

# Find Eigen
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${Eigen3_INCLUDE_DIRS})

add_library(wrapper SHARED wrapper.cpp)
set_target_properties(wrapper PROPERTIES PREFIX "")

target_link_libraries(wrapper Eigen3::Eigen)
target_link_libraries(wrapper ${Boost_LIBRARIES})




cmake_minimum_required(VERSION 3.0.0)
project(wtpframe)
set(CMAKE_BUILD_TYPE Release)
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif() 

# add_subdirectory(include)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wno-deprecated-declarations") # -Wall -pedantic

include_directories( "./include" ) 
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
add_executable (wtpframe wtpframe.h wtpframe.cpp) 
target_link_libraries (wtpframe Eigen3::Eigen mpfr gmp) 

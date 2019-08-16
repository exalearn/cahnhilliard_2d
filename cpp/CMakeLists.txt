cmake_minimum_required(VERSION 2.8.12)
project(ch2d)

include(ExternalProject)

find_package(Boost COMPONENTS program_options REQUIRED)
find_package(Threads)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

add_definitions( -DSRCDIR=\"${CMAKE_SOURCE_DIR}/\")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Werror=reorder")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -DBOOST_LOG_DYN_LINK -DNDEBUG")


#############################################
# Build external dependencies
#############################################

set( EXTERNAL_INSTALL_DIR ${CMAKE_BINARY_DIR}/external_install )

include_directories( ${CMAKE_SOURCE_DIR}
                     ${Boost_INCLUDE_DIRS}
)

add_library(ch_src-object OBJECT
	     src/stochastic_euler.hpp
	     src/chparams.cpp
	     src/utils_ch.cpp
	     src/cahnhilliard.cpp
	     src/cahnhilliard_thermal.cpp
	     src/cahnhilliard_thermal_nodiffusion.cpp
	     src/run_ch_solver.cpp
	     src/right_hand_side.cpp  
  )
add_library(ch_src-static STATIC $<TARGET_OBJECTS:ch_src-object>)
set_target_properties(ch_src-static PROPERTIES OUTPUT_NAME ch_src)


add_executable( ch2d
  src/driver.cpp
  )


target_link_libraries( ch2d
  ${Boost_LIBRARIES}   
  ch_src-static
)


##############################################
# PyBind11/Python build
##############################################
add_subdirectory(external/pybind11)

FIND_PACKAGE( PythonLibs )
INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_PATH} )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )
INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/src )

add_library(ch_src-shared SHARED $<TARGET_OBJECTS:ch_src-object>)
set_target_properties(ch_src-shared PROPERTIES OUTPUT_NAME ch_src)


pybind11_add_module(aligned_vector src/python/aligned_vector/py_aligned_vector.cpp)
set_target_properties(aligned_vector PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/ch2d/)


pybind11_add_module(cahnhilliard src/python/cahnhilliard/py_cahnhilliard.cpp)
target_link_libraries(cahnhilliard PRIVATE ch_src-shared)
set_target_properties(cahnhilliard PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/python/ch2d/)


add_custom_target( mvcmd_py ALL
  COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/src/python/__init__.py" "${CMAKE_SOURCE_DIR}/python/ch2d/__init__.py"
  DEPENDS "${CMAKE_SOURCE_DIR}/src/python/__init__.py"
  COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/src/python/driver.py" "${CMAKE_SOURCE_DIR}/python/driver.py"
  DEPENDS "${CMAKE_SOURCE_DIR}/src/python/driver.py"
)

##############################################
## Swig/Python build
##############################################
#
#FIND_PACKAGE( SWIG REQUIRED )
#INCLUDE( ${SWIG_USE_FILE} )
#
#FIND_PACKAGE( PythonLibs )
#INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_PATH} )
#INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )
#INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR}/src )
#	    
#SET( CMAKE_SWIG_FLAGS "-threads" )
#SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fPIC")
#
#SET_SOURCE_FILES_PROPERTIES( src/cahnhilliard.i PROPERTIES CPLUSPLUS ON )
#SET_SOURCE_FILES_PROPERTIES( src/cahnhilliard.i PROPERTIES SWIG_FLAGS "-threads" )
#
##SWIG_ADD_LIBRARY( cahnhilliard TYPE SHARED LANGUAGE python SOURCES src/cahnhilliard.i )
#SWIG_ADD_MODULE( cahnhilliard python src/cahnhilliard.i src/chparams.cpp src/utils_ch.cpp )
#
#SWIG_LINK_LIBRARIES( cahnhilliard
#                     ${PYTHON_LIBRARIES}
#		     ch_src
#)
#		   
## Copy the shared objects and Swig wrappers to cpp/swig/
#add_custom_target( mvcmd_so ALL
#  COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/build/_cahnhilliard.so" "${CMAKE_SOURCE_DIR}/swig/_cahnhilliard.so"
#  DEPENDS "${CMAKE_SOURCE_DIR}/build/_cahnhilliard.so"
#)
#add_custom_target( mvcmd_py ALL
#  COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/build/cahnhilliard.py" "${CMAKE_SOURCE_DIR}/swig/cahnhilliard.py"
#  DEPENDS "${CMAKE_SOURCE_DIR}/build/cahnhilliard.py"
#)

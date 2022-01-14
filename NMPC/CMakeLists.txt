cmake_minimum_required(VERSION 3.10.2)
project(nmpc)
 
set(CMAKE_CXX_STANDARD 11)
 
find_package(PythonLibs 2.7 REQUIRED)

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

add_executable(mpc mpc.cpp)
target_link_libraries(mpc ${PYTHON_LIBRARIES})
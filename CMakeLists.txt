set(CMAKE_CXX_FLAGS "-fpermissive")
cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_CXX_FLAGS "-fpermissive")
if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 17)
endif ()
project(etraj)
# Set source directory
set(SOURCE_DIR "src/")
set(INCLUDE_DIR "include/")
# Tell CMake that headers are also in SOURCE_DIR
include_directories(${SOURCE_DIR})
set(SOURCES "${SOURCE_DIR}vector.cpp"
            "${SOURCE_DIR}matrix.cpp"
            "${SOURCE_DIR}grid.cpp")

include_directories("${INCLUDE_DIR}nanoflann"
                    "${INCLUDE_DIR}nanoflann/examples"
                    "${INCLUDE_DIR}nanoflann/include")

#	fpermissive flag for perseus


# Generate Python module
add_subdirectory(include/pybind11)
add_executable(example ${SOURCE_DIR}main.cpp)
pybind11_add_module(etraj ${SOURCES} "src/bindings.cpp")
target_link_libraries(etraj PRIVATE ${MY_LIBRARIES})

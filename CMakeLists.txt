#set(CMAKE_CXX_FLAGS "-fpermissive")
cmake_minimum_required(VERSION 2.8.12)
if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 17)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif ()
if(USE_LOG_LEVEL_DEBUG)
  add_definitions(-DLOG_LEVEL_DEBUG)
endif()

# log level debug precompiler command
# this can be changed in setup.py
option(USE_LOG_LEVEL_DEBUG "Using debug mode for loggers" OFF)



MESSAGE("CMAKE FLAGS: " ${CMAKE_CXX_FLAGS})
project(etraj LANGUAGES CXX C Fortran)
include(FortranCInterface)

FortranCInterface_VERIFY(CXX)

FortranCInterface_HEADER(
  fc_mangle.h
  MACRO_NAMESPACE "FC_"
  SYMBOLS DSCAL DGESV
)

# Set source directory
set(SOURCE_DIR "src/")
set(INCLUDE_DIR "include/")
# Tell CMake that headers are also in SOURCE_DIR
include_directories(${SOURCE_DIR})
set(SOURCES "${SOURCE_DIR}linalg/vector.cpp"
            "${SOURCE_DIR}linalg/matrix.cpp"
            "${SOURCE_DIR}grid/grid.cpp"
            "${SOURCE_DIR}grid/ugrid.cpp"
            "${SOURCE_DIR}grid/kdtree.cpp"
            "${SOURCE_DIR}utils.cpp"
            "${SOURCE_DIR}geometry/geometry.cpp"
            "${SOURCE_DIR}fields/field.cpp"
            "${SOURCE_DIR}fields/scalarfield.cpp"
            "${SOURCE_DIR}fields/scalarfieldexample.cpp"
            "${SOURCE_DIR}fields/vectorfield.cpp"
            "${SOURCE_DIR}interpolation/interpolator.cpp"
            "${SOURCE_DIR}interpolation/radialbasis.cpp"
            "${SOURCE_DIR}interpolation/diffeq.cpp"
            "${SOURCE_DIR}integration/integrator.cpp"
            "${SOURCE_DIR}log/log.cpp"
            "${SOURCE_DIR}system/dynamicalsystem.cpp"
            "${INCLUDE_DIR}monomial/monomial.cpp")

include_directories("${INCLUDE_DIR}nanoflann"
                    "${INCLUDE_DIR}nanoflann/examples"
                    "${INCLUDE_DIR}nanoflann/include"
                    "${INCLUDE_DIR}monomial"
                    "${SOURCE_DIR}fields"
                    "${SOURCE_DIR}geometry"
                    "${SOURCE_DIR}grid"
                    "${SOURCE_DIR}interpolation"
                    "${SOURCE_DIR}integration"
                    "${SOURCE_DIR}linalg"
                    "${SOURCE_DIR}system"
                    "${SOURCE_DIR}log"
                    "/usr/include/x86_64-linux-gnu"
                    "/usr/lib/x86_64-linux-gnu")

#  lapacke must be linked first in order to
#  use the newer functions of the form,
#  LAPACKE_xxxx, rather than the older style,
#  xxxx_, which required the use of workspaces.
set (CMAKE_SHARED_LINKER_FLAGS "-L/usr/lib/x86_64-linux-gnu/lapack
                                -llapacke
                                -llapack
                                -lblas
                                -lgfortran")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${CMAKE_SOURCE_DIR}/include/morse_cmake/modules/find"
    "${CMAKE_SOURCE_DIR}/include/morse_cmake/modules")
include(MorseInit)

MESSAGE(STATUS "CMake source directory: " ${CMAKE_SOURCE_DIR})
MESSAGE(STATUS "CMake Module path is: " ${CMAKE_MODULE_PATH})
#find_library(blas libblas.a /usr/local/lib)
#MESSAGE(STATUS "BLAS location is: " ${blas})

# find BLAS, CBLAS, LAPACK, LAPACKE, and OPENBLAS
find_package( BLAS REQUIRED )
find_package( CBLAS REQUIRED )
find_package( LAPACK REQUIRED )
find_package( LAPACKE REQUIRED )
find_package( Boost REQUIRED )

# find CBLAS library
# find_path(BLAS_INCLUDE_DIRS cblas.h
#   /usr/include
#   /usr/local/include
#   $ENV{BLAS_HOME}/include)
MESSAGE(STATUS "BLAS    location is: " ${BLAS_INCLUDE_DIRS})
MESSAGE(STATUS "CBLAS   location is: " ${BLAS_INCLUDE_DIRS})
MESSAGE(STATUS "LAPACK  location is: " ${LAPACK_INCLUDE_DIRS})
MESSAGE(STATUS "LAPACKE location is: " ${LAPACKE_INCLUDE_DIRS})
MESSAGE(STATUS "BOOST   location is: " ${Boost_INCLUDE_DIRS})


include_directories(${BLAS_INCLUDE_DIRS}
                    ${CBLAS_INCLUDE_DIRS}
                    ${LAPACK_INCLUDE_DIRS}
                    ${LAPACKE_INCLUDE_DIRS}
                    ${Boost_INCLUDE_DIRS})
#-L/usr/lib/x86_64-linux-gnu/openblas -llapacke -llapack -lblas -lgfortran

# Generate Python module
add_subdirectory(include/pybind11)
#add_executable(example ${SOURCE_DIR}main.cpp)
pybind11_add_module(etraj ${SOURCES} "src/bindings.cpp")
target_link_libraries(etraj PUBLIC ${MY_LIBRARIES}
                                   ${BLAS_LIBRARIES}
                                   ${CBLAS_LIBRARIES}
                                   ${LAPACK_LIBRARIES}
                                   ${LAPACKE_LIBRARIES}
                                   ${Boost_LIBRARIES})

add_subdirectory("docs")

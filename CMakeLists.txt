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

# unit test option
option(UNIT_TEST "Conducting unit tests" OFF)

# Include useful macros
FILE( GLOB filenames "cmake/Macros/*.cmake" )
FOREACH( filename ${filenames} )
  INCLUDE( ${filename} )
ENDFOREACH()

# Set up the project etraj
MESSAGE("CMAKE FLAGS: " ${CMAKE_CXX_FLAGS})
project(etraj LANGUAGES CXX C Fortran)

# Configure FORTRAN interface
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)
FortranCInterface_HEADER(
  fc_mangle.h
  MACRO_NAMESPACE "FC_"
  SYMBOLS DSCAL DGESV
)

#  Setup external cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${CMAKE_SOURCE_DIR}/include/morse_cmake/modules/find"
    "${CMAKE_SOURCE_DIR}/include/morse_cmake/modules")
include(MorseInit)
MESSAGE(STATUS "CMake source directory: " ${CMAKE_SOURCE_DIR})
MESSAGE(STATUS "CMake Module path is: " ${CMAKE_MODULE_PATH})

# Set source directory
set(SOURCE_DIR "src/")
# Set the include directory
set(INCLUDE_DIR "include/")
# Tell CMake that headers are also in SOURCE_DIR
include_directories(${SOURCE_DIR})
# Set the include directories for the library
include_directories("${INCLUDE_DIR}nanoflann"
                    "${INCLUDE_DIR}nanoflann/examples"
                    "${INCLUDE_DIR}nanoflann/include"
                    "${INCLUDE_DIR}spdlog/include/"
                    "${INCLUDE_DIR}monomial"
                    "${SOURCE_DIR}fields"
                    "${SOURCE_DIR}geometry"
                    "${SOURCE_DIR}grid"
                    "${SOURCE_DIR}interpolation"
                    "${SOURCE_DIR}integration"
                    "${SOURCE_DIR}linalg"
                    "${SOURCE_DIR}system"
                    "${SOURCE_DIR}log"
                    "${SOURCE_DIR}utilities"
                    "/usr/include/x86_64-linux-gnu"
                    "/usr/lib/x86_64-linux-gnu"
)
# Set the sources for the library
set(SOURCES "${SOURCE_DIR}linalg/vector.cpp"
            "${SOURCE_DIR}linalg/matrix.cpp"
            "${SOURCE_DIR}grid/grid.cpp"
            "${SOURCE_DIR}grid/ugrid.cpp"
            "${SOURCE_DIR}grid/kdtree.cpp"
            "${SOURCE_DIR}utilities/random.cpp"
            "${SOURCE_DIR}utilities/utilities.cpp"
            "${SOURCE_DIR}geometry/geometry.cpp"
            "${SOURCE_DIR}fields/field.cpp"
            "${SOURCE_DIR}fields/scalarfield.cpp"
            "${SOURCE_DIR}fields/scalarfieldexample.cpp"
            "${SOURCE_DIR}fields/vectorfield.cpp"
            "${SOURCE_DIR}interpolation/interpolator.cpp"
            "${SOURCE_DIR}interpolation/localtaylor.cpp"
            "${SOURCE_DIR}interpolation/radialbasis.cpp"
            "${SOURCE_DIR}interpolation/diffeq.cpp"
            "${SOURCE_DIR}integration/integrator.cpp"
            "${SOURCE_DIR}log/log.cpp"
            "${SOURCE_DIR}system/dynamicalsystem.cpp"
            "${INCLUDE_DIR}monomial/monomial.cpp"
)


#  Setup BLAS, CBLAS, LAPACK, LAPACKE, etc. external libraries

#  lapacke must be linked first in order to
#  use the newer functions of the form,
#  LAPACKE_xxxx, rather than the older style,
#  xxxx_, which required the use of workspaces.
set (CMAKE_SHARED_LINKER_FLAGS "-L/usr/lib/x86_64-linux-gnu/lapack
                                -llapacke
                                -llapack
                                -lblas
                                -lgfortran")

# find BLAS, CBLAS, LAPACK, LAPACKE, and OPENBLAS
find_package( BLAS REQUIRED )
find_package( CBLAS REQUIRED )
find_package( LAPACK REQUIRED )
find_package( LAPACKE REQUIRED )
# find_package( Boost REQUIRED )

MESSAGE(STATUS "BLAS    location is: " ${BLAS_INCLUDE_DIRS})
MESSAGE(STATUS "CBLAS   location is: " ${BLAS_INCLUDE_DIRS})
MESSAGE(STATUS "LAPACK  location is: " ${LAPACK_INCLUDE_DIRS})
MESSAGE(STATUS "LAPACKE location is: " ${LAPACKE_INCLUDE_DIRS})
# MESSAGE(STATUS "BOOST   location is: " ${Boost_INCLUDE_DIRS})

# Include the directories just found
include_directories(${BLAS_INCLUDE_DIRS}
                    ${CBLAS_INCLUDE_DIRS}
                    ${LAPACK_INCLUDE_DIRS}
                    ${LAPACKE_INCLUDE_DIRS}
)

# generate the library libetraj
add_library(libetraj ${SOURCES})
# generate a test main.cpp file
add_executable(main ${SOURCE_DIR}main.cpp)
target_link_libraries(
  main PUBLIC libetraj
  ${BLAS_LIBRARIES}
  ${CBLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${LAPACKE_LIBRARIES}
  ${Boost_LIBRARIES}
)

# Generate Python module
add_subdirectory(include/pybind11)
pybind11_add_module(etraj ${SOURCES} "src/bindings.cpp")
target_link_libraries(
  etraj PUBLIC libetraj
  ${BLAS_LIBRARIES}
  ${CBLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${LAPACKE_LIBRARIES}
  ${Boost_LIBRARIES}
)

# generate documentation
add_subdirectory("docs")

# generate tests
if(UNIT_TEST)
  enable_testing()
  add_subdirectory(tests)
endif()

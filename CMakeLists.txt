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
set(SOURCES "${SOURCE_DIR}vector.cpp"
            "${SOURCE_DIR}matrix.cpp"
            "${SOURCE_DIR}grid.cpp")

include_directories("${INCLUDE_DIR}nanoflann"
                    "${INCLUDE_DIR}nanoflann/examples"
                    "${INCLUDE_DIR}nanoflann/include"
                    "/usr/include/x86_64-linux-gnu")

#	fpermissive flag for perseus

set (CMAKE_SHARED_LINKER_FLAGS "-L/usr/lib/x86_64-linux-gnu/openblas
                                -llapacke
                                -llapack
                                -lblas
                                -lgfortran")

find_library(blas libblas.a /usr/local/lib)
MESSAGE(STATUS "BLAS location is: " ${blas})

# Blas library
find_package( BLAS REQUIRED )
find_package( LAPACK REQUIRED )



# find CBLAS library
find_path(BLAS_INCLUDE_DIRS cblas.h
  /usr/include
  /usr/local/include
  $ENV{BLAS_HOME}/include)
MESSAGE(STATUS "CBLAS.h location is: " ${BLAS_INCLUDE_DIRS})
include_directories(${BLAS_INCLUDE_DIR})
#-L/usr/lib/x86_64-linux-gnu/openblas -llapacke -llapack -lblas -lgfortran

# Generate Python module
add_subdirectory(include/pybind11)
#add_executable(example ${SOURCE_DIR}main.cpp)
pybind11_add_module(etraj ${SOURCES} "src/bindings.cpp")
target_link_libraries(etraj PUBLIC ${MY_LIBRARIES}
                                   ${BLAS_LIBRARIES}
                                   ${LAPACK_LIBRARIES})

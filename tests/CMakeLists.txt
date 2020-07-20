#
# IF( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
#   # Ignore warnings related to the GNU statement expressions if we
#   # compile with clang.
#   ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wno-gnu-statement-expression" )
# ENDIF()

# Optimizes for a useful debugging experience. While the test cases
# might run slightly faster if "real" optimizations were enabled, I
# like this flag for debugging.
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Og" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-g" )
#
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wall" )
ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Wextra" )

# IF( NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
#   ENABLE_IF_SUPPORTED( CMAKE_CXX_FLAGS "-Werror" )
# ENDIF()

# Unit tests
set(BINARY ${CMAKE_PROJECT_NAME}_tests)

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp)

set(SOURCES ${TEST_SOURCE})

add_executable(test_vector                vector/test_vector.cpp)
add_executable(speed_test_dot_product     vector/speed_test_dot_product.cpp)
add_executable(test_kdtree                grid/test_kdtree.cpp)

add_test(vector               test_vector)
add_test(speed_dot_product    speed_test_dot_product)
add_test(kdtree               test_kdtree)

target_link_libraries(test_vector PUBLIC libetraj
  ${BLAS_LIBRARIES}
  ${CBLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${LAPACKE_LIBRARIES}
  ${Boost_LIBRARIES}
)
target_link_libraries(speed_test_dot_product PUBLIC libetraj
  ${BLAS_LIBRARIES}
  ${CBLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${LAPACKE_LIBRARIES}
  ${Boost_LIBRARIES}
)
target_link_libraries(test_kdtree PUBLIC libetraj
  ${BLAS_LIBRARIES}
  ${CBLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${LAPACKE_LIBRARIES}
  ${Boost_LIBRARIES}
)

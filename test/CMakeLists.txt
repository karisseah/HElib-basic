include_directories(${NTL_INCLUDE_DIRS})
add_library(NTL STATIC IMPORTED)
set_target_properties(NTL PROPERTIES IMPORTED_LOCATION ${NTL_LIBRARIES})

include_directories(${HELIB_INCLUDE_DIRS})
add_library(HElib STATIC IMPORTED)
set_target_properties(HElib PROPERTIES IMPORTED_LOCATION ${HELIB_LIBRARIES})

include_directories(${GTEST_INCLUDE_DIRS})
add_library(gtest STATIC IMPORTED)
set_property(TARGET gtest PROPERTY IMPORTED_LOCATION ${GTEST_LIBRARY_DIR}/libgtest.a)
add_library(gtest_main STATIC IMPORTED)
set_property(TARGET gtest_main PROPERTY IMPORTED_LOCATION ${GTEST_LIBRARY_DIR}/libgtest_main.a)

add_executable(BootstrappingTest BootstrappingTest.cpp)
target_link_libraries(BootstrappingTest gtest_main gtest ${CMAKE_THREAD_LIBS_INIT} HElib NTL -lgmp -pthread)
add_test(BootstrappingTest BootstrappingTest --gtest_output=xml:bootstrap-test.xml)
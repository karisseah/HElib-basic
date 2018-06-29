set(SOURCE Baseline.cpp encoding.cpp MatrixUtility.cpp)
set(HEADERS ../include/encoding.h ../include/MatrixUtility.h)
add_executable(HElib-basic ${SOURCE} ${HEADERS})

include_directories(${NTL_INCLUDE_DIRS})
add_library(NTL STATIC IMPORTED)
set_target_properties(NTL PROPERTIES IMPORTED_LOCATION ${NTL_LIBRARIES})

include_directories(${HELIB_INCLUDE_DIRS})
add_library(HElib STATIC IMPORTED)
set_target_properties(HElib PROPERTIES IMPORTED_LOCATION ${HELIB_LIBRARIES})

target_link_libraries(HElib-basic HElib NTL -lgmp -pthread)
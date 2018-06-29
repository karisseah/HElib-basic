INCLUDE(ExternalProject)

# Google Test
ExternalProject_Add(
        build-gtest
        # URL https://github.com/google/googletest/archive/master.zip
        URL ${CMAKE_CURRENT_SOURCE_DIR}/googletest-24062017.zip
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> -DCMAKE_BUILD_TYPE=release
)
ExternalProject_Get_Property(build-gtest INSTALL_DIR)
set(GTEST_INCLUDE_DIRS ${INSTALL_DIR}/include CACHE PATH "GTest include directory")
set(GTEST_LIBRARY_DIR "${INSTALL_DIR}/lib" CACHE PATH "GTest library folder path")

# NTL
ExternalProject_Add(
        build-ntl
        # URL http://www.shoup.net/ntl/ntl-10.5.0.tar.gz
        URL ${CMAKE_CURRENT_SOURCE_DIR}/ntl-10.5.0.tar.gz
        CONFIGURE_COMMAND cd src && ./configure NTL_THREADS=on NTL_THREAD_BOOST=on NTL_GMP_LIP=on PREFIX=<INSTALL_DIR>
        BUILD_COMMAND cd src && make -j${CORES_COUNT}
        INSTALL_COMMAND cd src && make install
        BUILD_IN_SOURCE 1
)
ExternalProject_Get_Property(build-ntl INSTALL_DIR)
set(NTL_INCLUDE_DIRS ${INSTALL_DIR}/include CACHE PATH "NTL include directory")
set(NTL_LIBRARIES "${INSTALL_DIR}/lib/libntl.a" CACHE PATH "NTL library path")

# HElib
set(HElib_CFLAGS "-std=c++11 -g -O2 -I${NTL_INCLUDE_DIRS} -DFHE_THREADS -DFHE_DCRT_THREADS")
ExternalProject_Add(
        build-HElib
        DEPENDS build-ntl
        # URL https://github.com/shaih/HElib/archive/master.zip
        #URL ${CMAKE_CURRENT_SOURCE_DIR}/Helib-09052018.zip
        URL ${CMAKE_CURRENT_SOURCE_DIR}/HElib-24062017.zip
        CONFIGURE_COMMAND ""
        BUILD_COMMAND cd src && make -j${CORES_COUNT} CFLAGS=${HElib_CFLAGS}
        INSTALL_COMMAND ""
        BUILD_IN_SOURCE 1
        BUILD_BYPRODUCTS <SOURCE_DIR>/src/fhe.a
)
ExternalProject_Get_Property(build-HElib SOURCE_DIR)
set(HELIB_INCLUDE_DIRS ${SOURCE_DIR}/src CACHE PATH "HElib include directory")
set(HELIB_LIBRARIES "${SOURCE_DIR}/src/fhe.a" CACHE PATH "HElib library path")
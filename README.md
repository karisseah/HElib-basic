# HElib-basic

To build the project, in project root folder:
mkdir build
cd build
cmake ../
make

If you are using CLion, in project root folder:
cd cmake-build-debug
rm 3rdparty
cd ..
cp build/3rdparty/ /cmake-build-debug/

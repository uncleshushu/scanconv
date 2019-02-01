if not exist build\windows\debug mkdir build\windows\debug
cd build\windows\debug
cmake -D CMAKE_BUILD_TYPE=Debug ..\..\..\ 
cmake --build . --config Debug
cd ..\..\..\
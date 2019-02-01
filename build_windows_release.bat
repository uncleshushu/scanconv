if not exist build\windows\release mkdir build\windows\release
cd build\windows\release
cmake -D CMAKE_BUILD_TYPE=Release ..\..\..\
MSBuild ALL_BUILD.vcxproj /p:Configuration=Release /v:minimal
cd ..\..\..\
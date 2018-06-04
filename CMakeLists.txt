cmake_minimum_required(VERSION 3.0.0)
project(usi_itp_ocl C)

set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 0)
set(PROJECT_VERSION_PATCH 0)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED on)   # requires C99

option(DEBUG "Debug version" ON)
option(RELEASE "Release version" OFF)
option(TIMING "Enable timing" ON)
option(USE_OPENCL "Use OpenCL" ON)

if(RELEASE)
    set(DEBUG OFF)    
endif(RELEASE)

configure_file(config.h.in ${PROJECT_SOURCE_DIR}/config.h)

if(USE_OPENCL)
    find_package(OpenCL REQUIRED)
    include_directories(PUBLIC ${OpenCL_INCLUDE_DIRS})  # set include directories for all targets
endif(USE_OPENCL)

if(${CMAKE_GENERATOR} MATCHES "Visual Studio")
    add_compile_options(/wd4996)
    if(DEBUG)
        file(COPY image.dat DESTINATION "${CMAKE_BINARY_DIR}/Debug")
        file(COPY interpolate.cl DESTINATION "${CMAKE_BINARY_DIR}/Debug")
        file(COPY evaluate.py DESTINATION "${CMAKE_BINARY_DIR}/Debug")
    endif(DEBUG)
    if(RELEASE)
        file(COPY image.dat DESTINATION "${CMAKE_BINARY_DIR}/Release")
        file(COPY interpolate.cl DESTINATION "${CMAKE_BINARY_DIR}/Release")
        file(COPY evaluate.py DESTINATION "${CMAKE_BINARY_DIR}/Release")
    endif(RELEASE)
else()
    file(COPY image.dat DESTINATION ${CMAKE_BINARY_DIR})
    file(COPY interpolate.cl DESTINATION ${CMAKE_BINARY_DIR})
    file(COPY evaluate.py DESTINATION ${CMAKE_BINARY_DIR})
endif()

add_executable(${PROJECT_NAME} main.c us_img.c)     # set the target

if(USE_OPENCL)
    target_link_libraries(${PROJECT_NAME} ${OpenCL_LIBRARIES})    # set link libraries for all targets    
endif(USE_OPENCL)

if(UNIX)
    target_link_libraries(${PROJECT_NAME} m)
endif(UNIX)




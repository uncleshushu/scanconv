cmake_minimum_required(VERSION 3.0.0)
project(usi_itp_ocl C)

set(PROJECT_VERSION_MAJOR 1)
set(PROJECT_VERSION_MINOR 3)
set(PROJECT_VERSION_PATCH 0)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED on)   # requires C11
add_compile_options(-Wall)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(NDEBUG "Disable debug print" ON)
option(TIMING "Enable timing" ON)
option(USE_OPENCL "Use OpenCL" ON)

if(${CMAKE_BUILD_TYPE} MATCHES "Debug" )
    set(NDEBUG OFF)
    remove_definitions(-DNDEBUG)
endif()

configure_file(config.h.in ${PROJECT_SOURCE_DIR}/config.h)

if(USE_OPENCL)
    find_package(OpenCL REQUIRED)
    include_directories(PUBLIC ${OpenCL_INCLUDE_DIRS})  # set include directories for all targets
endif(USE_OPENCL)


if(${CMAKE_GENERATOR} MATCHES "Visual Studio")
    add_compile_options(/wd4996)
endif()

add_executable(${PROJECT_NAME} main.c us_img.c)     # set the target

if(USE_OPENCL)
    target_link_libraries(${PROJECT_NAME} ${OpenCL_LIBRARIES})    # set link libraries for all targets    
endif(USE_OPENCL)

if(UNIX)
    target_link_libraries(${PROJECT_NAME} m)
endif(UNIX)


set(files2copy
    image.dat
    interpolate.cl
    evaluate.py
    evaluate_cpu.py
    profile.py)

if(${CMAKE_GENERATOR} MATCHES "Visual Studio")
    if(NDEBUG)
        set(copy_dest ${CMAKE_BINARY_DIR}"${CMAKE_BINARY_DIR}/Release")
    else(NDEBUG)
        set(copy_dest ${CMAKE_BINARY_DIR}"${CMAKE_BINARY_DIR}/Debug")
    endif(NDEBUG)
else()
    set(copy_dest ${CMAKE_BINARY_DIR})
endif()

foreach(f ${files2copy})
    file(COPY ${f} DESTINATION ${copy_dest})
endforeach(f)


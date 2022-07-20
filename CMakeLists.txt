cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)
cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0072 NEW)
cmake_policy(SET CMP0074 NEW)

project(psdr_cuda VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if( WIN32 )
    add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS -DNDEBUG)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -Wall -fPIC")
endif()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})

if( WIN32 )
    set(ENOKI_DIR ${CMAKE_SOURCE_DIR}/ext_win64/enoki)
    set(PYBIND11_ROOT ${CMAKE_SOURCE_DIR}/ext_win64/pybind11)
    set(PYTHON_ROOT C:/ProgramData/Miniconda3)
else()
    if ( NOT ENOKI_DIR )
        message(FATAL_ERROR "ENOKI_DIR not set!")
    endif()
    if ( NOT PYTHON_INCLUDE_PATH )
        message(FATAL_ERROR "PYTHON_INCLUDE_PATH not set!")
    endif()
endif()

set(PSDR_CUDA_FILE psdr_cuda)
add_compile_definitions(PSDR_CUDA_FILE="${PSDR_CUDA_FILE}")

option(BUILD_SHARED_LIBS "Build shared libraries" ON)

OPTION(CUDA_REMOVE_GLOBAL_MEMORY_SPACE_WARNING
       "Suppress the \"Advisory: Cannot tell what pointer points to, assuming global memory space\" warning nvcc makes." ON)

option(CUDA_GENERATE_DEPENDENCIES_DURING_CONFIGURE
       "Generate dependencies during configure time instead of only during build time." OFF)

find_package(CUDA 11.0 REQUIRED)
mark_as_advanced(CLEAR CUDA_64_BIT_DEVICE_CODE)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

if( WIN32 )
    set(OptiX_INSTALL_DIR "${CMAKE_SOURCE_DIR}/ext_win64/optix" CACHE PATH "Path to OptiX installed location.")
else()
    set(OptiX_INSTALL_DIR "${CMAKE_SOURCE_DIR}/ext/optix" CACHE PATH "Path to OptiX installed location.")
endif()
find_package(OptiX REQUIRED)


set(CUDA_GENERATED_OUTPUT_DIR "${CMAKE_BINARY_DIR}/lib/ptx")

if( WIN32 )
    include_directories(
        include/
        ${OptiX_INCLUDE}
        ${CMAKE_CURRENT_SOURCE_DIR}/cuda
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_BINARY_DIR}/include
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CUDA_INCLUDE_DIRS}
        ${PYTHON_ROOT}/include
        ${ENOKI_DIR}/include
        ${PYBIND11_ROOT}/include
    )
    link_directories(
        ${ENOKI_DIR}/lib
        ${PYTHON_ROOT}/libs
    )
else()
    include_directories(
        "include/"
        "${OptiX_INCLUDE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/cuda"
        ${CMAKE_CURRENT_SOURCE_DIR}
        "${CMAKE_BINARY_DIR}/include"
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CUDA_INCLUDE_DIRS}
        ${PYTHON_INCLUDE_PATH}
        ${ENOKI_DIR}/include
    )
    link_directories(
        ${ENOKI_DIR}/build
    )
endif()

function(build_ptx target_name_base target_name_var)
    set( target_name ${target_name_base} )
    set( ${target_name_var} ${target_name} PARENT_SCOPE )

    source_group("PTX Files"  REGULAR_EXPRESSION ".+\\.ptx$")
    source_group("CUDA Files" REGULAR_EXPRESSION ".+\\.cu$")

    CUDA_GET_SOURCES_AND_OPTIONS(source_files cmake_options options ${ARGN})
    CUDA_WRAP_SRCS( ${target_name} PTX generated_files ${source_files} ${cmake_options} OPTIONS ${options} )

    add_executable(${target_name}
        ${source_files}
        ${generated_files}
        ${cmake_options}
    )

    set_target_properties( ${target_name} PROPERTIES
        COMPILE_DEFINITIONS OPTIX_SAMPLE_NAME_DEFINE=${target_name}
    )
endfunction()

build_ptx( ptx target_name
    ./cuda/${PSDR_CUDA_FILE}.cu
    ./cuda/${PSDR_CUDA_FILE}.cpp
    ./cuda/${PSDR_CUDA_FILE}.h
    OPTIONS -rdc true
)

set(PYPSDR_LIBRARIES enoki-cuda enoki-autodiff cuda)

set(PSDR_INCLUDE_DIR include/psdr)
set(PSDR_SOURCE_DIR src)

set(PSDR_SOURCE_FILES

    ${PSDR_INCLUDE_DIR}/psdr.h
    ${PSDR_INCLUDE_DIR}/constants.h
    ${PSDR_INCLUDE_DIR}/fwd.h
    ${PSDR_INCLUDE_DIR}/macros.h
    ${PSDR_INCLUDE_DIR}/object.h
    ${PSDR_INCLUDE_DIR}/types.h
    ${PSDR_INCLUDE_DIR}/utils.h
    ${PSDR_INCLUDE_DIR}/core/ray.h
    ${PSDR_INCLUDE_DIR}/core/frame.h
    ${PSDR_INCLUDE_DIR}/core/intersection.h
    ${PSDR_INCLUDE_DIR}/core/bitmap.h
    ${PSDR_SOURCE_DIR}/core/bitmap.cpp
    ${PSDR_INCLUDE_DIR}/core/bitmap_loader.h
    ${PSDR_SOURCE_DIR}/core/bitmap_loader.cpp
    ${PSDR_INCLUDE_DIR}/core/sampler.h
    ${PSDR_SOURCE_DIR}/core/sampler.cpp
    ${PSDR_INCLUDE_DIR}/core/pmf.h
    ${PSDR_SOURCE_DIR}/core/pmf.cpp
    ${PSDR_INCLUDE_DIR}/core/cube_distrb.h
    ${PSDR_SOURCE_DIR}/core/cube_distrb.cpp

    ${PSDR_INCLUDE_DIR}/core/AQ_distrb.h
    ${PSDR_SOURCE_DIR}/core/AQ_distrb.cpp

    ${PSDR_INCLUDE_DIR}/core/tinyexr.h
    ${PSDR_INCLUDE_DIR}/core/miniz.h
    ${PSDR_SOURCE_DIR}/core/miniz.cpp

    ${PSDR_INCLUDE_DIR}/core/records.h
    ${PSDR_INCLUDE_DIR}/core/warp.h
    ${PSDR_INCLUDE_DIR}/core/transform.h
    ${PSDR_INCLUDE_DIR}/bsdf/bsdf.h
    ${PSDR_INCLUDE_DIR}/bsdf/diffuse.h
    ${PSDR_SOURCE_DIR}/bsdf/diffuse.cpp
    ${PSDR_INCLUDE_DIR}/bsdf/ggx.h
    ${PSDR_SOURCE_DIR}/bsdf/ggx.cpp
    ${PSDR_INCLUDE_DIR}/bsdf/roughconductor.h
    ${PSDR_SOURCE_DIR}/bsdf/roughconductor.cpp
    ${PSDR_INCLUDE_DIR}/edge/edge.h
    ${PSDR_INCLUDE_DIR}/emitter/emitter.h
    ${PSDR_INCLUDE_DIR}/emitter/area.h
    ${PSDR_SOURCE_DIR}/emitter/area.cpp
    ${PSDR_INCLUDE_DIR}/emitter/envmap.h
    ${PSDR_SOURCE_DIR}/emitter/envmap.cpp
    ${PSDR_INCLUDE_DIR}/sensor/sensor.h
    ${PSDR_SOURCE_DIR}/sensor/sensor.cpp
    ${PSDR_INCLUDE_DIR}/sensor/perspective.h
    ${PSDR_SOURCE_DIR}/sensor/perspective.cpp
    ${PSDR_INCLUDE_DIR}/shape/mesh.h
    ${PSDR_SOURCE_DIR}/shape/mesh.cpp
    ${PSDR_INCLUDE_DIR}/optix/ptx.h
    ${PSDR_SOURCE_DIR}/optix/ptx.cpp
    ${PSDR_INCLUDE_DIR}/scene/optix.h
    ${PSDR_INCLUDE_DIR}/scene/scene_optix.h
    ${PSDR_SOURCE_DIR}/scene/scene_optix.cpp
    ${PSDR_INCLUDE_DIR}/scene/scene.h
    ${PSDR_SOURCE_DIR}/scene/scene.cpp
    ${PSDR_INCLUDE_DIR}/scene/scene_loader.h
    ${PSDR_SOURCE_DIR}/scene/scene_loader.cpp
    ${PSDR_INCLUDE_DIR}/integrator/integrator.h
    ${PSDR_SOURCE_DIR}/integrator/integrator.cpp
    ${PSDR_INCLUDE_DIR}/integrator/field.h
    ${PSDR_SOURCE_DIR}/integrator/field.cpp
    ${PSDR_INCLUDE_DIR}/integrator/direct.h
    ${PSDR_SOURCE_DIR}/integrator/direct.cpp

    ${PSDR_INCLUDE_DIR}/bsdf/normalmap.h
    ${PSDR_SOURCE_DIR}/bsdf/normalmap.cpp

    ${PSDR_SOURCE_DIR}/psdr.cpp
)

set(CUDA_SOURCE_FILES
    ./cuda/kernel/util.cu
    ./cuda/kernel/NEE.cu
)

if( WIN32 )
    add_subdirectory(${PYBIND11_ROOT})
    cuda_add_library(cu_library STATIC ${CUDA_SOURCE_FILES})
    cuda_add_cublas_to_target(cu_library)
    pybind11_add_module(psdr_cuda ${PSDR_SOURCE_FILES})
    target_compile_options(psdr_cuda PRIVATE /O2 /wd4251 /MP)
    target_link_libraries(psdr_cuda PRIVATE cu_library ${PYPSDR_LIBRARIES} ${CUDA_LIBRARIES})
else()
    cuda_add_library(cu_library STATIC ${CUDA_SOURCE_FILES})
    cuda_add_cublas_to_target(cu_library)
    cuda_add_library(psdr_cuda MODULE ${PSDR_SOURCE_FILES})
    target_link_libraries(psdr_cuda cu_library ${PYPSDR_LIBRARIES} ${CUDA_LIBRARIES})
endif()


# target build
target_compile_definitions(psdr_cuda PRIVATE PTX_OUTPUT_DIR="${CUDA_GENERATED_OUTPUT_DIR}")

set_target_properties(psdr_cuda PROPERTIES SKIP_BUILD_RPATH FALSE)
set_target_properties(psdr_cuda PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
set_property(TARGET psdr_cuda PROPERTY CXX_STANDARD 17)
set_target_properties(psdr_cuda PROPERTIES PREFIX "")

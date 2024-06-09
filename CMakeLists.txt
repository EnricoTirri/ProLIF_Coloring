cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(ProLIF_Coloring VERSION 5.0)
enable_language(CXX)

#########################################################################
#### Tweak the default configuration
#########################################################################

# define an option for using native cpp code
set(TARGET_CPU_ARCHITECTURE "generic" CACHE STRING "The target cpu for architecture-specific optimizations")
set_property(CACHE TARGET_CPU_ARCHITECTURE PROPERTY STRINGS generic native)

# the idea is that we want to build the optimized version, unless the user says otherwise
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "The type of build" FORCE)
    message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified")
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif ()

# export also compile commands to have a better integration with the ide
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# enable additional warning flags
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -Wl,-z,defs -Wl,-z,now -Wl,-z,relro) #-Wpedantic
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    add_compile_options(-Wall -Wextra -Wshadow -Wdouble-promotion)
endif ()

# enable colored compiler output (with Ninja)
if ("${CMAKE_GENERATOR}" STREQUAL "Ninja")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        add_compile_options(-fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        add_compile_options(-fcolor-diagnostics)
    endif ()
endif ()

if (NOT DEFINED USECUDA)
    set(USECUDA 0)
endif ()

if (NOT DEFINED USEOMP)
    set(USEOPM 0)
endif ()

if(USECUDA AND USEOMP)
    message(FATAL_ERROR "You cannot build with both UseCuda and UseOpenMP")
endif ()


#########################################################################
#### Add external dependency
#########################################################################

# the rdkit to perform the heavy lifting
find_package(rdkit REQUIRED)

if (USEOMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    if (NOT DEFINED OMP_THREADS)
        set(OMP_THREADS 8)
    endif ()
    add_compile_definitions(OMP_NUM_THREADS=${OMP_THREADS})
endif ()

# define cuda options
if (USECUDA)
    enable_language(CUDA)
    if (NOT DEFINED CUDA_BLOCK_SIZE)
        set(CUDA_BLOCK_SIZE 128)
    endif ()
    add_compile_definitions(BLOCK_SIZE=${CUDA_BLOCK_SIZE})
endif ()


#########################################################################
#### Define the application sources
#########################################################################

# application headers
set(base_header_path "${CMAKE_CURRENT_SOURCE_DIR}/include-base")
set(extended_header_path "${CMAKE_CURRENT_SOURCE_DIR}/include-extended")

file(GLOB base_header_files "${base_header_path}/*.hpp")
file(GLOB extended_header_files "${extended_header_path}/*.hpp")

include_directories(${base_header_path})
include_directories(${extended_header_path})

# application sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
file(GLOB source_files "${source_path}/*.cpp")


if (USECUDA)
    file(GLOB impl_files "${CMAKE_CURRENT_SOURCE_DIR}/src-cuda/*.cu")
elseif (USEOMP)
    file(GLOB impl_files "${CMAKE_CURRENT_SOURCE_DIR}/src-omp/*.cpp")
else ()
    file(GLOB impl_files "${CMAKE_CURRENT_SOURCE_DIR}/src-normal/*.cpp")
endif ()

#########################################################################
#### Define the compilation step
#########################################################################

# define the target and its properties
add_executable(ProLIF_Coloring main.cpp ${base_header_files} ${extended_header_files} ${source_files} ${impl_files})

set_target_properties(ProLIF_Coloring
        PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
)

# define target properties if cuda is enabled
if (USECUDA)
    set_target_properties(ProLIF_Coloring
            PROPERTIES
            CMAKE_CUDA_STANDARD 17
            CUDA_SEPARABLE_COMPILATION ON
    )
endif ()

# enable link-time optimizations
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported)
if (ipo_supported)
    set_property(TARGET ProLIF_Coloring PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif ()

# since RDKit doesn't handle the dependencies in a correct way, we need to improvise
# NOTE: by using this way of importing RDKit, we need to manually include-base its dependencies
get_target_property(RDKIT_LIB_FULLPATH RDKit::RDGeneral LOCATION)
get_filename_component(RDKIT_LIB_DIRPATH "${RDKIT_LIB_FULLPATH}" DIRECTORY)
cmake_path(GET RDKIT_LIB_DIRPATH PARENT_PATH RDKIT_INSTALL_PREFIX)
target_include_directories(ProLIF_Coloring PUBLIC "${RDKIT_INSTALL_PREFIX}/include/rdkit")
target_link_directories(ProLIF_Coloring PUBLIC "${RDKIT_INSTALL_PREFIX}/lib")
target_link_directories(ProLIF_Coloring PUBLIC "${RDKIT_INSTALL_PREFIX}/lib64")
target_link_libraries(ProLIF_Coloring PUBLIC
        RDKitFileParsers
        RDKitGraphMol
        RDKitRDGeneral
        RDKitSmilesParse
        RDKitSubstructMatch
)
cmake_minimum_required(VERSION 3.13)

# Project definition:
# * Name: IChem
# * Version: 5.3.5
# * Languages: C and C++
project(IChem VERSION 5.3.5 LANGUAGES C CXX)

# CMake normally has NO build type set
# This forces a debug build
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type: Debug or Release" FORCE)
endif()


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# STATIC_ANALYSIS : enable/disable -fanalyzer and clang-tidy
# ICHEM_DEBUG : project debug
option(STATIC_ANALYSIS "Enable static analysis (clang-tidy / GCC -fanalyzer)" ON)
# -DICHEM_DEBUG=ON
option(ICHEM_DEBUG "Enable additional IChem debug code (ICHEM_DEBUG=1)" OFF)

# Debug : no optimization, full debug info
# Release : O3 + LTO + -march=native
set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto")

# Strict flags for development, can be commented when code is cleaner
set(STRICT_WARNING_FLAGS "")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(STRICT_WARNING_FLAGS
        -Wall
        -Wextra
        -Wpedantic
        -Wshadow
        -Wold-style-cast
        -Wzero-as-null-pointer-constant
        # We enable later:
        # -Wconversion
        # -Wsign-conversion
        # -Wnull-dereference
    )
endif()

# Static analysis configuration:
# If using GCC: append -fanalyzer to warning flags
# Clang-tidy is ran during the build if installed
# Toggle with: -DSTATIC_ANALYSIS=ON/OFF
if(STATIC_ANALYSIS)
    message(STATUS "Static analysis: ON")

    # GCC static analyzer
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(STATUS "  GCC: enabling -fanalyzer")
        list(APPEND STRICT_WARNING_FLAGS -fanalyzer)
    endif()

    # If Clang-Tidy found on the system
    find_program(CLANG_TIDY_EXE NAMES clang-tidy)
    if(CLANG_TIDY_EXE)
        message(STATUS "  clang-tidy: ${CLANG_TIDY_EXE}")
        set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE};-checks=clang-analyzer-*,modernize-*,bugprone-*,performance-*,readability-*,cppcoreguidelines-*,-clang-analyzer-alpha.*")
    else()
        message(STATUS "  clang-tidy is not found, skipping")
    endif()
else()
    message(STATUS "Static analysis: OFF")
endif()


set(SOURCES
    source/main.cpp
    source/sources/ICMole/coords.cpp
    source/sources/ICMole/edge.cpp
    source/sources/ICMole/graph.cpp
    source/sources/ICMole/vertex.cpp
    source/sources/ICMole/rigidbody.cpp
    source/sources/ICMole/residu.cpp
    source/sources/ICMole/physprop.cpp
    source/sources/ICMole/moletreat.cpp
    source/sources/ICMole/molecule.cpp
    source/sources/ICMole/global.cpp
    source/sources/ICMole/cycle.cpp
    source/sources/ICMole/complex.cpp
    source/sources/ICMole/chain.cpp
    source/sources/ICMole/bond.cpp
    source/sources/ICMole/atom.cpp
    source/sources/ICPars/molewriter.cpp
    source/sources/ICPars/molereader.cpp
    source/sources/ICPars/MOL2Parser.cpp
    source/sources/ICMole/grid.cpp
    source/sources/ICMole/fingerprint.cpp
    source/sources/ICMole/box.cpp
    source/sources/ICPars/PDBParser.cpp
    source/sources/ICCalcs/volsite.cpp
    source/sources/ICCalcs/pdbconvert.cpp
    source/sources/ICCalcs/interaction.cpp
    source/sources/ICRSA/RSA.cpp
    source/sources/ICRSA/PrimeGenerator.cpp
    source/sources/ICRSA/KeyPair.cpp
    source/sources/ICRSA/Key.cpp
    source/sources/ICRSA/EasyCrypt.cpp
    source/sources/ICRSA/BigInt.cpp
    source/sources/ICTools/switch.cpp
    source/sources/ICTools/atom_properties.cpp
    source/sources/ICTools/realign.cpp
    source/sources/ICTools/genKey.cpp
    source/sources/ICTools/IFP.cpp
    source/sources/ICTools/ifp_module.cpp
    source/sources/ICTools/BSAcalc.cpp
    source/sources/ICCalcs/grim.cpp
    source/sources/ICMole/similarity.cpp
    source/sources/ICMole/moleperceiv.cpp
    source/sources/ICMole/graphmatch.cpp
    source/sources/ICTools/ints.cpp
    source/sources/ICTools/grimcalc.cpp
    source/sources/ICCalcs/fragments.cpp
    source/sources/ICSVM/svm_main.cpp
    source/sources/ICSVM/svm_learn.c
    source/sources/ICSVM/svm_hideo.c
    source/sources/ICSVM/svm_common.c
    source/sources/ICTools/volcalc.cpp
    source/sources/ICTools/convert.cpp
    source/sources/ICTools/utils.cpp
    source/sources/ICTools/patch.cpp
    source/sources/ICTools/sims.cpp
    source/sources/ICTools/scoring.cpp
    source/sources/ICMole/pharmprop.cpp
    source/sources/ICPars/pharmwriter.cpp
    source/sources/ICPars/CHMParser.cpp
    source/sources/ICTools/phacalc.cpp
    source/sources/ICTools/convert_cav.cpp
    source/sources/ICTools/test.cpp
)

# Executable target IChem
add_executable(IChem ${SOURCES})

target_include_directories(IChem
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/source
        ${CMAKE_CURRENT_SOURCE_DIR}/source/headers/ICMole
        ${CMAKE_CURRENT_SOURCE_DIR}/source/headers/ICTools
        ${CMAKE_CURRENT_SOURCE_DIR}/source/headers/ICPars
        ${CMAKE_CURRENT_SOURCE_DIR}/source/headers/ICCalcs
        ${CMAKE_CURRENT_SOURCE_DIR}/source/headers/ICRSA
        ${CMAKE_CURRENT_SOURCE_DIR}/source/headers/ICSVM
        ${CMAKE_CURRENT_SOURCE_DIR}/source/external/nanoflann
)

# Enable with -DICHEM_DEBUG=ON
if(ICHEM_DEBUG)
    target_compile_definitions(IChem PRIVATE ICHEM_DEBUG=1)
endif()


# We have c files (SVM code), restrict STRICT_WARNING to C++ compilation units
# No strict warnings for now
# if(STRICT_WARNING_FLAGS)
#     target_compile_options(IChem PRIVATE
#         $<$<COMPILE_LANGUAGE:CXX>:${STRICT_WARNING_FLAGS}>
#     )
# endif()


# Treat warnings as errors in the debug build
# target_compile_options(IChem PRIVATE
#     $<$<AND:$<CONFIG:Debug>,$<COMPILE_LANGUAGE:CXX>>:-Werror>
# )

# find_package(Threads REQUIRED)
# target_link_libraries(IChem PRIVATE Threads::Threads)

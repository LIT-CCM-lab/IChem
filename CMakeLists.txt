cmake_minimum_required(VERSION 3.13)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -w")

# Optimization 
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto")

if (ICHEM_DEBUG)
    add_definitions(-DICHEM_DEBUG=1)
endif (ICHEM_DEBUG)

project(IChem VERSION 5.3.0)

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

include_directories(
    source
    source/headers/ICMole
    source/headers/ICTools
    source/headers/ICPars
    source/headers/ICCalcs
    source/headers/ICRSA
    source/headers/ICCalcs
    source/headers/ICSVM
    source/external/nanoflann
)
add_executable(IChem ${SOURCES})

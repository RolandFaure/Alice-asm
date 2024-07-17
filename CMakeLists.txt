cmake_minimum_required(VERSION 3.8.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(AliceAsm)


# Local header files here ONLY
set(TARGET_H
    assembly.h,
    basic_graph_manipulation.h,
    clipp.h,
    reduce_and_expand.h,
    rolling_hash.h,
    graphunzip.h
   )

# Local source files here
set(TARGET_SRC
    assembly.cpp,
    basic_graph_manipulation.cpp,
    main.cpp,
    reduce_and_expand.cpp,
    rolling_hash.cpp,
    graphunzip.cpp
    )


file(GLOB SOURCE "main.cpp" "assembly.cpp" "reduce_and_expand.cpp" "basic_graph_manipulation.cpp" "rolling_hash.cpp")
add_executable(aliceasm ${SOURCE})
target_compile_options (aliceasm PRIVATE -fopenmp)
target_compile_options (aliceasm PRIVATE -lz)
target_compile_options (aliceasm PRIVATE -g)
target_compile_options (aliceasm PRIVATE -lm)
target_compile_options (aliceasm PRIVATE -std=c++17)
target_compile_options (aliceasm PRIVATE -O3)
target_compile_options (aliceasm PRIVATE -march=x86-64)

# file (GLOB SOURCE_FASTG2GFA "fastg2gfa.c")
# add_executable(fastg2gfa ${SOURCE_FASTG2GFA})
# target_compile_options (fastg2gfa PRIVATE -g)
# target_compile_options (fastg2gfa PRIVATE -O3)
# target_compile_options (fastg2gfa PRIVATE -march=x86-64)
# target_link_libraries(fastg2gfa PRIVATE z)

file(GLOB SOURCE "graphunzip.cpp")
add_executable(graphunzip ${SOURCE})


find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(aliceasm PRIVATE OpenMP::OpenMP_CXX)
endif()




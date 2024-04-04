cmake_minimum_required(VERSION 3.8.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(AliceAsm)
#add_executable(${PROJECT_NAME} "main.cpp")

# Local header files here ONLY
set(TARGET_H
    basic_graph_manipulation.h,
    rolling_hash.h,
    clipp.h
   )

# Local source files here
set(TARGET_SRC
    reduce.cpp,
    basic_graph_manipulation.cpp,
    rolling_hash.cpp
    )

file(GLOB SOURCE "reduce.cpp" "basic_graph_manipulation.cpp" "rolling_hash.cpp")
add_executable(aliceasm ${SOURCE})
# target_include_directories(aliceasm PRIVATE ntHash/build/include)
# target_link_libraries(aliceasm PRIVATE ${CMAKE_SOURCE_DIR}/ntHash/build/lib/libnthash.a)
target_compile_options (aliceasm PRIVATE -fopenmp)
target_compile_options (aliceasm PRIVATE -g)
target_compile_options (aliceasm PRIVATE -lm)
target_compile_options (aliceasm PRIVATE -std=c++17)
target_compile_options (aliceasm PRIVATE -O3)
target_compile_options (aliceasm PRIVATE -march=x86-64)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(aliceasm PRIVATE OpenMP::OpenMP_CXX)
endif()


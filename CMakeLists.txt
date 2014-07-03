cmake_minimum_required (VERSION 2.6)
project (boxnet)
add_subdirectory(src)
add_subdirectory(doc)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99") # always use c99




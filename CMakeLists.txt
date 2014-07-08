cmake_minimum_required (VERSION 2.6)
project (boxnet)

# always use c99, all warnings, strict aliasing
add_definitions(-std=c99 -Wall -fstrict-aliasing)

add_subdirectory(src)
add_subdirectory(doc)


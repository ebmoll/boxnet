add_executable(how_to_use how_to_use.txt.c)
add_test(how_to_use how_to_use)
include_directories ("${PROJECT_SOURCE_DIR}/include")
target_link_libraries (how_to_use boxnet)


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99") # always use c99

add_library(boxnet boxnet.c)
include_directories ("${PROJECT_SOURCE_DIR}/include")


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99") # always use c99

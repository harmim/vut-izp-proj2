cmake_minimum_required(VERSION 3.6)
project(proj2)

set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_C_FLAGS "-std=c99 -Wall -Wextra -Werror -DDEBUG")

set(SOURCE_FILES proj2.c)
add_executable(proj2 ${SOURCE_FILES})

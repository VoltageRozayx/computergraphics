cmake_minimum_required(VERSION 3.29)
project(Lab2++)

set(CMAKE_CXX_STANDARD 23)

add_executable(Lab2++ main.cpp
        lodepng.cpp
        giflib/dgif_lib.c
        giflib/egif_lib.c
        giflib/getarg.c
        giflib/gif_err.c
        giflib/gif_font.c
        giflib/gif_hash.c
        giflib/gifalloc.c
        giflib/openbsd-reallocarray.c
        giflib/qprintf.c
        giflib/quantize.c
)
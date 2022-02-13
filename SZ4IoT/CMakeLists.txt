cmake_minimum_required(VERSION 3.13)
project(SZ_small C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(SZ_small
        adler32.c
        compress.c
        conf.h
        crc32.c
        crc32.h
        deflate.c
        deflate.h
        gzclose.c
        gzguts.h
        gzlib.c
        gzread.c
        gzwrite.c
        infback.c
        inffast.c
        inffast.h
        inffixed.h
        inflate.c
        inflate.h
        inftrees.c
        inftrees.h
        main.c
        sz.h
        trees.c
        trees.h
        uncompr.c
        zconf.h
        zlib.h
        zutil.c
        zutil.h)

target_link_libraries(SZ_small m)
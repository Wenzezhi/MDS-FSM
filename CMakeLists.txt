cmake_minimum_required(VERSION 3.20)

project(MDS_FSM LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(MDS_ENABLE_AVX2 "Enable AVX2 code paths" ON)
option(MDS_ENABLE_POPCNT "Enable POPCNT code paths" ON)
option(MDS_USE_MIMALLOC "Use mimalloc allocator when available" ON)
set(MDS_MIMALLOC_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/mimalloc")

file(GLOB_RECURSE MDS_CPP_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
)

add_executable(mds-fsm ${MDS_CPP_SOURCES})
target_include_directories(mds-fsm PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

if(MDS_USE_MIMALLOC)
    if(NOT EXISTS "${MDS_MIMALLOC_SOURCE_DIR}/CMakeLists.txt")
        message(FATAL_ERROR "Vendored mimalloc source not found at ${MDS_MIMALLOC_SOURCE_DIR}")
    endif()

    set(MI_BUILD_TESTS OFF CACHE BOOL "Disable mimalloc tests" FORCE)
    set(MI_BUILD_SHARED OFF CACHE BOOL "Build mimalloc shared library" FORCE)
    set(MI_BUILD_STATIC ON CACHE BOOL "Build mimalloc static library" FORCE)
    set(MI_BUILD_OBJECT OFF CACHE BOOL "Build mimalloc object library" FORCE)
    set(MI_INSTALL_TOPLEVEL OFF CACHE BOOL "Disable mimalloc install rules" FORCE)
    set(MI_OVERRIDE OFF CACHE BOOL "Disable mimalloc malloc override in the allocator build" FORCE)

    add_subdirectory("${MDS_MIMALLOC_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/mimalloc" EXCLUDE_FROM_ALL)

    if(TARGET mimalloc-static)
        target_compile_options(mimalloc-static PRIVATE
            $<$<NOT:$<C_COMPILER_ID:MSVC>>:-Wno-error=date-time>
        )
        target_link_libraries(mds-fsm PRIVATE mimalloc-static)
    elseif(TARGET mimalloc)
        target_compile_options(mimalloc PRIVATE
            $<$<NOT:$<C_COMPILER_ID:MSVC>>:-Wno-error=date-time>
        )
        target_link_libraries(mds-fsm PRIVATE mimalloc)
    else()
        message(FATAL_ERROR "MDS_USE_MIMALLOC is ON, but no mimalloc target is available")
    endif()

    target_sources(mds-fsm PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/mimalloc_override.cpp")
    target_compile_definitions(mds-fsm PRIVATE MDS_USE_MIMALLOC_RUNTIME=1)
endif()

if(MSVC)
    target_compile_options(mds-fsm PRIVATE
        /EHsc
        /permissive-
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/Oi>
        $<$<CONFIG:Release>:/Ot>
        $<$<CONFIG:Release>:/GL>
        $<$<CONFIG:Release>:/Gy>
        $<$<CONFIG:Release>:/Gw>
        $<$<CONFIG:Debug>:/Od>
        $<$<CONFIG:Debug>:/Zi>
        $<$<CONFIG:Debug>:/RTC1>
    )
    target_link_options(mds-fsm PRIVATE
        $<$<CONFIG:Release>:/LTCG>
        $<$<CONFIG:Release>:/OPT:REF>
        $<$<CONFIG:Release>:/OPT:ICF>
    )
else()
    target_compile_options(mds-fsm PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()

if(MDS_ENABLE_AVX2)
    if(MSVC)
        target_compile_options(mds-fsm PRIVATE /arch:AVX2)
    else()
        target_compile_options(mds-fsm PRIVATE -mavx2)
    endif()
endif()

if(MDS_ENABLE_POPCNT)
    if(NOT MSVC)
        target_compile_options(mds-fsm PRIVATE -mpopcnt)
    endif()
endif()

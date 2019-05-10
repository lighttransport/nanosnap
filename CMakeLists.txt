cmake_minimum_required(VERSION 3.5.1)

# ------------------------------------------------------------------------------
# ----------------------------------- Options ----------------------------------
# ------------------------------------------------------------------------------
option(NANOSNP_BUILD_EXAMPLE "Build NanoSNP example codes" ON)
option(NANOSNP_BUILD_TEST "Build NanoSNP test codes" ON)
message(STATUS "NanoSNP: Build examples: ${BUILD_EXAMPLE}, tests: ${BUILD_TEST}")


# ------------------------------------------------------------------------------
# ----------------------------------- Common -----------------------------------
# ------------------------------------------------------------------------------
project(nanosnp CXX C)
set (CMAKE_CXX_STANDARD 11)  # C++ 11

# Print make commands for debug
# set(CMAKE_VERBOSE_MAKEFILE 1)

# Set default build type
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# cmake modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/sanitizers)
find_package(Sanitizers) # Address sanitizer (-DSANITIZE_ADDRESS=ON)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Warning options
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Weveryting -Wall")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} --pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wundef -Wno-unknown-pragmas")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} /W4")
endif()


# core include
include_directories(${PROJECT_SOURCE_DIR}/include)

# Internal third party includes
# include_directories(
#     ${PROJECT_SOURCE_DIR}/third_party
# )

# ------------------------------- libnanosnp.so -------------------------------
add_library(nanosnp SHARED
    ${PROJECT_SOURCE_DIR}/src/signal.cc
    ${PROJECT_SOURCE_DIR}/src/audio_wav.cc
)
add_sanitizers(nanosnp)

if (WIN32)
    set_target_properties(nanosnp PROPERTIES
        OUTPUT_NAME libnanosnp
    )
endif()

if (NANOSNP_BUILD_EXAMPLE OR NANOSNP_BUILD_TEST)
    add_library(nanosnp_test STATIC
        ${PROJECT_SOURCE_DIR}/tests/main.cc
    )

    target_include_directories(nanosnp_test PRIVATE ${PROJECT_SOURCE_DIR}/third_party/doctest)
    add_sanitizers(nanosnp_test)

    if (WIN32)
        set_target_properties(nanosnp_test PROPERTIES
            OUTPUT_NAME libnanosnp_test
        )
    endif()
endif()

# ------------------------------ test executables ------------------------------
#if (BUILD_TEST)
#    add_executable(test_funcs
#        ${PROJECT_SOURCE_DIR}/tests/test_funcs.cpp
#    )
#    target_link_libraries(test_funcs nanosnp nanosnp_test)
#    target_link_libraries(test_funcs glad glfw ${OPENGL_LIBRARY})
#    add_sanitizers(test_funcs)
#endif()

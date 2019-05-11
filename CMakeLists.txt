cmake_minimum_required(VERSION 3.5.1)

# ------------------------------------------------------------------------------
# ----------------------------------- Options ----------------------------------
# ------------------------------------------------------------------------------
option(NANOSNAP_BUILD_EXAMPLES "Build NanoSNAP example codes" ON)
option(NANOSNAP_BUILD_TESTS "Build NanoSNAP test codes" ON)
message(STATUS "NanoSNAP: Build examples: ${BUILD_EXAMPLES}, tests: ${BUILD_TESTS}")


# ------------------------------------------------------------------------------
# ----------------------------------- Common -----------------------------------
# ------------------------------------------------------------------------------
project(nanosnap CXX C)
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

# ------------------------------- libnanosnap.so -------------------------------
add_library(nanosnap SHARED
    ${PROJECT_SOURCE_DIR}/src/signal.cc
    ${PROJECT_SOURCE_DIR}/src/audio_wav.cc
    ${PROJECT_SOURCE_DIR}/src/speech_features.cc
)
add_sanitizers(nanosnap)

if (WIN32)
    set_target_properties(nanosnap PROPERTIES
        OUTPUT_NAME libnanosnap
    )
endif()

# ------------------------------ test executables ------------------------------
if (NANOSNAP_BUILD_TESTS)
    add_executable(test_nanosnap
        ${PROJECT_SOURCE_DIR}/tests/main.cc
        ${PROJECT_SOURCE_DIR}/tests/test_audio.cc
        ${PROJECT_SOURCE_DIR}/tests/test_speech_features.cc
        ${PROJECT_SOURCE_DIR}/tests/test_medfilt1.cc
    )
    target_include_directories(test_nanosnap PRIVATE ${PROJECT_SOURCE_DIR}/third_party/doctest)
    target_link_libraries(test_nanosnap nanosnap)
    add_sanitizers(test_nanosnap)
endif()

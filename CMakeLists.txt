cmake_minimum_required(VERSION 3.5.1)


# ------------------------------------------------------------------------------
# ----------------------------------- Options ----------------------------------
# ------------------------------------------------------------------------------
option(NANOSNAP_BUILD_EXAMPLES "Build NanoSNAP example codes" ON)
option(NANOSNAP_BUILD_TESTS "Build NanoSNAP test codes" OFF)
message(STATUS "NanoSNAP: Build examples: ${NANOSNAP_BUILD_EXAMPLES}, tests: ${NANOSNAP_BUILD_TESTS}")

if (NANOSNAP_BUILD_TESTS)
  if (NOT IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/third_party/doctest/doctest")
    message(FATAL_ERROR "doctest submodule directory is missing!\n"
      "Please run \"git submodule update --init --recursive --depth 1\" to checkout submodules")
  endif ()
endif (NANOSNAP_BUILD_TESTS)

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

# core include
include_directories(${PROJECT_SOURCE_DIR}/include)

# Internal third party includes
# include_directories(
#     ${PROJECT_SOURCE_DIR}/third_party
# )

# ------------------------------- libnanosnap -------------------------------
add_library(nanosnap STATIC
    ${PROJECT_SOURCE_DIR}/include/nanosnap/nanosnap.h
    ${PROJECT_SOURCE_DIR}/src/ndarray.cc
    ${PROJECT_SOURCE_DIR}/src/signal.cc
    ${PROJECT_SOURCE_DIR}/src/audio_wav.cc
    ${PROJECT_SOURCE_DIR}/src/speech_features.cc
    ${PROJECT_SOURCE_DIR}/src/fft.cc
    ${PROJECT_SOURCE_DIR}/src/rng.cc
    ${PROJECT_SOURCE_DIR}/src/text-load-save.cc
    ${PROJECT_SOURCE_DIR}/src/image.cc
    ${PROJECT_SOURCE_DIR}/src/mel.cc
    ${PROJECT_SOURCE_DIR}/src/pocketfft.c
)
add_sanitizers(nanosnap)

# raise compiler warning level.
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  target_compile_options(nanosnap PRIVATE
     $<$<COMPILE_LANGUAGE:CXX>: -Weverything -Werror -Wno-c++98-compat -Wno-padded>
     )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GCC")
  target_compile_options(nanosnap PRIVATE
     $<$<COMPILE_LANGUAGE:CXX>: --pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wundef -Wno-unknown-pragmas>
     )
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  target_compile_options(nanosnap PRIVATE
     $<$<COMPILE_LANGUAGE:CXX>: /W4>
     )
endif()


#if (WIN32)
#    set_target_properties(nanosnap PROPERTIES
#        OUTPUT_NAME libnanosnap
#    )
#endif()

# ------------------------------ test executables ------------------------------
if (NANOSNAP_BUILD_TESTS)
    add_executable(test_nanosnap
        ${PROJECT_SOURCE_DIR}/tests/main.cc
        ${PROJECT_SOURCE_DIR}/tests/test_audio.cc
        ${PROJECT_SOURCE_DIR}/tests/test_speech_features.cc
        ${PROJECT_SOURCE_DIR}/tests/test_medfilt1.cc
        ${PROJECT_SOURCE_DIR}/tests/test_rfft.cc
        ${PROJECT_SOURCE_DIR}/tests/test_ifft.cc
        ${PROJECT_SOURCE_DIR}/tests/test_convolve.cc
        ${PROJECT_SOURCE_DIR}/tests/test_librosa_stft.cc
        ${PROJECT_SOURCE_DIR}/tests/test_librosa_istft.cc
        ${PROJECT_SOURCE_DIR}/tests/test_librosa_filters_mel.cc
        ${PROJECT_SOURCE_DIR}/tests/test_signal_get_window_hann.cc
        # it looks numpy.random uses slightly different configuration of Mersenne Twistter
        # and does not produce same result with std::mt19937. Disable random test for a while.
        #${PROJECT_SOURCE_DIR}/tests/test_random_uniform.cc
    )
    target_include_directories(test_nanosnap PRIVATE ${PROJECT_SOURCE_DIR}/src)
    target_include_directories(test_nanosnap PRIVATE ${PROJECT_SOURCE_DIR}/third_party/doctest)
    target_link_libraries(test_nanosnap nanosnap)
    add_sanitizers(test_nanosnap)

    # [VisualStudio]
    if (WIN32)
	# Set `test_nanosnap` as a startup project for VS IDE.
        if (CMAKE_VERSION VERSION_GREATER 3.6.0)
	    set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT test_nanosnap)
	endif ()

	# For easier debugging in VS IDE.
	if (CMAKE_VERSION VERSION_GREATER 3.8.0)
            set_target_properties(test_nanosnap PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
	endif ()

    endif()
else()
    # [VisualStudio]
    if (WIN32)
	# Set `nanosnap` as a startup project for VS IDE.
        if (CMAKE_VERSION VERSION_GREATER 3.6.0)
	    set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT nanosnap)
	endif ()
    endif()
endif()


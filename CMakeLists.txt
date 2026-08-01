cmake_minimum_required(VERSION 3.22)
project(NoiseRepellent VERSION 3.0.0 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable Position Independent Code (-fPIC) globally for all static libs linked into .so plugins
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Enable MSVC C11 Atomics support globally for Windows MSVC builds
if(MSVC)
    add_compile_options($<$<COMPILE_LANGUAGE:C>:/std:c11>)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/std:c++17>)
    # Required for MSVC to enable stdatomic.h support in MSVC C11 mode
    add_compile_options(/experimental:c11atomics) 
endif()

include(FetchContent)

# Fetch JUCE framework
message(STATUS "Fetching JUCE framework...")
FetchContent_Declare(
    JUCE
    GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
    GIT_TAG        8.0.4
)
FetchContent_MakeAvailable(JUCE)

# Add JUCE Audio Plugin target
juce_add_plugin(NoiseRepellent
    COMPANY_NAME "Luciano Dato"
    BUNDLE_ID "com.lucianodato.NoiseRepellent"
    LV2URI "https://github.com/lucianodato/noise-repellent"
    IS_SYNTH FALSE
    NEEDS_MIDI_INPUT FALSE
    NEEDS_MIDI_OUTPUT FALSE
    IS_MIDI_EFFECT FALSE
    EDITOR_WANTS_KEYBOARD_FOCUS FALSE
    COPY_PLUGIN_AFTER_BUILD FALSE
    PLUGIN_MANUFACTURER_CODE LdNR
    PLUGIN_CODE Nrep
    FORMATS VST3 AU LV2
    PRODUCT_NAME "Noise Repellent"
)

if(UNIX AND NOT APPLE)
    # Statically link GCC and C++ standard runtimes into the plugin library
    target_link_options(NoiseRepellent PRIVATE 
        -static-libgcc 
        -static-libstdc++
    )

    # --- Static FreeType for Linux ---
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(FT_DISABLE_ZLIB ON CACHE BOOL "" FORCE)
    set(FT_DISABLE_PNG ON CACHE BOOL "" FORCE)
    set(FT_DISABLE_HARFBUZZ ON CACHE BOOL "" FORCE)
    set(FT_DISABLE_BROTLI ON CACHE BOOL "" FORCE)

    FetchContent_Declare(
        freetype
        GIT_REPOSITORY https://gitlab.freedesktop.org/freetype/freetype.git
        GIT_TAG        VER-2-12-1
    )
    FetchContent_MakeAvailable(freetype)
endif()

if(APPLE)
    # Target modern C++ runtime provided by macOS (10.13+)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.13" CACHE STRING "Minimum macOS version")
endif()

# Include directories
target_include_directories(NoiseRepellent PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Source
)

# Source files
target_sources(NoiseRepellent PRIVATE
    Source/PluginProcessor.cpp
    Source/PluginProcessor.h
    Source/PluginEditor.cpp
    Source/PluginEditor.h
    Source/GUI/LookAndFeel.cpp
    Source/GUI/LookAndFeel.h
    Source/GUI/SpectralVisualizer.cpp
    Source/GUI/SpectralVisualizer.h
)

# Compile definitions
target_compile_definitions(NoiseRepellent PRIVATE
    JUCE_WEB_BROWSER=0
    JUCE_USE_CURL=0
    JUCE_VST3_CAN_REPLACE_VST2=0
)

# Global options for subprojects
set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "" FORCE)

# --- Static Universal FFTW3 for Linux and macOS ---
if(UNIX)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(ENABLE_FLOAT ON CACHE BOOL "" FORCE)     # Single-precision fftw3f
    set(ENABLE_THREADS ON CACHE BOOL "" FORCE)   # Thread support

    set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

    FetchContent_Declare(
        fftw3
        URL https://www.fftw.org/fftw-3.3.10.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(fftw3)

    set(FFTW3F_LIB fftw3f CACHE STRING "Path to fftw3f library" FORCE)
    set(FFTW3_LIB fftw3f CACHE STRING "Path to fftw3 library" FORCE)

    # Make fftw3.h available to subprojects
    include_directories(
        ${fftw3_SOURCE_DIR}/api
        ${fftw3_BINARY_DIR}
    )
endif()

if(APPLE)
    set(ENABLE_AVX OFF CACHE BOOL "" FORCE)
    set(USE_AVX OFF CACHE BOOL "" FORCE)
endif()

# --- Libspecbleach Source Resolution ---
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../libspecbleach/CMakeLists.txt")
    message(STATUS "Using local libspecbleach from ../libspecbleach")
    set(SPECBLEACH_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../libspecbleach")
    add_subdirectory("${SPECBLEACH_SRC}" "${CMAKE_CURRENT_BINARY_DIR}/libspecbleach_build")
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/subprojects/libspecbleach/CMakeLists.txt")
    message(STATUS "Using local libspecbleach from subprojects/")
    set(SPECBLEACH_SRC "${CMAKE_CURRENT_SOURCE_DIR}/subprojects/libspecbleach")
    add_subdirectory("${SPECBLEACH_SRC}" "${CMAKE_CURRENT_BINARY_DIR}/libspecbleach_build")
else()
    message(STATUS "Fetching libspecbleach v0.3.0 remotely...")
    FetchContent_Declare(
        libspecbleach
        GIT_REPOSITORY https://github.com/lucianodato/libspecbleach.git
        GIT_TAG        v0.3.0
    )
    FetchContent_MakeAvailable(libspecbleach)
    FetchContent_GetProperties(libspecbleach SOURCE_DIR SPECBLEACH_SRC)
endif()

# Explicitly add libspecbleach headers directory
target_include_directories(NoiseRepellent PRIVATE
    "${SPECBLEACH_SRC}/src"
)

# Fix libspecbleach AVX symbol issue on x86_64 architecture (Linux & macOS)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|amd64|AMD64)" OR "x86_64" IN_LIST CMAKE_OSX_ARCHITECTURES)
    target_sources(NoiseRepellent PRIVATE
        "${SPECBLEACH_SRC}/src/shared/denoiser_logic/processing/nlm_filter_avx.c"
    )
    set_source_files_properties(
        "${SPECBLEACH_SRC}/src/shared/denoiser_logic/processing/nlm_filter_avx.c"
        PROPERTIES COMPILE_FLAGS "-mavx -mavx2"
    )
endif()

# Hide internal static library symbols from DAWs
set_target_properties(NoiseRepellent PROPERTIES
    C_VISIBILITY_PRESET hidden
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
)

# Link JUCE Modules and DSP Dependencies
target_link_libraries(NoiseRepellent PRIVATE
    libspecbleach
    juce::juce_audio_utils
    juce::juce_audio_processors
    juce::juce_dsp
    juce::juce_gui_extra
    juce::juce_gui_basics
    juce::juce_graphics
    juce::juce_core
)

if(UNIX AND NOT APPLE)
    target_link_libraries(NoiseRepellent PRIVATE freetype)
endif()

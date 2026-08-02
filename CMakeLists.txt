cmake_minimum_required(VERSION 3.22)
project(NoiseRepellent VERSION 0.3.0 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable Position Independent Code (-fPIC) globally for static libs linked into shared object plugins
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Statically link MSVC C/C++ Runtime (/MT) on Windows
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

if(UNIX AND NOT APPLE)
    # Statically link GCC and C++ standard runtimes for ALL CMake targets (including JUCE VST3/LV2 targets)
    add_link_options(-static-libgcc -static-libstdc++)
endif()

# Enable MSVC C11 Atomics support globally for Windows MSVC builds
if(MSVC)
    add_compile_options($<$<COMPILE_LANGUAGE:C>:/std:c11>)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/std:c++17>)
    add_compile_options(/experimental:c11atomics) 
endif()

include(FetchContent)

# --- Static FreeType & Disable Fontconfig for Linux (MUST BE BEFORE JUCE) ---
if(UNIX AND NOT APPLE)
    # Disable JUCE Fontconfig search to eliminate libfontconfig.so
    add_compile_definitions(JUCE_USE_FONTCONFIG=0)
    set(juce_fontconfig_FOUND FALSE CACHE BOOL "" FORCE)

    # Disable external compression/font features in static FreeType
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

# Fetch JUCE framework
message(STATUS "Fetching JUCE framework...")
FetchContent_Declare(
    JUCE
    GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
    GIT_TAG        8.0.4
)
FetchContent_MakeAvailable(JUCE)

# --- Override JUCE's FreeType target to force our static library ---
if(UNIX AND NOT APPLE)
    if(TARGET PkgConfig::juce_freetype)
        set_target_properties(PkgConfig::juce_freetype PROPERTIES
            INTERFACE_LINK_LIBRARIES "freetype"
            INTERFACE_INCLUDE_DIRECTORIES "$<TARGET_PROPERTY:freetype,INTERFACE_INCLUDE_DIRECTORIES>"
        )
    else()
        add_library(PkgConfig::juce_freetype INTERFACE IMPORTED GLOBAL)
        set_target_properties(PkgConfig::juce_freetype PROPERTIES
            INTERFACE_LINK_LIBRARIES "freetype"
            INTERFACE_INCLUDE_DIRECTORIES "$<TARGET_PROPERTY:freetype,INTERFACE_INCLUDE_DIRECTORIES>"
        )
    endif()
endif()

# Add JUCE Audio Plugin target
juce_add_plugin(NoiseRepellent
    COMPANY_NAME "Luciano Dato"
    BUNDLE_ID "com.lucianodato.NoiseRepellent"
    LV2URI "https://github.com/lucianodato/noise-repellent#v0"
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

# --- Static FFTW3 for ALL Platforms (Linux, macOS, Windows) ---
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
    set(SPECBLEACH_SRC "${libspecbleach_SOURCE_DIR}")
endif()

# Explicitly add libspecbleach headers directory to NoiseRepellent
target_include_directories(NoiseRepellent PRIVATE
    "${SPECBLEACH_SRC}/src"
)

# Fix libspecbleach AVX missing symbol on x86_64 architecture
if(TARGET libspecbleach)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|amd64|AMD64)" OR "x86_64" IN_LIST CMAKE_OSX_ARCHITECTURES)
        target_sources(libspecbleach PRIVATE
            "${SPECBLEACH_SRC}/src/shared/denoiser_logic/processing/nlm_filter_avx.c"
        )
        set_source_files_properties(
            "${SPECBLEACH_SRC}/src/shared/denoiser_logic/processing/nlm_filter_avx.c"
            PROPERTIES COMPILE_FLAGS "-mavx -mavx2"
        )
    endif()
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

cmake_minimum_required(VERSION 3.22)
project(NoiseRepellent VERSION 0.3.1 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable Position Independent Code (-fPIC) globally for static libs linked into shared object plugins
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Statically link MSVC C/C++ Runtime (/MT) on Windows
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

if(UNIX AND NOT APPLE)
    # Statically link GCC runtimes AND hide all symbols inside static archives (-Wl,--exclude-libs,ALL)
    # This prevents host DAWs (like Ardour) from experiencing FreeType/FFTW dynamic linker symbol collisions
    add_link_options(-static-libgcc -static-libstdc++ -Wl,--exclude-libs,ALL)
endif()

if(MSVC)
    add_compile_options($<$<COMPILE_LANGUAGE:C>:/std:c11>)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/std:c++17>)
    add_compile_options(/experimental:c11atomics) 
endif()

# ==============================================================================
# Build Options (Default OFF for self-contained, portable release binaries)
# Downstream distro packagers can pass -DUSE_SYSTEM_*=ON to use system libraries
# ==============================================================================
option(USE_SYSTEM_FFTW       "Use system-installed FFTW3 library"       OFF)
option(USE_SYSTEM_FREETYPE   "Use system-installed FreeType library"   OFF)
option(USE_SYSTEM_SPECBLEACH "Use system-installed libspecbleach"      OFF)

include(FetchContent)

# ==============================================================================
# 1. FreeType Handling (Linux)
# ==============================================================================
if(UNIX AND NOT APPLE)
    if(USE_SYSTEM_FREETYPE)
        message(STATUS "Using system-installed FreeType library")
        find_package(Freetype REQUIRED)
    else()
        message(STATUS "Building static FreeType via FetchContent")
        set(JUCE_USE_FONTCONFIG OFF CACHE BOOL "" FORCE)
        add_compile_definitions(JUCE_USE_FONTCONFIG=0)

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
endif()

# ==============================================================================
# 2. JUCE Framework
# ==============================================================================
message(STATUS "Fetching JUCE framework...")
FetchContent_Declare(
    JUCE
    GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
    GIT_TAG        8.0.4
)
FetchContent_MakeAvailable(JUCE)

if(UNIX AND NOT APPLE AND NOT USE_SYSTEM_FREETYPE)
    # Redirect JUCE PkgConfig target to point to our static freetype library
    if(TARGET PkgConfig::juce_freetype)
        set_target_properties(PkgConfig::juce_freetype PROPERTIES
            IMPORTED_LOCATION ""
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

    if(TARGET PkgConfig::juce_fontconfig)
        set_target_properties(PkgConfig::juce_fontconfig PROPERTIES
            IMPORTED_LOCATION ""
            INTERFACE_LINK_LIBRARIES ""
            INTERFACE_INCLUDE_DIRECTORIES ""
        )
    endif()
endif()

# ==============================================================================
# 3. Audio Plugin Target Definition
# ==============================================================================
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
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.13" CACHE STRING "Minimum macOS version")
endif()

target_include_directories(NoiseRepellent PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Source
)

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

target_compile_definitions(NoiseRepellent PRIVATE
    JUCE_WEB_BROWSER=0
    JUCE_USE_CURL=0
    JUCE_VST3_CAN_REPLACE_VST2=0
)

# Global subproject options
set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "" FORCE)

# ==============================================================================
# 4. FFTW3 Dependency
# ==============================================================================
if(USE_SYSTEM_FFTW)
    message(STATUS "Using system-installed FFTW3 library")
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(FFTW3F REQUIRED IMPORTED_TARGET fftw3f)
    set(FFTW_TARGET PkgConfig::FFTW3F)
else()
    message(STATUS "Building static FFTW3 via FetchContent")
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(ENABLE_FLOAT ON CACHE BOOL "" FORCE)
    set(ENABLE_THREADS ON CACHE BOOL "" FORCE)

    set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

    FetchContent_Declare(
        fftw3
        URL https://www.fftw.org/fftw-3.3.10.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(fftw3)

    set(FFTW3F_LIB fftw3f CACHE STRING "Path to fftw3f library" FORCE)
    set(FFTW3_LIB fftw3f CACHE STRING "Path to fftw3 library" FORCE)

    include_directories(
        ${fftw3_SOURCE_DIR}/api
        ${fftw3_BINARY_DIR}
    )
    set(FFTW_TARGET fftw3f)
endif()

if(APPLE)
    set(ENABLE_AVX OFF CACHE BOOL "" FORCE)
    set(USE_AVX OFF CACHE BOOL "" FORCE)
endif()

# ==============================================================================
# 5. libspecbleach Dependency
# ==============================================================================
if(USE_SYSTEM_SPECBLEACH)
    message(STATUS "Using system-installed libspecbleach library")
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(SPECBLEACH REQUIRED IMPORTED_TARGET specbleach)
    set(SPECBLEACH_TARGET PkgConfig::SPECBLEACH)
else()
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

    target_include_directories(NoiseRepellent PRIVATE
        "${SPECBLEACH_SRC}/src"
    )

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
    set(SPECBLEACH_TARGET libspecbleach)
endif()

# ==============================================================================
# 6. Target Linkage & Symbol Visibility Properties
# ==============================================================================
set_target_properties(NoiseRepellent PROPERTIES
    C_VISIBILITY_PRESET hidden
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
)

target_link_libraries(NoiseRepellent PRIVATE
    ${SPECBLEACH_TARGET}
    ${FFTW_TARGET}
    juce::juce_audio_utils
    juce::juce_audio_processors
    juce::juce_dsp
    juce::juce_gui_extra
    juce::juce_gui_basics
    juce::juce_graphics
    juce::juce_core
)

if(UNIX AND NOT APPLE AND USE_SYSTEM_FREETYPE)
    target_link_libraries(NoiseRepellent PRIVATE Freetype::Freetype)
endif()

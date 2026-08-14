cmake_minimum_required(VERSION 3.22)
project(NoiseRepellent VERSION 0.3.2 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable Position Independent Code (-fPIC) globally for static libs linked into shared object plugins
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Statically link MSVC C/C++ Runtime (/MT) on Windows
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

# ==============================================================================
# Build Options (Default OFF for self-contained, portable release binaries)
# Downstream distro packagers can pass -DUSE_SYSTEM_*=ON to use system libraries
# ==============================================================================
option(USE_SYSTEM_FFTW       "Use system-installed FFTW3 library"       OFF)
option(USE_SYSTEM_FREETYPE   "Use system-installed FreeType library"   OFF)
option(USE_SYSTEM_SPECBLEACH "Use system-installed libspecbleach"      OFF)
option(USE_SYSTEM_JUCE       "Use system-installed JUCE library"       OFF)

if(UNIX AND NOT APPLE)
    # Statically link GCC runtimes AND hide all symbols inside static archives (-Wl,--exclude-libs,ALL)
    # This prevents host DAWs (like Ardour) from experiencing FreeType/FFTW dynamic linker symbol collisions
    # Only apply for portable builds (when NOT using system dependencies)
    if(NOT USE_SYSTEM_FFTW AND NOT USE_SYSTEM_FREETYPE AND NOT USE_SYSTEM_SPECBLEACH AND NOT USE_SYSTEM_JUCE)
        add_link_options(-static-libgcc -static-libstdc++ -Wl,--exclude-libs,ALL)
    endif()
endif()

if(MSVC)
    add_compile_options($<$<COMPILE_LANGUAGE:C>:/std:c11>)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:/std:c++17>)
    add_compile_options(/experimental:c11atomics)
    add_compile_options(/openmp:llvm)
endif()

include(FetchContent)

# ==============================================================================
# 1. FreeType Handling (Portable Static Build via FetchContent)
# ==============================================================================
if(UNIX AND NOT APPLE AND NOT USE_SYSTEM_FREETYPE AND NOT USE_SYSTEM_JUCE)
    message(STATUS "Building static FreeType 2.13.2 via FetchContent...")
    
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(FT_DISABLE_ZLIB ON CACHE BOOL "" FORCE)
    set(FT_DISABLE_PNG ON CACHE BOOL "" FORCE)
    set(FT_DISABLE_HARFBUZZ ON CACHE BOOL "" FORCE)
    set(FT_DISABLE_BROTLI ON CACHE BOOL "" FORCE)
    set(FT_DISABLE_BZIP2 ON CACHE BOOL "" FORCE)

    # Disable Fontconfig in JUCE C++
    set(JUCE_USE_FONTCONFIG OFF CACHE BOOL "" FORCE)
    add_compile_definitions(JUCE_USE_FONTCONFIG=0)

    FetchContent_Declare(
        freetype
        GIT_REPOSITORY https://gitlab.freedesktop.org/freetype/freetype.git
        GIT_TAG        VER-2-13-2
    )
    FetchContent_MakeAvailable(freetype)

    # --------------------------------------------------------------------------
    # THE PKG-CONFIG INTERCEPTOR
    # JUCE strictly demands pkg-config for freetype2. We trick it with a proxy 
    # script that returns empty flags so system libfreetype.so is never linked.
    # --------------------------------------------------------------------------
    find_program(REAL_PKG_CONFIG pkg-config)
    if(REAL_PKG_CONFIG)
        set(FAKE_PKG_CONFIG "${CMAKE_CURRENT_BINARY_DIR}/fake-pkg-config.sh")
        file(WRITE "${FAKE_PKG_CONFIG}"
        "#!/bin/sh
        IS_TARGET=0
        for arg in \"$@\"; do
            if [ \"$arg\" = \"freetype2\" ] || [ \"$arg\" = \"fontconfig\" ]; then
                IS_TARGET=1
            fi
        done

        if [ \"$IS_TARGET\" = \"1\" ]; then
            for arg in \"$@\"; do
                if [ \"$arg\" = \"--modversion\" ]; then
                    echo \"2.13.2\"
                    exit 0
                fi
            done
            # Return 0 (success) but output nothing for --cflags and --libs
            exit 0
        fi

        # Pass all other requests (alsa, webkit2gtk, etc.) to the real pkg-config
        exec \"${REAL_PKG_CONFIG}\" \"$@\"
        ")
        execute_process(COMMAND chmod +x "${FAKE_PKG_CONFIG}")
        # Force CMake to use our interceptor script instead of the real pkg-config
        set(PKG_CONFIG_EXECUTABLE "${FAKE_PKG_CONFIG}" CACHE FILEPATH "Fake pkg-config" FORCE)
    endif()
endif()

# ==============================================================================
# 2. JUCE Framework
# ==============================================================================
if(USE_SYSTEM_JUCE)
    message(STATUS "Using system-installed JUCE library...")
    find_package(JUCE CONFIG REQUIRED)
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../JUCE/CMakeLists.txt")
    message(STATUS "Using local JUCE from ../JUCE")
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../JUCE" "${CMAKE_CURRENT_BINARY_DIR}/juce_build")
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/subprojects/JUCE/CMakeLists.txt")
    message(STATUS "Using local JUCE from subprojects/JUCE")
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/subprojects/JUCE" "${CMAKE_CURRENT_BINARY_DIR}/juce_build")
else()
    message(STATUS "Fetching JUCE framework...")
    FetchContent_Declare(
        JUCE
        GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
        GIT_TAG        8.0.4
    )
    FetchContent_MakeAvailable(JUCE)
endif()

# ==============================================================================
# 3. Bind Static FreeType to JUCE
# ==============================================================================
if(UNIX AND NOT APPLE)
    if(USE_SYSTEM_FREETYPE)
        # Clear any previously-cached fake pkg-config interceptor to allow real pkg-config discovery
        unset(PKG_CONFIG_EXECUTABLE CACHE)
        # Discover system FreeType package
        find_package(Freetype REQUIRED)
    elseif(NOT USE_SYSTEM_JUCE)
        # The interceptor gave JUCE an empty PkgConfig target.
        # Now we inject our static freetype build directly into JUCE's graphics module.
        if(TARGET juce_graphics)
            target_link_libraries(juce_graphics INTERFACE freetype)
        endif()
    endif()
endif()

# ==============================================================================
# 4. Audio Plugin Target Definition
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
# 5. FFTW3 Dependency
# ==============================================================================
# Disable AVX on Apple platforms BEFORE fetching/configuring FFTW3
if(APPLE)
    set(ENABLE_AVX OFF CACHE BOOL "" FORCE)
    set(USE_AVX OFF CACHE BOOL "" FORCE)
endif()

if(USE_SYSTEM_FFTW)
    message(STATUS "Using system-installed FFTW3 library")
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(FFTW3F REQUIRED IMPORTED_TARGET fftw3f)
    set(FFTW_TARGET PkgConfig::FFTW3F)
else()
    message(STATUS "Building static FFTW3 via FetchContent...")
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

# ==============================================================================
# 6. libspecbleach Dependency
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
        message(STATUS "Fetching libspecbleach v0.3.1 remotely...")
        FetchContent_Declare(
            libspecbleach
            GIT_REPOSITORY https://github.com/lucianodato/libspecbleach.git
            GIT_TAG        v0.3.1
        )
        FetchContent_MakeAvailable(libspecbleach)
        set(SPECBLEACH_SRC "${libspecbleach_SOURCE_DIR}")
    endif()

    target_include_directories(NoiseRepellent PRIVATE
        "${SPECBLEACH_SRC}/src"
    )

    if(TARGET libspecbleach)
        set(AVX_C_FILE "${SPECBLEACH_SRC}/src/shared/denoiser_logic/processing/nlm_filter_avx.c")
        if(EXISTS "${AVX_C_FILE}")
            if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|amd64|AMD64|x86-64)" OR "x86_64" IN_LIST CMAKE_OSX_ARCHITECTURES)
                target_sources(libspecbleach PRIVATE "${AVX_C_FILE}")
                # Apply compiler-specific AVX flags
                if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
                    set_source_files_properties("${AVX_C_FILE}" PROPERTIES COMPILE_FLAGS "-mavx -mavx2")
                elseif(MSVC)
                    set_source_files_properties("${AVX_C_FILE}" PROPERTIES COMPILE_FLAGS "/arch:AVX2")
                endif()
            endif()
        endif()
    endif()
    set(SPECBLEACH_TARGET libspecbleach)
endif()

# ==============================================================================
# 7. Target Linkage & Symbol Visibility Properties
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

if(UNIX AND NOT APPLE)
    if(USE_SYSTEM_FREETYPE)
        target_link_libraries(NoiseRepellent PRIVATE Freetype::Freetype)
    elseif(NOT USE_SYSTEM_JUCE)
        target_link_libraries(NoiseRepellent PRIVATE freetype)
    endif()
endif()

# ==============================================================================
# 8. Installation Rules
# ==============================================================================
include(GNUInstallDirs)

if(APPLE)
    set(DEFAULT_VST3_DIR "/Library/Audio/Plug-Ins/VST3")
    set(DEFAULT_AU_DIR   "/Library/Audio/Plug-Ins/Components")
    set(DEFAULT_LV2_DIR  "/Library/Audio/Plug-Ins/LV2")
elseif(WIN32)
    if(DEFINED ENV{CommonProgramFiles})
        file(TO_CMAKE_PATH "$ENV{CommonProgramFiles}" COMMON_PROGRAM_FILES)
        set(DEFAULT_VST3_DIR "${COMMON_PROGRAM_FILES}/VST3")
        set(DEFAULT_LV2_DIR  "${COMMON_PROGRAM_FILES}/LV2")
    else()
        set(DEFAULT_VST3_DIR "${CMAKE_INSTALL_PREFIX}/vst3")
        set(DEFAULT_LV2_DIR  "${CMAKE_INSTALL_PREFIX}/lv2")
    endif()
else() # Linux / Unix
    set(DEFAULT_VST3_DIR "${CMAKE_INSTALL_LIBDIR}/vst3")
    set(DEFAULT_LV2_DIR  "${CMAKE_INSTALL_LIBDIR}/lv2")
endif()


set(INSTALL_VST3_DIR "${DEFAULT_VST3_DIR}" CACHE PATH "Installation directory for VST3 plugins")
set(INSTALL_AU_DIR   "${DEFAULT_AU_DIR}"   CACHE PATH "Installation directory for AU plugins")
set(INSTALL_LV2_DIR  "${DEFAULT_LV2_DIR}"  CACHE PATH "Installation directory for LV2 plugins")

if(TARGET NoiseRepellent_VST3)
    install(DIRECTORY "$<TARGET_PROPERTY:NoiseRepellent_VST3,JUCE_PLUGIN_ARTEFACT_FILE>"
        DESTINATION "${INSTALL_VST3_DIR}"
        USE_SOURCE_PERMISSIONS
        COMPONENT VST3
    )
endif()

if(TARGET NoiseRepellent_AU)
    install(DIRECTORY "$<TARGET_PROPERTY:NoiseRepellent_AU,JUCE_PLUGIN_ARTEFACT_FILE>"
        DESTINATION "${INSTALL_AU_DIR}"
        USE_SOURCE_PERMISSIONS
        COMPONENT AU
    )
endif()

if(TARGET NoiseRepellent_LV2)
    install(DIRECTORY "$<TARGET_PROPERTY:NoiseRepellent_LV2,JUCE_PLUGIN_ARTEFACT_FILE>"
        DESTINATION "${INSTALL_LV2_DIR}"
        USE_SOURCE_PERMISSIONS
        COMPONENT LV2
    )
endif()



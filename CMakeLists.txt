cmake_minimum_required(VERSION 3.22)
project(NoiseRepellent VERSION 3.0.0 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
    LV2URI "http://lucianodato.net/plugins/nrepellent"
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

# Add libspecbleach static library dependency via CMake
set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../libspecbleach/CMakeLists.txt")
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../libspecbleach" "${CMAKE_CURRENT_BINARY_DIR}/libspecbleach_build")
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/subprojects/libspecbleach/CMakeLists.txt")
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/subprojects/libspecbleach" "${CMAKE_CURRENT_BINARY_DIR}/libspecbleach_build")
else()
    FetchContent_Declare(
        libspecbleach
        GIT_REPOSITORY https://github.com/lucianodato/libspecbleach.git
        GIT_TAG        v0.3.0
    )
    FetchContent_MakeAvailable(libspecbleach)
endif()

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

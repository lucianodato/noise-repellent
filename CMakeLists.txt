cmake_minimum_required(VERSION 3.22)
project(NoiseRepellent VERSION 3.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
    FORMATS VST3 AU LV2 Standalone
    PRODUCT_NAME "Noise Repellent"
)

# Include directories
target_include_directories(NoiseRepellent PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/Source
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/subprojects/libspecbleach/include
    ${CMAKE_CURRENT_SOURCE_DIR}/../libspecbleach/include
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
    src/signal_crossfade.c
    src/signal_crossfade.h
    src/noise_profile_state.c
    src/noise_profile_state.h
)

# Compile definitions
target_compile_definitions(NoiseRepellent PRIVATE
    JUCE_WEB_BROWSER=0
    JUCE_USE_CURL=0
    JUCE_VST3_CAN_REPLACE_VST2=0
    JUCE_DISPLAY_SPLASH_SCREEN=0
)

find_library(FFTW3F_LIB fftw3f REQUIRED)
find_library(MATH_LIB m)
find_library(OMP_LIB omp HINTS /opt/homebrew/opt/libomp/lib /usr/local/lib)

# Link JUCE Modules and DSP Dependencies
target_link_libraries(NoiseRepellent PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../libspecbleach/build/src/libspecbleach.a
    ${FFTW3F_LIB}
    ${MATH_LIB}
    ${OMP_LIB}
    juce::juce_audio_utils
    juce::juce_audio_processors
    juce::juce_gui_extra
    juce::juce_gui_basics
    juce::juce_graphics
    juce::juce_core
)

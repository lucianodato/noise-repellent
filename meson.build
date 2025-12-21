project('nrepellent.lv2', 'c',
        version: '0.2.4',
        default_options: ['default_library=shared', 'c_std=c17', 'warning_level=2'])

# Note: Meson 0.60.0 or newer is recommended

# Get the host operating system and cpu architecture
current_os = host_machine.system()
current_arch = build_machine.cpu_family()
cc = meson.get_compiler('c')

# Install folder
lv2_directory = get_option('lv2dir')
if lv2_directory == ''
  if current_os == 'darwin'
    lv2_directory = '/Library/Audio/Plug-Ins/LV2'
  elif current_os == 'windows'
    lv2_directory = join_paths(get_option('prefix'), 'bin', 'lv2')
  else
    lv2_directory = join_paths(get_option('libdir'), 'lv2')
  endif
endif
install_folder = join_paths(lv2_directory, meson.project_name())

# Sources to compile
common_src = ['src/signal_crossfade.c']
noise_repellent_src = ['plugins/nrepellent.c', 'src/noise_profile_state.c']
noise_repellent_adaptive_src = 'plugins/nrepellent-adaptive.c'

# Dependencies for noise repellent
lv2_dep = dependency('lv2', required: true)
libspecbleach_dep = dependency('libspecbleach',
                               fallback: ['libspecbleach', 'libspecbleach_dep'],
                               default_options: ['default_library=static'],
                               required: true)
m_dep = cc.find_library('m', required: true)
all_dep = [lv2_dep, libspecbleach_dep, m_dep]

# Shared c_args for libraries
lib_c_args = ['-fvisibility=hidden']

# Compiler warning flags
custom_warning_level = get_option('custom_warning_level')
if custom_warning_level >= 1
  lib_c_args += cc.get_supported_arguments([
    '-Wall',
    '-Wextra',
  ])
endif

if custom_warning_level >= 2
  lib_c_args += cc.get_supported_arguments([
    '-Wpedantic',
    '-Wno-unused-parameter',
    '-Wno-unused-variable',
  ])
endif

# Treat warnings as errors if requested
if get_option('treat_warnings_as_errors')
  lib_c_args += ['-Werror']
endif

# Architecture-specific optimizations
# x86 and x86_64 optimizations (excluding macOS)
if (current_arch == 'x86' or current_arch == 'x86_64') and current_os != 'darwin'
  lib_c_args += cc.get_supported_arguments([
    '-msse',
    '-msse2',
    '-mfpmath=sse',
    '-ffast-math',
    '-fomit-frame-pointer',
    '-fno-finite-math-only'
  ])
endif

# ARM64 optimizations (NEON)
if current_arch == 'aarch64'
  lib_c_args += cc.get_supported_arguments([
    '-mfpu=neon-fp-armv8',
  ])
endif

# Sanitizers (for debug builds)
if get_option('enable_sanitizers')
  if get_option('buildtype') == 'debug'
    sanitizer_args = []
    if get_option('sanitize_address')
      sanitizer_args += ['-fsanitize=address']
    endif
    if get_option('sanitize_undefined')
      sanitizer_args += ['-fsanitize=undefined']
    endif
    lib_c_args += cc.get_supported_arguments(sanitizer_args)
  else
    warning('Sanitizers are only enabled in debug builds')
  endif
endif

# Get configurable frame size
frame_size_ms = get_option('default_frame_size_ms')
lib_c_args += ['-DDEFAULT_FRAME_SIZE_MS=' + frame_size_ms.to_string()]

# Configure extension for shared object
if current_os == 'darwin'  # mac
  extension = '.dylib'
elif current_os == 'windows'  # windows
  extension = '.dll'
else  # unix like (linux, bsd)
  extension = '.so'
endif

# Build of the shared object
library('nrepellent',
        common_src,
        noise_repellent_src,
        c_args: lib_c_args,
        name_prefix: '',
        dependencies: all_dep,
        install: true,
        install_dir: install_folder)

library('nrepellent-adaptive',
        common_src,
        noise_repellent_adaptive_src,
        c_args: lib_c_args,
        name_prefix: '',
        dependencies: all_dep,
        install: true,
        install_dir: install_folder)

# Getting version from project configuration or from git tags
version_array = meson.project_version().split('.')

git = find_program('git', native: true, required: false)
if version_array.length() == 0 and git.found()
  git_describe = run_command('git', 'describe', check: false)
  if git_describe.returncode() == 0
    version_array = git_describe.stdout().strip().split('-')[0].split('.')
  endif
endif

# Configure manifest.ttl
manifest_conf = configuration_data()
manifest_conf.set('LIB_EXT', extension)

manifest_ttl = configure_file(
  input: 'lv2ttl/manifest.ttl.in',
  output: 'manifest.ttl',
  configuration: manifest_conf,
  install: true,
  install_dir: install_folder
)

data_conf = configuration_data()
data_conf.set('MAJOR_VERSION', version_array[0])
data_conf.set('MINOR_VERSION', version_array[1])
data_conf.set('MICRO_VERSION', version_array[2])

# Configure nrepellent.ttl
nrepel_ttl = configure_file(
  input: join_paths('lv2ttl', 'nrepellent.ttl.in'),
  output: 'nrepellent.ttl',
  configuration: data_conf,
  install: true,
  install_dir: install_folder
)

# Configure nrepellent#stereo.ttl
nrepel_ttl_stereo = configure_file(
  input: join_paths('lv2ttl', 'nrepellent#stereo.ttl.in'),
  output: 'nrepellent#stereo.ttl',
  configuration: data_conf,
  install: true,
  install_dir: install_folder
)

# Configure nrepellent-adaptive.ttl
nrepel_ttl_adaptive = configure_file(
  input: join_paths('lv2ttl', 'nrepellent-adaptive.ttl.in'),
  output: 'nrepellent-adaptive.ttl',
  configuration: data_conf,
  install: true,
  install_dir: install_folder
)

# Configure nrepellent-adaptive#stereo.ttl
nrepel_ttl_adaptive_stereo = configure_file(
  input: join_paths('lv2ttl', 'nrepellent-adaptive#stereo.ttl.in'),
  output: 'nrepellent-adaptive#stereo.ttl',
  configuration: data_conf,
  install: true,
  install_dir: install_folder
)

# Code formatting target
run_target('format',
  command: [
    'sh', '-c',
    'find "@SOURCE_ROOT@" -name "*.[ch]" | grep -v "/build/" | grep -v "/subprojects/" | xargs clang-format -i && echo "Code formatting complete!"'
  ]
)

# Summary
summary({
  'Version': meson.project_version(),
  'C Standard': get_option('c_std'),
  'Warning Level': custom_warning_level,
  'Warnings as Errors': get_option('treat_warnings_as_errors'),
  'Sanitizers': get_option('enable_sanitizers'),
  'Default Frame Size': frame_size_ms.to_string() + 'ms',
  'Install Directory': install_folder,
}, section: 'Configuration')

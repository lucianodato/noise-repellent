project('libnrepel', 'c', version: '0.1.0',default_options: ['c_std=c99'])

#sources to compile
noise_repellent_lib = [
    'src/shared/spectral_utils.c',
    'src/shared/spectral_features.c',
    'src/shared/noise_profile.c',
    'src/noisemodel/noise_estimator.c',
    'src/denoiser/transient_detector.c',
    'src/denoiser/masking_estimator.c',
    'src/denoiser/spectral_smoother.c',
    'src/denoiser/gain_estimator.c',
    'src/denoiser/spectral_denoiser.c',
    'src/stft/stft_windows.c',
    'src/stft/stft_processor.c',
    'src/noise_repellent.c',
]

# public header
inc = include_directories('include')

m_dep = meson.get_compiler('c').find_library('m', required : true)
fftw_dep = dependency('fftw3f', required : true)
dep = [m_dep, fftw_dep]

#get the host operating system and configure install path and shared object extension
current_os = host_machine.system()
current_arch = host_machine.cpu_family()
cflags = []

#compiler flags
if get_option('buildtype') == 'release'
    # Add x86_64 optimization where appropriate (not for ARM, not for riscv)
    if current_arch != 'aarch64' and current_arch != 'riscv64' and current_os != 'darwin' 
        cflags += ['-ffast-math', '-msse','-msse2','-mfpmath=sse']
        
    endif
    # Add osx multiarch flags when appropriate
    if current_os == 'darwin'
        cflags += ['-mrecip']
    endif
endif

libnrepel = static_library('nrepel',
  noise_repellent_lib,
  include_directories : inc,
  c_args: cflags,
  dependencies: dep,
  install : true)

libnrepel_dep = declare_dependency(include_directories : inc,
  link_with : libnrepel)


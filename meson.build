project('nrepellent.lv2','c',version: '0.1.6',default_options: ['c_std=c99'])

#install folder
if get_option('buildtype') == 'release'
    install_folder = join_paths(get_option('libdir'), 'lv2', meson.project_name())
else
    install_folder = join_paths('.lv2', meson.project_name())
endif

#sources to compile
noise_repellent_sources = [
    'src/spectral_utils.c',
    'src/transient_detector.c',
    'src/masking_estimator.c',
    'src/spectrum_smoother.c',
    'src/noise_estimator.c',
    'src/gain_estimator.c',
    'src/fft_denoiser.c',
    'src/stft_processor.c',
    'src/noise_repellent.c',
    'src/plugin.c'
]

#dependencies for noise repellent
m_dep = meson.get_compiler('c').find_library('m', required : true)
fftw_dep = dependency('fftw3f', required : true)
lv2_dep = dependency('lv2', required : true)
nr_dep = [m_dep,fftw_dep,lv2_dep]

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

# Extension for library by os
if current_os == 'darwin' #mac
    extension = '.dylib'
else #unix like    
    extension = '.so'
endif 
if current_os == 'windows' #windows
    extension = '.dll'
    exe_wrapper = true
endif

#build of the shared object
library('nrepellent',
    noise_repellent_sources,
    name_prefix: '',
    dependencies: nr_dep,
    c_args: cflags, 
    install: true,
    install_dir : install_folder
)
	
#Getting version from project configuration or from git tags
version_array = meson.project_version().split('.')

git = find_program('git', native : true, required : false)
if version_array.length() == 0 and git.found()
    version_array = run_command('git', 'describe').stdout().strip().split('-')[0].split('.')
endif

#Configure manifest.ttl and install
manifest_conf = configuration_data()
manifest_conf.set('LIB_EXT', extension)
manifest_conf.set('MAJOR_VERSION', version_array[0])
manifest_conf.set('MINOR_VERSION', version_array[1])
manifest_conf.set('MICRO_VERSION', version_array[2])

manifest_ttl = configure_file(
    input : 'lv2ttl/manifest.ttl.in',
    output : 'manifest.ttl',
    configuration : manifest_conf,
    install : true,
	install_dir : install_folder
)

#Configure nrepellent.ttl and install
nrepel_ttl = configure_file(
    input : join_paths('lv2ttl', 'nrepellent.ttl.in'),
    output : 'nrepellent.ttl',
    configuration : manifest_conf,
    install : true,
	install_dir : install_folder
)
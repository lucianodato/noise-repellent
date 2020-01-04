project('noise-repellent','c',default_options: 'c_std=c99')

#shared object name
lv2_name = 'nrepellent'

#source to compile
src = 'src/noise_repellent.c'

#dependencies for noise repellent
m_dep = meson.get_compiler('c').find_library('m', required : true)
fftw_dep = dependency('fftw3f', required : true)
lv2_dep = dependency('lv2', required : true)
nr_dep = [m_dep,fftw_dep,lv2_dep]

#get the build operating system and configure install path and shared object extension
current_os = build_machine.system()
current_arch = build_machine.cpu_family()

#compiler global flags
add_global_arguments('-Wno-unused-function', '-ffast-math','-fomit-frame-pointer','-fno-finite-math-only', language : 'c')

#specific cflags
cflags = []

if current_os == 'darwin' #mac
    cflags += ['-mrecip', '-march=i386', '-march=x86-64','-mmacosx-version-min=10.5','-DMAC_OS_X_VERSION_MAX_ALLOWED=105']
endif

# Add x86_64 optimization where appropriate (not for ARM)
if current_arch != 'aarch64'
    cflags += ['-msse','-msse2','-mfpmath=sse']
endif

# Extension for library by os
if current_os == 'darwin' #mac
    extension = '.dylib'
else #unix like    
    extension = '.so'
endif #windows
if current_os == 'windows' or current_os == 'cygwin'
    extension = '.dll'
endif

#install folder
install_folder = 'nrepellent.lv2'

#build of the shared object
shared_library(lv2_name,src,name_prefix: '',dependencies: nr_dep,c_args: cflags,install: true,install_dir : install_folder)

#Configure manifest ttl in order to replace the correct shared object extension
manifest_conf = configuration_data()
manifest_conf.set('LIB_EXT', extension)
manifest_ttl = configure_file(input : 'lv2ttl/manifest.ttl.in',output : 'manifest.ttl',configuration : manifest_conf)

#add manifest.ttl and nrepel.ttl to be installed with the shared object
install_data([manifest_ttl, 'lv2ttl/nrepellent.ttl'],install_dir : install_folder)

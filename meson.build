project('noise-repellent','c',default_options: 'c_std=c99')

#shared object name
lv2_name = 'nrepellent'

#source to compile
src = 'src/plugin.c'

#dependencies for noise repellent
cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : true)
fftw_dep = dependency('fftw3f', required : true)
lv2_dep = dependency('lv2', required : true)
nr_dep = [m_dep,fftw_dep,lv2_dep]

#compiler optimization flags
cflags = ['-msse','-msse2','-mfpmath=sse','-ffast-math','-fomit-frame-pointer','-fno-finite-math-only']

#get the build operating system and configure install path and shared object extension
current_os = build_machine.system()
if current_os == 'darwin' #mac
    i_path = '/Library/Audio/Plug-Ins/LV2/nrepellent.lv2'
    extension = '.dylib'
else #unix like
    i_path = '/usr/local/lib/lv2/nrepellent.lv2'
    extension = '.so'
endif #windows
if current_os == 'windows'
    i_path = '%COMMONPROGRAMFILES%/LV2/nrepellent.lv2'
    extension = '.dll'
endif

#build of the shared object
shared_library(lv2_name,src,name_prefix: '',dependencies: nr_dep,c_args: cflags,install: true,install_dir : i_path)

#Configure manifest ttl in order to replace the correct shared object extension
manifest_conf = configuration_data()
manifest_conf.set('LIB_EXT',extension)
configure_file(input : 'lv2ttl/manifest.ttl.in',output : 'manifest.ttl',configuration : manifest_conf)

#add manifest.ttl and nrepel.ttl to be installed with the shared object
install_data(['build/manifest.ttl', 'lv2ttl/nrepellent.ttl'],install_dir : i_path)
project('noise-repellent','c',default_options: 'c_std=c99')

#shared object name
lv2_name = 'nrepel'

#source to compile
src = 'src/nrepel.c'

#dependencies for noise repellent
cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : true)
fftw_dep = dependency('fftw3f', required : true)
lv2_dep = dependency('lv2', required : true)
nr_dep = [m_dep,fftw_dep,lv2_dep]

#compiler optimization flags
if meson.get_compiler('c').get_id() == 'clang'
  add_global_arguments('-mrecip', language : 'c')
endif
cflags = ['-msse','-msse2','-mfpmath=sse','-ffast-math','-fomit-frame-pointer','-fno-finite-math-only']
cflags += ['-Wno-unused-function']

#install folder
install_folder = 'nrepel.lv2'

#get the build operating system and configure install path and shared object extension
current_os = build_machine.system()
if current_os == 'darwin' #mac
    extension = '.dylib'
else #unix like
    #i_path = '/usr/local/lib/lv2/nrepel.lv2'
    extension = '.so'
endif #windows
if current_os == 'windows'
    extension = '.dll'
endif

#build of the shared object
shared_library(lv2_name,src,name_prefix: '',dependencies: nr_dep,c_args: cflags,install: true,install_dir : install_folder)

#Configure manifest ttl in order to replace the correct shared object extension
manifest_conf = configuration_data()
manifest_conf.set('LIB_EXT',extension)
manifest_ttl = configure_file(input : 'lv2ttl/manifest.ttl.in',output : 'manifest.ttl',configuration : manifest_conf)

#add manifest.ttl and nrepel.ttl to be installed with the shared object
install_data([manifest_ttl, 'lv2ttl/nrepel.ttl'],install_dir : install_folder)

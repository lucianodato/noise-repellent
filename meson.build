project(
    'noise-repellent',
    'c',
    default_options: 'c_std=c99'
)

#shared object name
lv2_name = 'nrepel'

#source to compile
src = 'src/nrepel.c'

#dependencies for noise repellent
fftw_dep = dependency('fftw3f')
lv2_dep = dependency('lv2')
nr_dep = [
    fftw_dep,
    lv2_dep
]

#compiler optimization flags
cflags = [
    '-msse',
    '-msse2',
    '-mfpmath=sse',
    '-ffast-math',
    '-fomit-frame-pointer',
    '-fno-finite-math-only'
]

#linking flags (besides dependencies)
largs = '-lm'

#install destination
i_path = '/usr/local/lib/lv2/nrepel.lv2'

#build of the shared object
shared_library(
    lv2_name,
    src,
    name_prefix: '',
    dependencies: nr_dep,
    c_args: cflags,
    link_args: largs,
    install: true,
    install_dir : i_path
)

#get the build operating system
current_os = build_machine.system()
if current_os == 'darwin'
  extension = '.dylib' #mac
else
  extension = '.so' #linux and others
endif
if current_os == 'windows'
  extension = '.dll' #windows
endif

#Configure manifest ttl in order to replace the correct shared object extension
manifest_conf = configuration_data()
manifest_conf.set('LIB_EXT',extension)
configure_file(input : 'lv2ttl/manifest.ttl.in',
  output : 'manifest.ttl',
  configuration : manifest_conf)

#add manifest.ttl and nrepel.ttl to be installed with the shared object
install_data(
    ['build/manifest.ttl', 'lv2ttl/nrepel.ttl'],
    install_dir : i_path
)

project(
    'noise-repellent',
    'c',
    default_options: 'c_std=c99'
)

#shared object name
lv2_name = 'nrepel'

#source to compile
src = 'src/nrepel.c'

#install destination
i_path = '/usr/local/lib/lv2/nrepel.lv2'

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
largs = [
    '-lm'
]

#Ensamble ttl files in order to be copied to install folder
#meson.add_install_script('../generate_ttls.sh')

#build configuration
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
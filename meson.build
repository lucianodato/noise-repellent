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
largs = [
    '-lm'
]

#install destination
i_path = '/usr/local/lib/lv2/nrepel.lv2'

#Ensamble ttl files in order to be copied to install folder
#version = run_command('git', 'describe', '--tags HEAD 2>/dev/null | sed "s/-g.*$$//;s/^v//" || echo "LV2"').stdout().strip()
#lib_ext = get_option(name_suffix)
#ttl_dir = 'lv2ttl'

#sed "s/@LV2NAME@/lv2_name/;s/@LIB_EXT@/lib_ext/" \
#    ttl_dir/manifest.ttl.in > i_path/manifest.ttl

#sed "s/@VERSION@/lv2:microVersion $(LV2MIC) ;lv2:minorVersion $(LV2MIN) ;/g" \
#    ttl_dir/lv2_name.ttl.in > i_path/lv2_name.ttl

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
project('noise-repellent', 'c',
    default_options : 'c_std=c99',
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

#build configuration
shared_library(lv2_name,
    src,
    dependencies: nr_dep,
    c_args: cflags,
    name_prefix: ''
)


project('nrepel.lv2','c',version: '0.1.6',default_options: ['c_std=c99'])

sord_validate = find_program('sord_validate', native : true, required : false)
cp = find_program('cp')
clone = [cp, '@INPUT@', '@OUTPUT@']

#source to compile
src = 'src/nrepel.c'

#dependencies for noise repellent
cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : true)
fftw_dep = dependency('fftw3f', required : true)
lv2_dep = dependency('lv2', required : true)
nr_dep = [m_dep,fftw_dep,lv2_dep]

#compiler default flags
add_global_arguments('-fomit-frame-pointer','-fno-finite-math-only','-Wno-unused-function',language : 'c')

#install folder
install_folder = join_paths(get_option('libdir'), 'lv2', meson.project_name())

#get the host operating system and configure install path and shared object extension
current_os = host_machine.system()
current_arch = host_machine.cpu_family()
current_compiler = meson.get_compiler('c')
cflags = []

# Add x86_64 optimization where appropriate (not for ARM, not for riscv)
if current_arch != 'aarch64' and current_arch != 'riscv64' and current_os != 'darwin' 
    cflags += ['-ffast-math', '-msse','-msse2','-mfpmath=sse']
endif

# Add osx multiarch flags when appropriate
if current_os == 'darwin'
    cflags += ['-mrecip']
endif

# Extension for library by os
if current_os == 'darwin' #mac
    extension = '.dylib'
else #unix like    
    extension = '.so'
endif #windows
if current_compiler.contains('mingw')
    extension = '.dll'
endif

#build of the shared object
library('nrepel',
    src,
    name_prefix: '',
    dependencies: nr_dep,
    c_args: cflags, 
    install: true,
    install_dir : install_folder
)
	
#Getting version from project configuration or from git tags
version_array = meson.project_version().split('.')
if version_array.length() == 0
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

#Configure nrepel.ttl and install
nrepel_ttl = configure_file(
    input : join_paths('lv2ttl', 'nrepel.ttl.in'),
    output : 'nrepel.ttl',
    configuration : manifest_conf,
    install : true,
	install_dir : install_folder
)

# testing ttls
if sord_validate.found()
	test('LV2 validation', sord_validate,
		args : [run_command('find','./lv2 -name "*.ttl"').stdout(), run_command('find','. -name "*.ttl"').stdout()])
endif

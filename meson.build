project('nrepel.lv2','c',default_options: 'c_std=c99')

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
add_global_arguments('-ffast-math','-fomit-frame-pointer','-fno-finite-math-only','-Wno-unused-function',language : 'c')

#install folder
install_folder = join_paths(get_option('libdir'), 'lv2', meson.project_name())

#get the build operating system and configure install path and shared object extension
current_os = build_machine.system()
current_arch = build_machine.cpu_family()
cflags = []

# Add x86_64 optimization where appropriate (not for ARM)
if current_arch != 'aarch64'
    cflags += ['-msse','-msse2','-mfpmath=sse']
endif

# Add osx multiarch flags when appropriate
if current_os == 'darwin'
    cflags += ['-mrecip',
            '-march=i386',
            '-march=x86-64',
            '-mmacosx-version-min=10.5',
            '-DMAC_OS_X_VERSION_MAX_ALLOWED=105']
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

#build of the shared object
shared_library('nrepel',
        src,
        name_prefix: '',
        dependencies: nr_dep,
        c_args: cflags, 
        install: true,
        install_dir : install_folder
        )

#Configure manifest ttl in order to replace the correct shared object extension and install
manifest_conf = configuration_data()
manifest_conf.set('LIB_EXT', extension)
manifest_ttl = configure_file(
    input : 'lv2ttl/manifest.ttl.in',
    output : 'manifest.ttl',
    configuration : manifest_conf,
    install : true,
	install_dir : install_folder
    )

#add nrepel.ttl to be installed with the shared object
nrepel_ttl = custom_target('nrepel_ttl',
	input : join_paths('lv2ttl', 'nrepel.ttl'),
	output : 'nrepel.ttl',
	command : clone,
	install : true,
	install_dir : install_folder)

# testing ttls
if sord_validate.found()
	test('LV2 validation', sord_validate,
		args : [run_command('find','./lv2 -name "*.ttl"').stdout(), run_command('find','. -name "*.ttl"').stdout()])
endif
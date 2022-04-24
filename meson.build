project('nrepellent.lv2','c',version: '0.2.2',default_options: ['default_library=shared','c_std=c99'])

#install folder
if(get_option('buildtype') == 'release')
    lv2_directory = join_paths(get_option('libdir'), 'lv2')
else
    lv2_directory = '.lv2'
endif
install_folder = join_paths(lv2_directory, meson.project_name())

# sources to compile
common_src = ['src/signal_crossfade.c']
noise_repellent_src = ['plugins/nrepellent.c', 'src/noise_profile_state.c']
noise_repellent_adaptive_src = 'plugins/nrepellent-adaptive.c'

#dependencies for noise repellent
lv2_dep = dependency('lv2', required: true)
libspecbleach_dep = dependency('libspecbleach', fallback : ['libspecbleach', 'libspecbleach_dep'], default_options: ['default_library=static'], required: true)
m_dep = meson.get_compiler('c').find_library('m', required: true)
all_dep = [lv2_dep,libspecbleach_dep,m_dep]

#get the host operating system and configure install path and shared object extension
current_os = host_machine.system()

if current_os == 'darwin' #mac
    extension = '.dylib'
else #unix like    
    extension = '.so'
endif 
if current_os == 'windows' #windows
    extension = '.dll'
endif

#build of the shared object
library('nrepellent',
    common_src,
    noise_repellent_src,
    name_prefix: '',
    dependencies: all_dep,
    install: true,
    install_dir: install_folder
)

library('nrepellent-adaptive',
    common_src,
    noise_repellent_adaptive_src,
    name_prefix: '',
    dependencies: all_dep,
    install: true,
    install_dir: install_folder
)
	
#Getting version from project configuration or from git tags
version_array = meson.project_version().split('.')

git = find_program('git', native: true, required: false)
if version_array.length() == 0 and git.found()
    version_array = run_command('git', 'describe').stdout().strip().split('-')[0].split('.')
endif

#Configure manifest.ttl
manifest_conf = configuration_data()
manifest_conf.set('LIB_EXT', extension)

manifest_ttl = configure_file(
    input: 'lv2ttl/manifest.ttl.in',
    output: 'manifest.ttl',
    configuration: manifest_conf,
    install: true,
	install_dir: install_folder
)

data_conf = configuration_data()
data_conf.set('MAJOR_VERSION', version_array[0])
data_conf.set('MINOR_VERSION', version_array[1])
data_conf.set('MICRO_VERSION', version_array[2])

#Configure nrepellent.ttl
nrepel_ttl = configure_file(
    input: join_paths('lv2ttl', 'nrepellent.ttl.in'),
    output: 'nrepellent.ttl',
    configuration: data_conf,
    install: true,
	install_dir: install_folder
)

#Configure nrepellent.ttl
nrepel_ttl = configure_file(
    input: join_paths('lv2ttl', 'nrepellent#stereo.ttl.in'),
    output: 'nrepellent#stereo.ttl',
    configuration: data_conf,
    install: true,
	install_dir: install_folder
)

#Configure nrepellent-adaptive.ttl
nrepel_ttl_adaptive = configure_file(
    input: join_paths('lv2ttl', 'nrepellent-adaptive.ttl.in'),
    output: 'nrepellent-adaptive.ttl',
    configuration: data_conf,
    install: true,
	install_dir: install_folder
)

#Configure nrepellent-adaptive#stereo.ttl
nrepel_ttl_adaptive_stereo = configure_file(
    input: join_paths('lv2ttl', 'nrepellent-adaptive#stereo.ttl.in'),
    output: 'nrepellent-adaptive#stereo.ttl',
    configuration: data_conf,
    install: true,
	install_dir: install_folder
)
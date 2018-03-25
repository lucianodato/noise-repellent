Name:           noise-repellent
Version:        0.1.4
Release:        1%{?dist}
Summary:        An lv2 plugin for broadband noise reduction

License:        LGPLv3
URL:            https://github.com/lucianodato/noise-repellent
Source0:        https://github.com/lucianodato/noise-repellent/archive/0.1.4.tar.gz

BuildRequires:  fftw3
BuildRequires:  lv2

%description
Noise-repellent is a spectral noise reduction lv2 plugin based on spectral subtraction with
a bunch of improvements to avoid artifacts. 

%prep
%autosetup


%build
%make_build


%install
rm -rf $RPM_BUILD_ROOT
%make_install


%files
%license https://raw.githubusercontent.com/lucianodato/noise-repellent/master/LICENSE
%doc https://raw.githubusercontent.com/lucianodato/noise-repellent/master/README.md



%changelog
* Sun Mar 25 2018 Luciano <lucianodato@gmail.com>
- 

# .spec file for building an rpm of coset
# written by Paul D. Boyle, University of Western Ontario

Summary: A program to perform left coset decompositions using Flack's algorithms
Name: coset
Version: 1.3.0
Release: 1
License: GPL
Group: Applications/Scientific/Crystallography
Source: http://xray.uwo.ca/COSET/coset-1.3.0.tar.gz
URL:http://xray.uwo.ca/COSET/
Packager: Paul D. Boyle <pboyle@uwo.ca>
Prefix: /usr/local
BuildRoot:%{_tmppath}%{name}-%{version}-root

%_tmppath /tmp

%description
COSET is a program which performs left coset decompositions according
to the algorithms given in Acta Cryst. (1987), A43, 564-568, by
H. D. Flack.  The program is used for deriving twin laws for merohedrally
or pseudomerohedrally twinned crystals.  The program can also create SHELX
.ins file which incorporate the appropriate TWIN instructioni, and can, if
a locally available copy of SHELXL is available, test the twin law using
a least-squares refinement. See the README.txt file for more information.

%prep
%setup

%build
make

%install
make install

%files
%doc /usr/local/doc/COSET_README.txt
/usr/local/bin/coset

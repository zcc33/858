#!/usr/local/bin/perl
$ENV{'LD_LIBRARY_PATH'}="/usr/local/stow/matlabr2010b/lib/matlabr2010b/tomlab/shared";
$ENV{'TOMLAB_LICENSE_FILE'}= "/usr/local/lib/matlabr2010b/tomlab/tomlab.lic";
$cmd = "matlabr2010b -nodisplay -r \"x = pwd; cd /usr/local/stow/matlabr2010b/lib/matlabr2010b/tomlab; startup; cd (x);\"";
system($cmd);

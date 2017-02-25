#!/usr/bin/perl
#
# Utility to replace a specific string in a set of Fortran files
# Alexei Pankin, 2002
#
$location1='../pedestal';
$location2='.';
opendir(THISDIR, $location1) || die "Unable to read the directory!\n";
@allfiles=grep(/(dtransp\.f$|\.f90$|\.F$|\.F90$|\.for$|\.FOR$)/, readdir(THISDIR));
closedir(THISDIR);
foreach $filename (@allfiles) {
  open(FFILE1, "<$location1\/$filename") || print "unable to open file $filename for reading\n";
  open(FFILE2, ">$location2\/$filename") || print "unable to open file $filename for writting\n";
  while (<FFILE1>) {
    #(s/call reset[ir]\((\w+),55,0\.0\)/\1 = 0.0/i);
    (s/call include \'\.\.\/com\//include \'/i);
    printf(FFILE2 $_); 
  };
  close(FFILE1);
  close(FFILE2);
};

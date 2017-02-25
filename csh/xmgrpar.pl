#!/usr/bin/perl
#
# Utility to replace a specific string in a set of Fortran files
# Alexei Pankin, 2003
#
$location1=$ENV{'CODESYSDIR'}.'/data';
$location2='.';
opendir(THISDIR, $location1) || die "Unable to read the directory!\n";
@allfiles=grep(/plot\.par$/, readdir(THISDIR));
closedir(THISDIR);
foreach $filename (@allfiles) {
  open(FFILE1, "<$location1\/$filename") || print "unable to open file $filename for reading\n";
  open(FFILE2, ">$location2\/$filename") || print "unable to open file $filename for writting\n";
  while (<FFILE1>) {
#    #(s/call reset[ir]\((\w+),55,0\.0\)/\1 = 0.0/i);
    (s/BALDUR simulation/BALDUR simulations of discharge $ARGV[0]/i);
    (s/______/ $ARGV[1]/i);
    (s/sim 1/$ARGV[2]/i);
    (s/sim 2/$ARGV[3]/i);
    (s/sim 3/$ARGV[4]/i);
    (s/sim 4/$ARGV[5]/i);
    (s/sim 5/$ARGV[6]/i);
    (s/sim 6/$ARGV[7]/i);
    (s/sim 7/$ARGV[8]/i);
    (s/sim 8/$ARGV[9]/i);
    (s/sim 9/$ARGV[10]/i);
    printf(FFILE2 $_); 
  };
  close(FFILE1);
  close(FFILE2);
}

#!/bin/csh
# 
# Alexei Pankin, Lehigh Univ.
# pankin@fusion.physics.lehigh.edu
#
# to run the script use the command
#
# ./findlhmode.sh ID run_id1 run_id2 ...
#
# where ID is id used for set of the runs (discharge number,
# for example), run_id1 and run_id2 are specific run id is series of
# simulations (for example a01, a02, etc).
set shot = $argv[1]
shift
foreach dr ($argv[*])
set dt1=`date`
set runid = $shot$dr
echo "Checking $runid case ... " `date`
cd $dr
@ hmode1=`fgrep -c '1  (t' j$runid`
@ hmode2=`fgrep -c '1  (P' j$runid`
@ lmode1=`fgrep -c '0  (t' j$runid`
@ lmode2=`fgrep -c '0  (P' j$runid`
echo 'Guzdar model L-mode/H-mode ' $lmode1/$hmode1 ' cases' > lhmode.log
echo 'Emprical model L-mode/H-mode ' $lmode2/$hmode2 ' cases' >> lhmode.log
echo `fgrep 'Mode =' j$runid` | sed 's/time/#time/g' | tr "#" "\n"  >> lhmode.log
cd ../
end


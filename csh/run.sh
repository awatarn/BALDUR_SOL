#!/bin/csh
# 
# Alexei Pankin, Lehigh Univ.
# pankin@fusion.physics.lehigh.edu
#
# to run the script use the command 
# 
# ./run.sh ID run_id1 run_id2 ...
#
# where ID is id used for set of the runs (discharge number, 
# for example), run_id1 and run_id2 are specific run id is series of
# simulations (for example a01, a02, etc). 
#
# For example, to run four test cases use the command 
#
# ./run.sh test a01 a02 a03 a04 &
#
# The input file for each run should be in a separate directory. 
# The name of input file should have the structure IDrun_id.
#
if (! -e notes) then
 echo > notes
endif
if (! -e run.log) then
  echo > run.log
endif  
set shot = $argv[1]
shift
foreach dr ($argv[*])
set dt1=`date`
set runid = $shot$dr
echo "Running $runid case ... " `date` >> run.log 
echo "${runid}" >> notes
cd $dr
cp $runid jobxdat 
cp $CODESYSDIR/data/for22 .
$CODESYSDIR/bin/xbaldur
rm for22
echo "Finished ... " `date` >> ../run.log
mv jobxlpt j$runid
rm jobxdat
cp /home/simruns/sim_aw/script/callSlist.sh .
set dt2=`date`
echo "Start time:" $dt1 > ../notes.tmp
echo "Finish time:" $dt2 >> ../notes.tmp
echo "XBALDUR v." `fgrep '1-1/2-D BALDUR' j$runid | sed 's#.*version##'` >> ../notes.tmp
echo "------------" >> ../notes.tmp
head -5 $runid >> ../notes.tmp
@ skip=`fgrep -c 'te-axis' j$runid` 
echo `fgrep 'te-axis' j$runid | sed -n ${skip}~1p`  >> ../notes.tmp
echo `fgrep 'te-avg.' j$runid | sed -n ${skip}~1p`  >> ../notes.tmp
echo `fgrep 'q-axis' j$runid | sed -n ${skip}~1p`   >> ../notes.tmp
echo `fgrep 'Walpha=' j$runid | sed -n ${skip}~1p`>> ../notes.tmp
echo `fgrep 'Peloss=' j$runid | sed -n ${skip}~1p`>> ../notes.tmp 
echo `fgrep 'Piloss=' j$runid | sed -n ${skip}~1p`>> ../notes.tmp
echo `fgrep 'Pheat=' j$runid | sed -n ${skip}~1p` >> ../notes.tmp
echo `fgrep 'ITER89-P,' j$runid | sed -n ${skip}~1p`>> ../notes.tmp
echo `fgrep 's ne-bar=' j$runid | sed -n ${skip}~1p`>> ../notes.tmp
tail -5 j$runid >> ../notes.tmp
echo "====================================================================================" >> ../notes.tmp
cd ../
/home/simruns/sim_aw/script/callSlist.sh $shot $dr
mail -s "BALDUR run you submitted at $dt1 is completed at $dt2" $MAIL_ADDRESS < notes.tmp
cat notes.tmp >> notes
rm -f notes.tmp
end


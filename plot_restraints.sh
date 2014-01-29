#!/bin/bash 
echo "#Num	total	fix	slvnt_rad	slvnt_pol	shell solute" >restraints.dat
echo "#Num	el	vdW	bond	angle	proper	improper">solute.dat

if [ -f restraints2.dat ];then
rm restraints2.dat
fi
if [ -f solute2.dat ];then
rm solute2.dat
fi
#Enter the name for the files you wan to plot here. If a file is not found it will be ignored
FILELIST=""
touch restraints2.dat solute2.dat

for step in ${FILELIST} 
do
	if [ -f ${step}.log ];then
	grep "^restraints" ${step}.log >>restraints2.dat 
	grep "^solute " ${step}.log >>solute2.dat
	fi
done

#This part looks for files from a FEP or EVB run (fep_100.log to fep_000.log) and tries to add them to the plot.
#Again, if a file is missing it will be ignored
for (( c=100; c>=0; c-- )); do
NUM=${c}
if [[ ${c} -lt 100 ]];then
NUM=0${c}
fi
if [[ ${c} -lt 10 ]];then
NUM=00${c}
fi
MODULUS=$(( ${c}%2 ))
if [[ $MODULUS = 0 ]];then
 if [ -f fep_${NUM}.log ] ;then
 grep "^restraints" fep_${NUM}.log >>restraints2.dat
 grep "^solute " fep_${NUM}.log >>solute2.dat
 fi
fi
done
nl restraints2.dat >>restraints.dat
nl solute2.dat >>solute.dat

#Pipes the command to Gnuplot
#Make sure you are working locally or are using a ssh -X session
echo "set terminal x11 enhanced size 1920,1200
set size 1,1
set xlabel 'Step [x500]'
set ylabel 'Energy [kcal/mol]'
plot \"restraints.dat\" using 1:3 with lines ls 1 title \"total\" , \\
\"restraints.dat\" using 1:4 with lines ls 2 title \"restraints\" , \\
\"restraints.dat\" using 1:5 with lines ls 3 title \"solv-rad\" , \\
\"restraints.dat\" using 1:6 with lines ls 4 title \"solv-pol\" , \\
\"restraints.dat\" using 1:7 with lines ls 5 title \"shell\" " | gnuplot -presist -

echo "set terminal x11 enhanced size 1920,1200
set size 1,1
set xlabel 'Step [x500]'
set ylabel 'Energy [kcal/mol]'
plot \"solute.dat\" using 1:5 with lines ls 1 title \"bond\", \\
\"solute.dat\" using 1:6 with lines ls 2 title \"angle\", \\
\"solute.dat\" using 1:7 with lines ls 3 title \"proper\", \\
\"solute.dat\" using 1:8 with lines ls 4 title \"improper\" " | gnuplot -presist -

#rm restraints2.dat solute2.dat

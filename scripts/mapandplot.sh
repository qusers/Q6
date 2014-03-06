#!/bin/bash -l
#Paul Bauer, 2014-01-21
#Script to evaluate a Q EVB run and obtain reaction free energies, activation energies
#Distances of reacting atoms and an evaluation of the sampling rate
#Released under the Beer-Ware License

#You can iterate over directories and other stuff
#if the run is in only one directory, set the LIST to 'empty'
#Iterate over the directories needed and fill the list
LIST=""
for ((dirnum=1;dirnum<=10;dirnum++))
do
LIST="${LIST} frame-${dirnum}"
done

#Logical variables
TRUE='true'
FALSE='false'

#If you just want to do the anaylsis again, set RERUN to true
RERUN=${FALSE}
#Path to Qfep binary, needed in any case
qfep='/home/p/pabau/pfs/Q_PFS_MPIFIX/qfep5'
#Function that generates the temporary awk file needed to get the distances
#Currently five is the maximum, need to find a way to change this to be independent
#Also generates the awk file for the mapping analysis and the energy printout
#That is just a table and the statistics so I can copy the numbers directly
#Hacksih in the way I like it :D
#So they don't have to be copied any more :P
function make_temp_awk() {
echo "#!/usr/bin/awk -f
{

if (NR==1) {
COUNT=1
MCOUNT=0
}

        STOR[COUNT]=\$1
	STORA[COUNT]=\$2
        if (COUNT==${NUMDIST}) {
                COUNT=0
                print \" \",MCOUNT,\" \",STORA[1],\" \",STORA[2],\" \",STORA[3],\" \",STORA[4],\" \",STORA[5]
                MCOUNT++
                }
	COUNT++
}">$1
echo "#!/usr/bin/awk -f

{
        if ( startnow==1 ) {
                print \" \",\$2,\"      \",\$4,\"  \",\$8,\"      \",\$3
        }

        if ( \$2==\"bin\" && \$3==\"energy\" && \$4==\"gap\" ) {
                startnow++
                print \"#        energygap       dG(norm)        r_xy    dG(nonnorm)\"
        }
}">getenergy.awk
echo "#!/usr/bin/awk -f

{
        if ( \$2==\"Part\" && \$3==\"3:\" ) {
                startnow=2
        }


        if ( startnow==1 ) {
                print \" \",\$1,\"      \",\$3,\"  \",\$7
        }

        if ( \$2==\"Lambda(1)\" && \$4==\"Energy\" && \$5==\"gap\" ) {
                startnow=1
                print \"#        Lambda  energygap       pts\"
        }

}
">getlambda.awk
#Now using the energies printed by getenergy.awk to get the reaction free enrgy and activation energy
echo "#!/usr/bin/awk -f
{
	if (NR==1) {
		COUNTONE=1
		COUNTTWO=1
		MINARRAY[COUNTONE]=0
		MINIMAS[COUNTONE]=3000000000
		MAXIMAS[COUNTTWO]=-3000000000
		MAXARRAY[COUNTTWO]=0
		ENERGY[NR]=300000000
		DGSTAR[COUNTONE]=-10000
		DGFREE[COUNTONE]=10000
	}
	if (NR>=2) {
		ENERGY[NR]=\$2
		if ((ENERGY[NR] > ENERGY[NR-1]) && (MINARRAY[COUNTONE]==0)) {
			MINIMAS[COUNTONE]=ENERGY[NR-1]
#			print \"COUNTONE\", COUNTONE,\"Value\",MINIMAS[COUNTONE]
			MINARRAY[COUNTONE]=1
			if (COUNTONE>1) {
				DGFREE[COUNTTWO]=MINIMAS[COUNTONE]-MINIMAS[COUNTONE-1]
				COUNTTWO++
				MAXARRAY[COUNTTWO]=0
			}
		}
		if ((ENERGY[NR] < ENERGY[NR-1]) && (MINARRAY[COUNTONE]==1) && (MAXARRAY[COUNTTWO]==0)) {
			MAXIMAS[COUNTTWO]=ENERGY[NR-1]
#			print \"COUNTTWO\", COUNTTWO,\"Value\",MAXIMAS[COUNTTWO]
			MAXARRAY[COUNTTWO]=1
			DGSTAR[COUNTONE]=MAXIMAS[COUNTTWO]-MINIMAS[COUNTONE]
			COUNTONE++
			MINARRAY[COUNTONE]=0
		}
	}
}
END {
	print \"Number of minimas = \",COUNTONE-1
	print \"Number of maximas = \",COUNTTWO-1
	for (tmpcount=1;tmpcount<=COUNTONE-1;tmpcount++) {
		print \"dG* Minima \",tmpcount,\" Maximum \",tmpcount,\" = \",DGSTAR[tmpcount]
	}
	for (tmpcount=1;tmpcount<=COUNTTWO-1;tmpcount++) {
		print \"dGFree Minimum \",tmpcount,\" Minimum \", tmpcount+1,\" = \", DGFREE[tmpcount]
	}
}">getdeltag.awk
}
#This function makes the static inputs that are always the same for all Gnuplot scripts
#So we only need to write the stuff once and not a number of times
#It takes the filename of the script as input to know where to write to and the terminal type
#And also how many plots should be combined
function gnuplot_static() {
tnum=$( echo "scale=2;1/$3 " | bc -q)
echo "set terminal $2 enhanced size 1920,1200
set size 1,1
set origin 0.0, 0.0
set multiplot layout $3,1
set size 1,${tnum}
set style line 1 lt 1 lc 1 lw 5 pt 6
set style line 2 lt 1 lc 2 lw 5 pt 6
set style line 3 lt 1 lc 3 lw 5 pt 6
set style line 4 lt 1 lc 4 lw 5 pt 6
set style line 5 lt 1 lc 5 lw 5 pt 6
set style line 6 lt 1 lc 6 lw 5 pt 6
set key inside left top enhanced autotitles box linetype -1 linewidth 1.000">$1
}

#This function will generate the plot of the energies of the run from the energy*dat files
#It takes the name of the file to write to as the first argument, the plot coordinates as the second and third,
#the label names as fourth and fifth and the lines to use as sixth argument. The actual name of the file
#(energy or lambda) is the seventh and final argument
#The function iterates over the number of files mapped before and adds them to the plot
#If the last file is reached, the plot command is finished
#If only one file is being plotted, a one line command is issued only
#The linestyle is changed using the counter in LINEC
function gnuplot_plot() {
if [[ ${7} != "energy" ]];then
LINES=""
else
LINES="with lines"
fi
echo "set xlabel '$4'
set ylabel '$5'
set origin $2, $3">>$1
LINEC=1
for i in ${LIST} ;do
	
if [[ ${COUNTER} -eq 1 ]];then
	echo "plot '${7}_${i}.dat' using 1:$6 title '${i}'  ${LINES} ls ${LINEC} ">>$1
else
	if [[ ${LINEC} -eq 1 ]];then
		echo "plot '${7}_${i}.dat' using 1:$6 title '${i}'  ${LINES} ls ${LINEC} ,\\">>$1
	elif [[ ${LINEC} -eq ${COUNTER} ]];then
		echo " '${7}_${i}.dat' using 1:$6 title '${i}'  ${LINES} ls ${LINEC}">>$1
	else
		echo " '${7}_${i}.dat' using 1:$6 title '${i}'  ${LINES} ls ${LINEC},\\">>$1
	fi
	LINEC=$(( ${LINEC} + 1 ))
fi
done
}
#Special function for the distance plotting is needed, as we do not know how
# many distances are there in the beginning, so we check that dynamically
# arguments are the same as above, but with the extra one for the number of distances and without
# the fixed column number
function gnuplot_plot_distance() {
echo "set xlabel '$4'
set ylabel '$5'
set origin $2, $3">>$1
LINEC=1
LINEN=1
for i in ${LIST} ;do
if [[ ${COUNTER} -eq 1 ]] && [[ ${NUMDIST} -eq 1 ]];then
	echo "plot '${6}_${i}.dat' using 1:2 title '${i} distance 1' ${LINES} ls ${LINEC} ">>$1
else 
	for ((numc=1;numc<=${NUMDIST};numc++))
	do
	if [[ ${numc} -eq 1 ]] && [[ ${LINEC} -eq 1 ]] ; then
	echo "plot '${6}_${i}.dat' using 1:$(( ${numc} + 1)) title '${i} distance ${numc}' ${LINES} ls ${LINEC} ,\\">>$1
	elif [[ ${numc} -eq ${NUMDIST} ]] && [[ ${LINEN} -eq ${COUNTER} ]]; then
	echo "'${6}_${i}.dat' using 1:$(( ${numc} + 1)) title '${i} distance ${numc}' ${LINES} ls ${LINEC}">>$1
	else
	echo "'${6}_${i}.dat' using 1:$(( ${numc} + 1)) title '${i} distance ${numc}' ${LINES} ls ${LINEC} ,\\">>$1
	fi
	LINEC=$(( ${LINEC} + 1 ))
	done
fi
LINEN=$(( ${LINEN} + 1 ))		
done

}
#This function generates the final result table for the free energies obtained from the mapping
#It also does some statistics, depending on the number of replicas available
function make_energy_table() {
#Printing the header
echo "Results of the free energy calculation"
echo "Only showing the first minima, the others are in the summary.dat files"
echo "#Run	dGFree		dG*"
AFVALUE=0
ASVALUE=0
STFVALUE=0
STSVALUE=0
for ((ncount=1;ncount<=${COUNTER};ncount++))
do
echo "${ncount}	${DGFREESTOR[${ncount}]}		${DGSTARSTOR[${ncount}]}"
done
if [ ${COUNTER} -ge 2 ] ; then
for ((ncount=1;ncount<=${COUNTER};ncount++))
do
AFVALUE=$( echo "scale=2;${AFVALUE} + ${DGFREESTOR[${ncount}]} " | bc -q)
ASVALUE=$( echo "scale=2;${ASVALUE} + ${DGSTARSTOR[${ncount}]} " | bc -q)
done
AFVALUE=$( echo "scale=2;${AFVALUE}/${COUNTER} " | bc -q)
ASVALUE=$( echo "scale=2;${ASVALUE}/${COUNTER} " | bc -q)
echo "Average free energy			=	${AFVALUE}"
echo "Average activation energy		=	${ASVALUE}"
fi
if [ ${COUNTER} -ge 3 ] ; then
for ((ncount=1;ncount<=${COUNTER};ncount++))
do
TMPVALUE=$( echo "scale=5;${AFVALUE} - ${DGFREESTOR[${ncount}]} " | bc -q)
TMPVALUE=$( echo "scale=5;${TMPVALUE} * ${TMPVALUE} " | bc -q)
STFVALUE=$( echo "scale=5;${TMPVALUE} + ${STFVALUE} " | bc -q)
TMPVALUE=$( echo "scale=5;${ASVALUE} - ${DGSTARSTOR[${ncount}]} " | bc -q)
TMPVALUE=$( echo "scale=5;${TMPVALUE} * ${TMPVALUE} " | bc -q)
STSVALUE=$( echo "scale=5;${TMPVALUE} + ${STSVALUE} " | bc -q)
done
DEGREEFREDOM=$(( ${COUNTER} - 1))
STFVALUE=$( echo "scale=5;sqrt(${STFVALUE} / ${DEGREEFREDOM}) " | bc -q)
STSVALUE=$( echo "scale=5;sqrt(${STSVALUE} / ${DEGREEFREDOM}) " | bc -q)
echo "Free energy std. deviation		=	${STFVALUE}"
echo "Activation energy std. deviation	=	${STSVALUE}"
fi
}




#To see how many files are plotted, needed to generate Gnuplot input later on
COUNTER=0
for i in ${LIST} ;do
if [[ ${LIST} != "empty" ]] ; then
cp map_files.inp ${i}
cd ${i}
fi
if [[ ${RERUN} == ${FALSE} ]] ; then
echo "Running qfep in ${PWD}"
${qfep} <map_files.inp>map_files.out
fi
echo "qfep is done, analysing data"
#Hackish way to get the distances directly from the fep files by iterating over them
#Same way as the iteration goes in the frame generation and running
if [ -f distances.dat ] ; then
rm distances.dat
fi
checknum=0
for ((num=100;num>=0;num--))
do
NEWMOD=$(( ${num}%2 ))
if [[ ${NEWMOD} = 0 ]] ; then
if [[ ${num} -ge 100 ]] ; then
ENUM=${num}
elif [[ ${num} -ge 10 ]] ; then
ENUM=0${num}
else
ENUM=00${num}
fi
grep "dist. between Q-atoms" fep_${ENUM}.log |  awk 'NF { print ( $(NF) ) }' | nl >>distances.dat
sed -i s/"="//g distances.dat
fi
done
if [[ ${checknum} -eq 0 ]] ; then
tail -50 fep_${ENUM}.log | grep "dist. between Q-atoms" |  awk 'NF { print ( $(NF) ) }' | nl >tdis.dat
NUMDIST=`awk 'END { print $1 }' tdis.dat ` 
#echo "NUMDIST = ${NUMDIST}"
checknum=1337
rm tdis.dat
fi
make_temp_awk tmp.awk
awk -f tmp.awk <distances.dat>distance_${i}.dat
awk -f getenergy.awk <map_files.out>energy_${i}.dat
awk -f getlambda.awk <map_files.out>lambda_${i}.dat
awk -f getdeltag.awk <energy_${i}.dat>summary_${i}.dat



rm distances.dat *.awk
#Moves the files back for plotting everything together
COUNTER=$(( ${COUNTER} + 1 ))
NUMMIN=`awk '($3=="minimas") {print $5}' summary_${i}.dat`
NUMMAX=`awk '($3=="maximas") {print $5}' summary_${i}.dat`
for ((cc=1;cc<=${NUMMIN};cc++))
do
DGSTARSTOR[${COUNTER}]=`awk '(($2=="Minima") && ($4=="Maximum")) {print $7}' summary_${i}.dat | head -n 1 | tail -n 1`
#echo "${COUNTER}	${cc}"
#echo "${DGSTARSTOR[${COUNTER},${cc}]}"
done
for ((cc=1;cc<=${NUMMAX};cc++))
do
DGFREESTOR[${COUNTER}]=`awk '(($2=="Minimum") && ($4=="Minimum")) {print $7}' summary_${i}.dat | head -n 1 | tail -n 1`
#echo "${COUNTER}        ${cc}"
#echo "${DGFREESTOR[${COUNTER},${cc}]}"
done
if [[ ${LIST} != "empty" ]] ; then
mv energy_${i}.dat lambda_${i}.dat distance_${i}.dat ../
cd ..
fi
done
#Now comes the real part, generating the Gnuplot scripts after making sure that no old files are there
#Otherwise there might be problems with the new generated files...
if [[ -f *plot*.gplt ]];then
rm *plot*.gplt
fi
#This loop now generates the energy and lambda plots
#If you want to print the plots to pdf or ps, just change the terminal here from x11 to pdf/ps
for plots in energy lambda distance;do
if [[ ${plots} == "energy" ]];then
xlabel="energy gap [kcal/mol]"
ylabel[1]="dG - normalised [kcal/mol]"
ylabel[2]="Distance of reaction coordinate [Å]"
elif [[ ${plots} == "lambda" ]];then
xlabel="Lambda"
ylabel[1]="energy gap [kcal/mol]"
ylabel[2]="Points per bins"
else
xlabel="Point"
ylabel[1]="Distance"
ylabel[2]="Distance"
fi
if [[ ${plots} != "distance" ]] ; then
gnuplot_static ${plots}plot.gplt x11 2
gnuplot_plot ${plots}plot.gplt '0.0' '0.0' "${xlabel}" "${ylabel[1]}" '2' ${plots}
gnuplot_plot ${plots}plot.gplt '0.0' '0.5' "${xlabel}" "${ylabel[2]}" '3' ${plots}
else
gnuplot_static ${plots}plot.gplt x11 1
gnuplot_plot_distance ${plots}plot.gplt '0.0' '0.0' "${xlabel}" "${ylabel[1]}" ${plots}
fi

gnuplot -persist <${plots}plot.gplt &

done

make_energy_table




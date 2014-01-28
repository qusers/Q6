#!/bin/bash -l
#Paul Bauer, 2014-01-21
#Script to evaluate a Q EVB run and obtain reaction free energies, activation energies
#Distances of reacting atoms and an evaluation of the sampling rate
#Released under the Beer-Ware License

#Path to Qfep binary, needed in any case
qfep='Qfep5'

#This function makes the static inputs that are always the same for all Gnuplot scripts
#So we only need to write the stuff once and not a number of times
#It takes the filename of the script as input to know where to write to and the terminal type
function gnuplot_static() {
echo "set terminal $2 enhanced size 1920,1200
set size 1,1
set origin 0.0, 0.0
set multiplot layout 2,1
set size 1,0.5
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
if [[ ${7} == "lambda" ]];then
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
	echo "plot '${7}_${i}.dat' using 1:$6 title '${i}'  ${LINES} ls ${START} ">>$1
else
	if [[ ${LINEC} -eq 1 ]];then
		echo "plot '${7}_${i}.dat' using 1:$6 title '${i}'  ${LINES} ls ${START} ,\\">>$1
	elif [[ ${LINEC} -eq ${COUNTER} ]];then
		echo " '${7}_${i}.dat' using 1:$6 title '${i}'  ${LINES} ls ${START}">>$1
	fi
	LINEC=$(( ${LINEC} + 1 ))
fi
done
}
#You can iterate over directories and other stuff, this is just one local test so no loop is executed
#You will need the two awk scripts and a mapping input in the current directory
LIST="empty"
#To see how many files are plotted, needed to generate Gnuplot input later on
COUNTER=0
for i in ${LIST} ;do
cp get*.awk map_files.inp ${i}
cd ${i}
${qfep} <map_files.inp>map_files.out
awk -f getenergy.awk <map_files.out>energy_${i}.dat
awk -f getlambda.awk <map_files.out>lambda_${i}.dat
#Moves the files back for plotting everything together
mv energy_${i}.dat lambda_${i}.dat ../

cd ..
COUNTER=$(( ${COUNTER} + 1 ))

done
#Now comes the real part, generating the Gnuplot scripts after making sure that no old files are there
#Otherwise there might be problems with the new generated files...
if [[ -f *plot*.gplt ]];then
rm *plot*.gplt
fi
#This loop now generates the energy and lambda plots
#If you want to print the plots to pdf or ps, just change the terminal here from x11 to pdf/ps
for plots in energy lambda;do
if [[ ${plots} == "energy" ]];then
xlabel="energy gap [kcal/mol]"
ylabel[1]="dG - normalised [kcal/mol]"
ylabel[2]="Distance of forming bond [Å]"
else
xlabel="Lambda"
ylabel[1]="energy gap [kcal/mol]"
ylabel[2]="Points per bins"
fi
gnuplot_static ${plots}plot.gplt x11
gnuplot_plot ${plots}plot.gplt '0.0' '0.0' "${xlabel}" "${ylabel[1]}" '2' ${plots}
gnuplot_plot ${plots}plot.gplt '0.0' '0.5' "${xlabel}" "${ylabel[2]}" '3' ${plots}

gnuplot -persist <${plots}plot.gplt &

done


#!/bin/sh

#set directories
export RESULTS_DIR=./results
export CEL_DIR=./CEL
export BIN_DIR=./bin
export SRC_DIR=./src

#Gresham added directories
export PEAKS_DIR=./PEAKS
export SGR_DIR=./SGR
export IGB_DIR=./IGB
export NORMALIZED_DIR=./NORMALIZED

norm=$3
strand=$4

if [ $# -lt 2 ]; then
	echo "./tiling_predict [ cel file ] Train option [yes/no] Norm option [SI(default)/RMA1/RMA2] Strand option [reverse(default)/forward";
	exit 0
fi

if [ $# -lt 3 ]; then
	norm=SI
fi

if [ $# -lt 4 ]; then
	strand=reverse
fi


if [ $norm = "RMA2" ]; then
	echo "RMA2 requires 5 FY hybridizations used in training (FY3-2.CEL, FY3-3.CEL, FY3-4.CEL, FY3-5.CEL, FY3-6.CEL )!!!";
fi

#delete previous results files
rm results_chr*
#rm full_results_chr*

#take name of Cel file file as input
echo $1

#compile SNPScanner
g++ -O -m64 -o tiling $SRC_DIR/SNPScanner.cpp

#run default prediction
./tiling $1 $2 $norm $strand

#concatenate and write full results file to results directory
results=$RESULTS_DIR/results_${1}
full_results=full_results_${1}

cat results_chr1 results_chr2 results_chr3 results_chr4 results_chr5 results_chr6 results_chr7 results_chr8 results_chr9 results_chr10 results_chr11 results_chr12 results_chr13 results_chr14 results_chr15 results_chr16 > $results

#take name of Cel file file as input
echo "producing summary files for $1 results...........";

#find peaks compared to default/model using standard criteria and write to the PEAKS directory and then sort them according to score

echo "Attempting to identify peaks and write file to PEAKS directory.....";

perl ./bin/compare_strains.pl $RESULTS_DIR/results_$1 default 5 1 yes 1 7 > $PEAKS_DIR/results_$1_peaks5_1_yes_1_7

echo "Attempting to sort peaks file...";

#sort +2nr $PEAKS_DIR/results_$1_peaks5_1_yes_1_7 > $PEAKS_DIR/results_$1_peaks5_1_yes_1_7_sorted
sort $PEAKS_DIR/results_$1_peaks5_1_yes_1_7 > $PEAKS_DIR/results_$1_peaks5_1_yes_1_7_sorted

echo "number of peaks found is:";

echo wc $PEAKS_DIR/results_$1_peaks5_1_yes_1_7

#rename the normalized intensity file <STRAIN_NORM> to match the sample name and write to the NORMALIZED directory

echo "renaming STRAIN_NORM according to sample name..........";

mv STRAIN_NORM $NORMALIZED_DIR/$1_NormInt

#make a .sgr file from the intensity file and write to the SGR directory                                                
echo "making sgr file from intensity file..............";

perl ./bin/makeSgrFromAffyMerged.pl $NORMALIZED_DIR/$1_NormInt

#make a .sgr file from the results file and write to the SGR diectory

echo "making sgr file from results file........."; 

perl ./bin/makeSGR.pl $RESULTS_DIR/results_$1

#map the snp calls in the peak file to the genome features and write to the PEAKS directory

#echo "mapping calls to genome features..........";

perl ./bin/snpScannerSummary.pl $PEAKS_DIR/results_$1_peaks5_1_yes_1_7_sorted 

#make an IGB bookmarks file and write to the IGB directory 

#echo "making an IGB bookmarks file..............";

#perl ./bin/igbBookmarks.pl $PEAKS_DIR/results_$1_peaks5_1_yes_1_7_sorted

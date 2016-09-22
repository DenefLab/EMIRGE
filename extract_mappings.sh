#!/bin/bash
set -e 
for i in */ ; do
	cd $i
	# print working directory to the screen
	echo $i
	# Check if /iter.65 exists
	if [ -d "iter.65" ]; then
		# Extract mapped number of single reads
		samtools flagstat iter.65/bowtie.iter.65.PE.bam | sed -n 7p | awk -v pwd="${PWD##*/}" '{print pwd, $1}' >> /scratch/vdenef_fluxm/rprops/Emirge/processing/mapped.reads.txt
		# Extract the sequences names per sample and store in long format .txt file
		grep ">" iter.65/iter.65.cons.renamed.fasta | awk -v pwd="${PWD##*/}" '{print pwd, $1, $2, $3, $4}' >> /scratch/vdenef_fluxm/rprops/Emirge/processing/read.info.txt
		# go to previous work directory
	fi
        cd - 
done
# Remove all 'Prior=', 'NormPrior' and 'Length=' labels from the file
sed -i -e 's/NormPrior=//g' /scratch/vdenef_fluxm/rprops/Emirge/processing/read.info.txt 
sed -i -e 's/Prior=//g' /scratch/vdenef_fluxm/rprops/Emirge/processing/read.info.txt 
sed -i -e 's/Length=//g' /scratch/vdenef_fluxm/rprops/Emirge/processing/read.info.txt
sed -i -e 's/>//g' /scratch/vdenef_fluxm/rprops/Emirge/processing/read.info.txt

#!/bin/bash
set -e
for i in */ ; do
        cd $i
		# print working directory to the screen
        echo $i
		# Check if /iter.65 exists
		if [ -d "iter.65" ]; then
		# make new fasta file which concatenates all the individual files
		cat iter.65/iter.65.cons.renamed.fasta >> /scratch/vdenef_fluxm/rprops/Emirge/processing/total.emirge.renamed.fasta
		# go to previous work directory
fi
cd -
done

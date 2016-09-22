#!/bin/bash
set -e
for i in */ ; do
        cd $i
		# print working directory to the screen
        echo $i
		# Check if /iter.65 exists
		if [ -d "iter.65" ]; then
		# Run rename script to included prior abundances and rank sequences according to abundance
		emirge_rename_fasta.py iter.65	> iter.65/iter.65.cons.renamed.fasta
		fi
		# go to previous work directory
        cd -
done

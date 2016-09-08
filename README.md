# EMIRGE
Set of scripts/guidelines for EMIRGE analysis of metagenomic data.

## The following script has been tested with the following dependencies

```R
bowtie v1.0.0
usearch v8.1
python-anaconda2/201607
emirge v0.60.3
mothur v1.38.1 
biopython v1.60
ncbi-blast v2.2.29
```

## Workflow
Copy metaG data from nfs to scratch. This code copies only the fastq.gz in from the respective folders to the new location. 
```R 
rsync -a --include '*/' --include '*.fastq.gz' --exclude '*' /nfs/vdenef-lab/Shared/Sequence_data/CSP_LM13/LM13_JGI_MetaG /scratch/vdenef_fluxm/rprops/metaG --progress
```

Unzip all .gz fastq in the folders. -r stands for recursive, i.e. enter all directories in the provided path and unzip stuff. Best to run this separately in a pbs script if possible, because it automatically parallelizes this over multiple cores, or run it as a background process with <code>nohup</code>.
```R
gunzip -r /scratch/vdenef_fluxm/rprops/metaG/
nohup gunzip -r /scratch/vdenef_fluxm/rprops/metaG/ &
```
Extract forward and reverse fastq from original fasta (script comb_to_rever_forw_fastq.pbs).
```R
qsub comb_to_rever_forw_fastq.pbs
```
Copy the fastq of your samples to the correct directory (example for one sample).

```R
cp /nfs/vdenef-lab/Shared/Sequence_data/CSP_LM13/LM13_JGI_MetaG/Fa13.BD.MM110.DN/Fa13.BD.MM110.DN.10Mpairs.PE.2.fastq /scratch/vdenef_fluxm/rprops/metaG

cp /nfs/vdenef-lab/Shared/Sequence_data/CSP_LM13/LM13_JGI_MetaG/Fa13.BD.MM110.DN/Fa13.BD.MM110.DN.10Mpairs.PE.1.fastq /scratch/vdenef_fluxm/rprops/metaG
```

Cluster silva_v123 NR99 database to 97 % with usearch. In the original article they also prune the database from sequences with a length less than 1200 and more than 1900.

```R
usearch -cluster_fast SILVA_123_SSURef_Nr99_tax_silva2.fasta -id 0.97 -centroids SSURef_NR97_123_for_emirge.fasta -uc SSURef_NR97_123_for_emirge.clusters.uc
```

Replace non-standard base characters in reference database (script in same directory as database). This was run with python-anaconda2/latest - could also work with the 201607 but not sure yet

```R
python fix_nonstandard_chars.py < SSURef_NR97_123_for_emirge.fasta > /scratch/vdenef_fluxm/rprops/databases/SSURef_NR97_123_for_emirge2.fasta</code>
```

Create bowtie-index for your reference database
```R
bowtie-build SSURef_NR97_123_for_emirge2.fasta SSU_candidate_db_btindex
```
Make sure these scripts are in the directory where all your sample folders are. Run EMIRGE (see batch scripts for parallelizing). Make sure that --Phred33 is present (this is mandatory for new fastq files from Illumina). Takes approx. 20-40h on 10 cores per sample for 65 iterations. Running with -j 1.0 (thus 100 % identity) which is different from the 0.97 that others use (unique seqs). There is also little point in parallelizing beyond 10 - 12 cores because of the dependency on usearch (uses 8 cores)
```R
bash -x Batch_01.sh
bash -x Batch_02.sh
bash -x Batch_03.sh
bash -x Batch_04.sh
bash -x Batch_05.sh
bash -x Batch_06.sh
```
## Process EMIRGE output
Reformat fasta header to sort sequences according to abundances. Usage: emirge_rename_fasta.py [options] <iter.DIR> > renamed.fasta
Make sure biopython/1.60 is loaded. Make sure that the files mapped.reads.txt, total.emirge.cons.fasta and reads.info.txt don't exist yet.
```R
bash -x rename.sh
```

Count number of mapping reads from bowtie bam files in final iteration folders. This is needed to estimate the number of reads for each sequence

```R
bash -x extract_mappings.sh
```

Concatenate all the fasta files in one fasta for taxonomic classification

```R
bash -x concat.emirge.sh
```

## Classify the sequences according to the TaxASS pipeline for freshwater communities (<a href="https://github.com/McMahonLab/TaxAss">https://github.com/McMahonLab/TaxAss</a> )

### Step 1: Remove all redundant information from the first line in the fasta file. Retain sample names only. Fasta file must not be aligned.

```R
awk '{print $1}' total.emirge.renamed.fasta > total.emirge.renamed2.fasta
mv total.emirge.renamed2.fasta total.emirge.renamed.fasta
```
### Step 2: make BLAST database file (blast)

```R
makeblastdb -dbtype nucl -in FreshTrain18Aug2016.fasta -input_type fasta -parse_seqids -out FWonly_18Aug2016custom.db
```

### Step 3: Run blast and reformat output blast file
```R
blastn -query total.emirge.renamed.fasta -task megablast -db /scratch/vdenef_fluxm/rprops/Emirge/BLAST/FWonly_18Aug2016custom.db -out custom.blast -outfmt 11 -max_target_seqs 5

blast_formatter -archive custom.blast -outfmt "6 qseqid pident length qlen qstart qend" -out otus.custom.blast.table
```
### Step 4: Correct BLAST pident using custom script
This accounts for sequence length differences

```R
Rscript calc_full_length_pident.R otus.custom.blast.table otus.custom.blast.table.modified
```

###Step 5: Filter BLAST results

```R
Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.above.97 97 TRUE

Rscript filter_seqIDs_by_pident.R otus.custom.blast.table.modified ids.below.97 97 FALSE
```
###Step 6: Make plots to evaluate blast run
```R
mkdir plots

Rscript plot_blast_hit_stats.R otus.custom.blast.table.modified 97 plots
```
###Step 7: recover sequence IDs left out of blast (python, bash)

```R
python find_seqIDs_blast_removed.py total.emirge.renamed.fasta otus.custom.blast.table.modified ids.missing
cat ids.below.97 ids.missing > ids.below.97.all
```
###Step 8: create fasta files of desired sequence IDs (python)

```R
python create_fastas_given_seqIDs.py ids.above.97 total.emirge.renamed.fasta otus.above.97.fasta

python create_fastas_given_seqIDs.py ids.below.97.all total.emirge.renamed.fasta otus.below.97.fasta
```

###Step 9: extract unique sequences
```R
mothur "#unique.seqs(fasta=otus.below.97.fasta)"

mothur "#unique.seqs(fasta=otus.above.97.fasta)"
```
###Step 10: classify sequences
```R
mothur "#classify.seqs(fasta=otus.below.97.unique.fasta, template=/scratch/vdenef_fluxm/rprops/databases/silva.nr_v123.align, taxonomy=/scratch/vdenef_fluxm/rprops/databases/silva.nr_v123.tax, method=wang, probs=T, processors=10, cutoff=80)"

mothur "#classify.seqs(fasta=otus.above.97.unique.fasta,template=/scratch/vdenef_fluxm/rprops/databases/FreshTrain18Aug2016.fasta,  taxonomy=/scratch/vdenef_fluxm/rprops/databases/FreshTrain18Aug2016.taxonomy, method=wang, probs=T, processors=10, cutoff=0)"
```
Step 11: combine taxonomy files and names files
```R
cat otus.above.97.unique.FreshTrain18Aug2016.wang.taxonomy otus.below.97.unique.nr_v123.wang.taxonomy > otus.97.taxonomy

cat otus.above.97.names otus.below.97.names > otus.97.names
```
Step 12: Run R script to create Sequence table (same format as OTU table)
```R
Rscript emirge_format.R mapped.reads.txt otus.97.taxonomy read.info.txt otus.97.names
```

library('splitstackshape')
library('dplyr')
library('tidyr')
library("phyloseq")
# library('qdap')
library('easyGgplot2')

### Get working directory and read in files
userprefs <- commandArgs(trailingOnly = TRUE)
mapped.reads <- userprefs[1]
taxonomy <- userprefs[2]
read.info <- userprefs[3]
names <- userprefs[4]

### Data import and formatting
sample_size <- read.table(file=mapped.reads, col.names=c("Sample","Mapped.reads"))
tax.table <- read.table(file=taxonomy,stringsAsFactors = FALSE)
tax.table <- data.frame(cSplit(indt = tax.table, splitCols = "V2", sep = ";", drop = TRUE))
# tax.table <- genX(tax.table, "(", ")")
rownames(tax.table) <- tax.table[,1]; tax.table <- tax.table[,-1]; colnames(tax.table) <- c("kingdom","phylum","class","order","lineage","clade","tribe")
tax.table <- tax_table(as.matrix(tax.table))
tax.table[is.na(tax.table)] <- "Unclassified"
count_seq <- read.table(file=read.info, col.names=c("Sample","Sequence","Prior","Length","NormPrior"),stringsAsFactors = FALSE)
name_seq <- read.table(file=names,stringsAsFactors = FALSE)

### Remove all sequences (count_seq) which are not in the name file
name_seq_temp <- data.frame(cSplit(indt = name_seq, splitCols = "V2", sep = ",", drop = TRUE,type.convert = FALSE))
#remove NA and find uniques
uniques <- data.frame(Sequence=as.character(na.omit(unique(unlist(name_seq_temp)))),stringsAsFactors = FALSE)
#remove sequences from count
count_seq <- semi_join(x=count_seq,y=uniques,by="Sequence")

### Find sequences that are identical
name_seq <- name_seq[!(name_seq[,1]==name_seq[,2]),]
name_seq <- data.frame(cSplit(indt = name_seq, splitCols = "V2", sep = ",", drop = TRUE,type.convert = FALSE))

### Format to two columns
init <- nrow(name_seq)
name_seq <- as.matrix(name_seq)
name_seq_final <- as.matrix(name_seq[,1:2])
for(i in 2:ncol(name_seq)){
  name_seq_final <- rbind(name_seq_final[,1:2],cbind(name_seq[(1:init),c(1,i)]))
}
### remove NA rows
name_seq_final <- data.frame(name_seq_final[!is.na(name_seq_final[,2]),],stringsAsFactors=FALSE); colnames(name_seq_final) <- c("unique","replicate")

### Loop through second column and if sequence name is found, change name to the unique one in first column
for(i in 1:nrow(name_seq_final)){
  count_seq$Sequence[which(count_seq$Sequence == name_seq_final$replicate[i])] <- name_seq_final$unique[i]
}

count_seq <- cbind(count_seq,counts.prior=rep(1,nrow(count_seq)),counts.norm=rep(1,nrow(count_seq)))

### Multiply priors by mapped reads
for(i in 1:nrow(count_seq)){
  count_seq$counts.prior[which(count_seq$Sample == sample_size$Sample[i])] <- floor(count_seq$Prior[which(count_seq$Sample == sample_size$Sample[i])] * sample_size$Mapped.reads[i])
  count_seq$counts.norm[which(count_seq$Sample == sample_size$Sample[i])] <- floor(count_seq$NormPrior[which(count_seq$Sample == sample_size$Sample[i])] * sample_size$Mapped.reads[i])
}
count_seq$Sample <- factor(count_seq$Sample)

### Make plot of length distribution
png("length_distribution.png",width=6,height=5,res=500,units="in",pointsize=12)
ggplot2.histogram(data=count_seq, xName='Length',
                  fill="white", color="black",
                  linetype="longdash",binwidth=round(nrow(count_seq)/500,0),addMeanLine=TRUE, meanLineColor="red",
                  meanLineType="dashed", meanLineSize=1)+
  theme_bw() + labs(y="Frequency")
dev.off()

### Final sequence table requires reformatting from long to wide format
data_wide <- spread(count_seq[,c(1,2,6)],Sequence,counts.prior)
data_wide[is.na(data_wide)] <- 0
data_wide[,1] <- as.character(data_wide[,1])
rownames(data_wide) <- data_wide[,1]; data_wide <- data_wide[,-1]
data_wide <- data.frame(t(data_wide))

# Output sequence table
write.table(x=data_wide,file="otu_table.txt")
###################################################################################################################


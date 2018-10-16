---
title: "R lesson"
output: html_document
# Data Load from my desktop
Read data fang_et_at_genotypes AND SNP_position files

List the packages that need to be installed at the start of the README. JZ

Genotype<-read.delim(file.choose(), header=T)
SNPs<-read.delim(file.choose(),header = T)
fang_et_al<read.delim(file.choose(),header =T)

I am not sure what to do to get started. I have tried putting the Fang and SNP data into my working directory but it has not worked. JZ
#DATA INSPECTION:
head() to inspect the files. 
head(Genotypes, 2, n=5)
head(SNP_position, n=5)
Both files are structured as a list
typeof(Genotype)
typeof(SNPs)

Equally we have that both files are a data.frame according to the function class. For Genotype file we have all variables as factors, while for SNP_position we have some variables as a integer.

class(Genotype)
class(SNPs)
str(Genotype)
str(SNPs)

The Genotype file has 2782 number of rows and 986 number of columns, while the SNP_position has 983 rows with 15 columns.
dim(Genotype)
dim(SNPs)

The name of each column per file is given by colnames function.
colnames(Genotype)
colnames(SNPs)

Use the levels() function to investigate the levels.
levels(SNPs$Chromosome)
levels(Genotype$Group)

All of these commands worked for me. JZ

### DATA PROCESSING

#SPLIT DATA INTO GROUPS
maize <- subset(Genotypes, Group == "ZMMIL" | Group == "ZMMLR" | Group =="ZMMMR")

teosinte <- subset(Genotypes, Group == "ZMPBA" | Group == "ZMPIL" | Group == "ZMPJA")

tripsacum <- subset(Genotypes, !Group == "ZMPBA" & !Group == "ZMPIL" & !Group == "ZMPJA" & !Group == "ZMMIL" & !Group == "ZMMLR" & !Group=="ZMMMR") The results shows the total number of lines in the three files is the same as the totla number in the original file.

Spelling error in total above. I am having trouble with these but it should work. JZ

####Data transposing
The t function can be used to transpose data
transposed_maize <- as.data.frame(t(maize))

transposed_teosinte <- as.data.frame(t(teosinte))

The new dataframe do not include the header in the first column. To solve this problem, the function library(tibble) is loaded.

if(!require("tibble")) install.packages("tibble")
library(tibble)

names(transposed_maize) <- lapply(transposed_maize[1, ], as.character)
transposed_maize <- transposed_maize[-1,]

names(transposed_teosinte) <- lapply(transposed_teosinte[1, ], as.character)
transposed_teosinte <- transposed_teosinte[-1,]

transposed_maize <- rownames_to_column(transposed_maize, var="SNP_ID")
transposed_teosinte <- rownames_to_column(transposed_teosinte, var="SNP_ID")

In the data.frame snp_position the column 3 is for chromosome, the position of the chromosome is the column 4, while the column 1 (SNP_ID) is the common column with the transposed files. Then, we proceed to create a new file with only the information required from the snp_position.txt file:
snp_ID_chro_pos <- SNP_position[, c("SNP_ID", "Chromosome", "Position"

Once we have the new data.frame we proceed to order the data.frames of the transposed files as follows:
if(!require("dplyr")) install.packages("dplyr")
library(dplyr)

transposed_teosinte <- arrange(transposed_teosinte, SNP_ID) 
transposed_maize <- arrange(transposed_maize, SNP_ID) 
snp_ID_chro_pos <- arrange(snp_ID_chro_pos, SNP_ID) 

The files are sorted by snp_ID, therefore the SNP POSITION and chromosome can be joined to it.

teosinte_all_join <- merge(snp_ID_chro_pos,transposed_teosinte, by.x="SNP_ID", by.y="SNP_ID", all = TRUE)

maize_all_join <- merge(snp_ID_chro_pos,transposed_maize, by.x="SNP_ID", by.y="SNP_ID", all = TRUE)

The two rows that are not informative can now be removed.
teosinte_all_join <- subset(teosinte_all_join, !(teosinte_all_join$Chromosome == "unknown") | !(teosinte_all_join$Chromosome =="multiple"))

maize_all_join <- subset(maize_all_join, !(maize_all_join$Chromosome == "unknown") | !(maize_all_join$Chromosome == "multiple"))

After the two uninformative columns are removed,each dataframe is sorted out by position in an ascending order. Different from UNIX here we do not have to have in order the chromosomes to be able to extract the information. Sorting by position will give us after the split a new files with the position sorted in an ascending order.

teosinte_all_join <- arrange(teosinte_all_join, Position)
maize_all_join <- arrange(maize_all_join, Position)

# creating a loop to create 10 files per chromosome
teosinte_split <- split(teosinte_all_join, teosinte_all_join$Chromosome)
nom <- names(teosinte_split)
for(myname in nom){
  savename = paste0('Teosinte_Chr_', myname, '.txt')
write.table(teosinte_split[[myname]], file = savename, quote = FALSE, sep = "\t", row.names = FALSE)
}

maize_split <- split(maize_all_join, maize_all_join$Chromosome)
nom <- names(maize_split)
for(myname in nom){
  savename = paste0('Maize_Chr_', myname, '.txt')
write.table(maize_split[[myname]], file = savename, quote = FALSE, sep = "\t", row.names = FALSE)
}

Now, the default has the  missing command ? and has to be changed.

Good loop commands JZ. 

# Replacing the missing code ? with - . 
To do this, the data frame has to be ordered in a decreasing order and report that the values in the newd data fram are now characters in R.

rev_teosinte_all_join <- arrange(teosinte_all_join, desc(Position))
rev_maize_all_join <- arrange(maize_all_join, desc(Position))

rev_teosinte_all_join[] <- lapply(rev_teosinte_all_join, as.character)
rev_teosinte_all_join[rev_teosinte_all_join == '?/?'] <- '-/-'

rev_maize_all_join[] <- lapply(rev_maize_all_join, as.character)
rev_maize_all_join[rev_maize_all_join == '?/?'] <- '-/-'

After this is done, files per chromosome can be created.

rev_teosinte_split <- split(rev_teosinte_all_join, rev_teosinte_all_join$Chromosome)
nom <- names(rev_teosinte_split)
for(myname in nom){
  savename = paste0('Rev_teosinte_Chr_', myname, '.txt')
write.table(rev_teosinte_split[[myname]], file = savename, quote = FALSE, sep = "\t", row.names = FALSE)}

rev_maize_split <- split(rev_maize_all_join, rev_maize_all_join$Chromosome)
nom <- names(rev_maize_split)
for(myname in nom){
  savename = paste0('Rev_Maize_Chr_', myname, '.txt')
write.table(rev_maize_split[[myname]], file = savename, quote = FALSE, sep = "\t", row.names = FALSE)}

### PART 2
library(dplyr)

Transpose and merge original data
fang_et_al<-as.data.frame(t(fang_et_al))
joined.data<-merge(SNPs, fang_et_al, by.x="SNP_ID",by.y="row.names",all=TRUE)

# plot of SNPs per chromosome
library(ggplot2)
joined.data$Chromosome<-factor(joined.data$Chromosome, levels = c("1","2","3","4","5","6","7","8","9","10","unknown","multiple","NA"))
ggplot(joined.data)+ geom_bar(aes(joined.data$Chromosome))+xlab("Chromosome") +ylab("Total Number of SNPs")

I got rid of the NA JZ. 

# Tidying the data
library(reshape2)
genotype.info <- colnames(fang_et_al)[-c(1:3)]
fang_et_al.tidy<-melt(fang_et_al,measure.vars = genotype.info)

# Recording missing data as NA
fang_et_al.tidy[]<- lapply(fang_el_al.tidy, as.character)
fang_et_al.tidy[fang_et_al.tidy =='?/?'] <- 'NA'

# Classifying Genotype SNPs as homozygotes or heterozygotes. Dataframe sorted based on Group and Species_ID
library(plyr)
fang_et_al.tidy$hom.het <- (fang_et_al.tidy$value=="A/A"|fang_et_al.tidy$value=="C/C"|fang_et_al.tidy$value=="G/G"|fang_et_al.tidy$value=="T/T")
fang_et_al.class.sorted<-arrange(fang_et_al.tidy,Sample_ID,Group)
counts <- ddply(fang_et_al.class.sorted,c("Sample_ID"),summarise,total_homozygous=sum(hom.het,na.rm=TRUE),total_heterozygous=sum(!hom.het,na.rm = TRUE), total_NA=sum(is.na(hom.het)))
counts.combined<-melt(counts,measure.vars = c("total_homozygous","total_heterozygous","total_NA"))

#Ploting proportion of homozygous, heterozygous and missing sites per Group
counts.group<-ddply(fang_et_al.class.sorted,c("Group"),summarise,total_homozygous=sum(hom.het,na.rm=TRUE),total_heterozygous=sum(!hom.het,na.rm = TRUE), total_NA=sum(is.na(hom.het)))
counts.group.combined<-melt(counts.group,measure.vars = c("total_homozygous","total_heterozygous","total_NA"))
ggplot(counts.group.combined,aes(x=Group,y=value,fill=variable))+geom_bar(stat="identity",position = "stack")

### My own visulaization
maize.hom.het <- subset(fang_et_al, Group == "ZMMIL" | Group == "ZMMLR" | Group=="ZMMMR")
maize.hom.het.transposed <- as.data.frame(t(maize.hom.het))
joined.maize<-merge(SNPs, maize.hom.het.transposed, by.x="SNP_ID",by.y="row.names",all=TRUE)
maize.info <- colnames(maize.hom.het)[-c(1:3)]
maize.tidy<-melt(maize.hom.het,measure.vars = maize.info)
maize.tidy[]<- lapply(maize.tidy, as.character)
maize.tidy[maize.tidy=='?/?'] <- 'NA'
maize.tidy$hom.het <- (maize.tidy$value=="A/A"|maize.tidy$value=="C/C"|maize.tidy$value=="G/G"|maize.tidy$value=="T/T")
maize.class.sorted<-arrange(maize.tidy,Sample_ID,Group)
counts.maize.group<-ddply(maize.class.sorted,c("Group"),summarise,total_homozygous=sum(hom.het,na.rm=TRUE),total_heterozygous=sum(!hom.het,na.rm = TRUE), total_NA=sum(is.na(hom.het)))
maize.counts.group.combined<-melt(counts.maize.group,measure.vars = c("total_homozygous","total_heterozygous","total_NA"))
ggplot(maize.counts.group.combined,aes(x=Group,y=value,fill=variable))+geom_bar(stat="identity",position = "stack")

# My own viualization of common types of nucleotides across sites.
ggplot(fang_et_al.tidy)+ geom_bar(aes(fang_et_al.tidy$value))+xlab("Nucleotide") +ylab("Number of observations")
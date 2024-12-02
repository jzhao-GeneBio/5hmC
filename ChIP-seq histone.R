####################################################################
########### part I.1: find histone peaks overlapped with DhMRs
# 'H3K9acDomains.csv' downloaded from https://www.synapse.org/Synapse:syn17016212
histone_name=read.csv('H3K9acDomains.csv',sep='\t') 
dim(histone_name) # 26384
histone_name=histone_name[!histone_name[,'Chr']%in%c('X','Y'),]
histone_name=histone_name[,c('Chr','Start','End','DomainID')]

### convert hg19 genomic coordinate to hg38 version
#https://genome.ucsc.edu/cgi-bin/hgLiftOver
position38=read.delim('H3K9acDomains_hg38.bed',header=F)
dim(position38) # 25732
sum(position38[,4]%in%histone_name[,4]) # 25732
length(unique(position38[,4])) # 25725
freq=data.frame(table(position38[,4]))
table(freq[,2])
histone_name=histone_name[histone_name[,4]%in%position38[,4],]
sum(position38[,4]==histone_name[,4]) # check whether peak id match
histone_name[,2]=position38[,2]
histone_name[,3]=position38[,3]
#https://bioconductor.org/packages/devel/workflows/vignettes/liftOver/inst/doc/liftov.html
#library(rtracklayer)
### convert hg19 genomic coordinate to hg38 version
write.table(histone_name,'peaks_histone.txt',row.names=F,col.names=F,sep="\t",quote = FALSE)




####### part I.2: run 'bedtools.sbatch' (output: overlap_histone.txt)
#!/bin/sh

module load bedtools
bedtools intersect -wa -wb -a peaks_5hmC.txt -b peaks_histone.txt > overlap_histone.txt






####### part II: correlation between histone and 5hmC
peaks=read.delim("peaks_5hmC.txt",header=F)
match=read.delim('overlap_histone.txt',header=F)
match=match[match[,'5hmc']%in%peaks[,'5hmc'],]
dim(match)
match=match[,c('5hmc','histone')]

histone=read.table('~/ROSMAP/H3K9ac/readcounts_fpkm.txt', header = TRUE)
histone=histone[rownames(histone)%in%unique(match[,'histone']),]
dim(histone) 

hmC=read.table('5hmC_AD.txt', header = TRUE)
### match participant ID
id=read.delim('~/ROSMAP/5hmcData/Hydro_dataset_655_basic_10-11-2020.Neat1060.txt')
sum(rownames(hmC)%in%id[,'projid']) 
id=id[id[,'projid']%in%rownames(hmC),]
sum(id_histone%in%id[,'projid']) 
id_histone[!id_histone%in%as.numeric(id[,'projid'])]
id=id[id[,'projid']%in%id_histone,]
dim(id)
hmC=hmC[rownames(hmC)%in%id[,'id'],]
dim(hmC)
hmC=hmC[match(id[,'id'],rownames(hmC)),]
histone=histone[id_histone%in%as.numeric(id[,'projid']),]
id_histone=as.numeric(rownames(histone))
histone=histone[match(as.numeric(id[,'projid']),id_histone),]
dim(histone)

corr=rep(0,nrow(match))
p=rep(0,nrow(match))
for (i in 1:nrow(match)){
  corr[i]=cor(hmC[,colnames(hmC)==match[i,'5hmc']],histone[,colnames(histone)==match[i,'histone']],method='spearman')
  p[i]=cor.test(hmC[,colnames(hmC)==match[i,'5hmc']],histone[,colnames(histone)==match[i,'histone']],method='spearman')$p.value
  print(i)
}
match=cbind(match,corr,p)
match=data.frame(match)
dim(match)
write.csv(match,'corr_AD.csv')

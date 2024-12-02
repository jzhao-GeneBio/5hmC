library("data.table")
library("readr") 
file='~/ROSMAP/5hmcData/5hmC data.txt'
# read the first row
col_name <- scan(file,what = "",sep = "\t",nlines = 1,skip = 0,quiet = TRUE,strip.white = TRUE)
length(col_name) 
dt_gene <- read_delim(file,skip=1,col_names = F)
dim(dt_gene)  
length(unique(as.matrix(dt_gene[,1]))) #  197765
colnames(dt_gene)=c('peakID',col_name)
ID=colnames(dt_gene)[2:ncol(dt_gene)] # col_name

dim(dt_gene) 
dt_gene=as.matrix(dt_gene) 
dt_gene=t(dt_gene)
dim(dt_gene)
colname=dt_gene[1,]
dt_gene=dt_gene[-1,]
colnames(dt_gene)=colname
dt_gene=cbind(ID,dt_gene)
colnames(dt_gene)[1]='projid'
rownames(dt_gene)=dt_gene[,1]
dt_gene=dt_gene[,-1]
id=rownames(dt_gene)
dt_gene=as.matrix(dt_gene)
dt_gene=as.numeric(dt_gene)
dim(dt_gene)=c(1050,197765)

####### PCA 
x=dt_gene # row: sample, column: 5hmC
dim(x)
library(factoextra)
pcs=prcomp(x,center = TRUE,  scale = TRUE)
var_pc=get_eig(pcs)[1:100,'variance.percent'] 
pcs=pcs$x
dim(pcs) 
rownames(pcs)=id
basic=pcs[,1:100]
write.table(basic,'PCs_100.txt',row.names=T,col.names=T,sep="\t",quote = FALSE)
write.table(var_pc,'PCs_100_VarExplained.txt',row.names=F,col.names=T,sep="\t",quote = FALSE)

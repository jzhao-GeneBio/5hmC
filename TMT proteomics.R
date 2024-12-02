####################################################################
########### part I: find proteins overlapped with DhMRs
expression=read.table("position_protein.txt", header = F)
dim(expression)
ad=read.delim('~/ROSMAP/program/regional plot/peaks_AD.txt',header=F)
# within 50kb of each 5hmC region
ad[,2]=as.numeric(ad[,2])-50*1000
ad[,3]=as.numeric(ad[,3])+50*1000
overlap=NULL
for (i in 1:nrow(ad)){
  expression_sub=expression[expression[,1]==ad[i,1],]
  ind_overlap=(expression_sub[,3]>ad[i,2])*(expression_sub[,2]<ad[i,3])
  n_overlap=sum(ind_overlap)
  if (n_overlap>0){
    overlap=rbind(overlap,cbind(ad[i,],expression_sub[ind_overlap==1,])) # have warning message when one 5hmC peak matches with multiple gene expression
  }
  print(i)
}
dim(overlap) 
write.table(overlap,'overlap_protein 5hmC.txt',row.names=F,col.names=F,sep="\t",quote = FALSE)




####################################################################
############# part II: match 5hmC with protein and calculate correlations
protein=read.table('protein.csv',sep=',',header=T)
pair=read.delim("overlap_protein 5hmC.txt",header=F)
dim(pair)
protein=protein[,colnames(protein)%in%pair[,'protein']]

hmC=read.table('5hmC_AD.protein.txt', header = TRUE)
hmC=hmC[,colnames(hmC)%in%pair[,'dhmr']]
dim(hmC)

hmC=hmC[as.numeric(rownames(hmC))%in%as.numeric(rownames(protein)),]
protein=protein[as.numeric(rownames(protein))%in%as.numeric(rownames(hmC)),]
hmC=hmC[match(as.numeric(rownames(protein)),as.numeric(rownames(hmC))),]

corr=rep(0,nrow(pair))
p=rep(0,nrow(pair))
for (i in 1:nrow(pair)){
  corr[i]=cor(as.numeric(hmC[,colnames(hmC)==pair[i,'dhmr']]),as.numeric(protein[,colnames(protein)==pair[i,'protein']]),method='spearman',use='pairwise.complete.obs')
  p[i]=cor.test(as.numeric(hmC[,colnames(hmC)==pair[i,'dhmr']]),as.numeric(protein[,colnames(protein)==pair[i,'protein']]),method='spearman',use='pairwise.complete.obs')$p.value
}
fdr=p.adjust(p,method="fdr")
print('significant after multiple testing')
sum(fdr<0.05)
match=cbind(pair[,c('dhmr','protein')],corr,p)
match=data.frame(match)
dim(match)
write.csv(match,'corr_AD.csv',row.names=F)

library(RUVSeq)
library(readxl)
basic <- read_excel("~/ROSMAP/dataset_1009_basic_01-30-2023.xlsx")
basic <- as.data.frame(basic) 
basic2 <- basic[,c("projid",'AD',"msex","age_death","pmi",'educ')]
rosmap_5hmc <- read.table("~/ROSMAP/5hmcData/5hmC_data.txt", header = TRUE)
x=basic2[,c('AD','msex','age_death','pmi','educ')] 
design <- model.matrix(~AD+factor(msex)+age_death+pmi+educ,data = x)
y=rosmap_5hmc # p*n dimensional 5hmC data
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
dim(y) # p*n
dim(design) 
head(design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")
controls <- rownames(rosmap_5hmc)
seqRUVr <- RUVr(rosmap_5hmc, controls, k=10, res)
UV=pData(seqRUVr)
dim(UV)
saveRDS(UV, "AD_unwanted variation.RDS")



############################ regress out unwanted variation ############################
a=read.delim('~/5hmcData/clinical data.txt')
a=data.frame(a)
### read 5hmC data 
hmC=readRDS('5hmC_AD.RDS')
outcome=merge(a[,c('projid','study','ad_reagan','age_death','msex','pmi','educ')],hmC,by.x="projid",by.y="projid")
### read unwanted variation 
UV=readRDS("AD_unwanted variation.RDS")
outcome=merge(outcome,UV,by.x="projid",by.y="projid")
n_peak=ncol(outcome)-7 # the first 7 columns in 'outcome' data are covariates
outcome[,'age_death']=as.numeric(outcome[,'age_death'])
outcome[,'educ']=as.numeric(outcome[,'educ'])
outcome[,'pmi']=as.numeric(outcome[,'pmi'])
for (i in 8:(7+n_peak)){
  outcome[,'hmc']=as.numeric(outcome[,i])
  fit=lm(hmc~UV1+UV2+UV3+UV4+UV5+UV6+UV7+UV8+UV9,outcome)
  outcome[,i]=residuals(fit)
}
write.table(outcome,'residual_AD.csv',sep=',',row.names = F, col.names=T)

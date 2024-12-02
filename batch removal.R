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
HF=pData(seqRUVr)
dim(HF)
saveRDS(HF, "AD_unwanted variation.RDS")

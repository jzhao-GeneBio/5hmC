################## part I: meta-analysis
output=read.csv("~/5hmC_AD/output/glm_map_AD.csv",header=F)
colnames(output)=c('hmc',"beta","se","t","p")

output.2=read.csv("~/5hmC_AD/output/glm_ros_AD.csv",header=T)
colnames(output.2)=c('hmc',"beta","se","t","p")

output.2=output.2[match(output[,1],output.2[,1]),]

### meta-analysis
library(meta) 
res.meta=matrix(data = NA, nrow = nrow(output), ncol = 10)
for (i in 1:nrow(output)){
  yi=c(output[i,'beta'],output.2[i,'beta']) 
  vi=c(output[i,'se'],output.2[i,'se']) 
  out=metagen(yi,vi)
  res.meta[i,]<-c(out$TE.fixed, out$seTE.fixed, out$pval.fixed)
  print(i)
}
colnames(res.meta)<-c("TE.fixed", "seTE.fixed", "pval.fixed")
new=data.frame(cbind(output[,1],res.meta))
dim(new) 
write.csv(new,file='~/5hmC_AD/output/meta_AD.csv', row.names=FALSE)

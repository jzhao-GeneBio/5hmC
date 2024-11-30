a=read.delim('~/5hmcData/Hydro_dataset_655_basic_10-11-2020.Neat1060.txt')
#colnames(a)
a=data.frame(a)
table(a[,'ad_reagan'])
### read 5hmC data 
hmC=readRDS('5hmC.RDS')

sum(hmC[,'projid']%in%a[,'projid']) 
outcome=merge.data.frame(a[,c('projid','study','amylsqrt','tangles','age_death','msex','pmi','educ')],hmC,by.x = "projid",by.y="projid")
dim(outcome)                           

outcome[,'age_death']=as.numeric(outcome[,'age_death'])
outcome[,'educ']=as.numeric(outcome[,'educ'])
outcome[,'pmi']=as.numeric(outcome[,'pmi'])
outcome[,'amylsqrt']=as.numeric(outcome[,'amylsqrt'])
outcome[,'tangles']=as.numeric(outcome[,'tangles'])

dt.dhmr=read.delim("DhMR.identified.txt")
dhmr=dt.dhmr[,'DhMR']
keep=c('age_death','msex','educ','pmi','amylsqrt','tangles',dhmr)
outcome=outcome[,colnames(outcome)%in%keep]
output=matrix(0,length(dhmr),2)

j=1
for (i in 7:ncol(outcome)){ 
  outcome[,'dhmr']=as.numeric(outcome[,i])
  fit0=lm(amylsqrt~age_death+msex+educ+pmi,outcome)
  fit=lm(amylsqrt~dhmr+age_death+msex+educ+pmi,outcome)
  output[j,1]=summary(fit)$r.squared-summary(fit0)$r.squared
  j=j+1
}

j=1
for (i in 7:ncol(outcome)){ 
  outcome[,'dhmr']=as.numeric(outcome[,i])
  fit0=lm(tangles~age_death+msex+educ+pmi,outcome)
  fit=lm(tangles~dhmr+age_death+msex+educ+pmi,outcome)
  output[j,2]=summary(fit)$r.squared-summary(fit0)$r.squared
  j=j+1
}
output=cbind(colnames(outcome)[7:ncol(outcome)],output)
write(t(output),'R2.csv',ncol=ncol(output),sep=',')

############################ AD diagnosis ############################
a=read.delim('~/5hmcData/Hydro_dataset_655_basic_10-11-2020.Neat1060.txt')
#colnames(a)
a=data.frame(a)
#head(a)
a=a[a[,'study']=='MAP',] # train: MAP; test: ROS
#a=a[a[,'study']=='ROS',]
table(a[,'ad_reagan'])
### read 5hmC data 
hmC=readRDS('5hmC_AD.RDS')
hmC[1:2,1:3]

sum(hmC[,'projid']%in%a[,'projid']) 
outcome=merge.data.frame(a[,c('projid','study','ad_reagan','age_death','msex','pmi','educ')],hmC,by.x = "projid",by.y="projid")
dim(outcome)                         
table(outcome[,'ad_reagan'])
table(outcome[,'msex'])
table(outcome[,'study'])
dim(outcome)

colnames(outcome)[ncol(outcome)]
n_peak=ncol(outcome)-7 # the first 7 columns in 'outcome' data are covariates
print('No. of peaks')
outcome[,'age_death']=as.numeric(outcome[,'age_death'])
outcome[,'educ']=as.numeric(outcome[,'educ'])
outcome[,'pmi']=as.numeric(outcome[,'pmi'])
output=matrix(0,n_peak,4)
for (i in 8:(7+n_peak)){
  outcome[,'peak']=as.numeric(outcome[,i])
  fit=glm(ad_reagan~peak+age_death+msex+educ+pmi,outcome,family=binomial())
  ctable <- coef(summary(fit))
  output[j,]=ctable[rownames(ctable)=='peak',]
  j=j+1
}
dim(output)
output=cbind(colnames(outcome)[start:end],output)
colnames(output)=c('peak',"beta","se","t","p")
write.table(output,'glm_map_AD.csv',sep=',',row.names = F, col.names=T)
#write.table(output,'glm_ros_AD.csv',sep=',',row.names = F, col.names=T)







############################ Abeta ############################
a=read.delim('~/5hmcData/Hydro_dataset_655_basic_10-11-2020.Neat1060.txt')
#colnames(a)
a=data.frame(a)
#head(a)
a=a[a[,'study']=='MAP',] # train: MAP; test: ROS
#a=a[a[,'study']=='ROS',]
### read 5hmC data 
hmC=readRDS('5hmC_amylsqrt.RDS')
hmC[1:2,1:3]

sum(hmC[,'projid']%in%a[,'projid']) 
outcome=merge.data.frame(a[,c('projid','study','amylsqrt','age_death','msex','pmi','educ')],hmC,by.x = "projid",by.y="projid")
dim(outcome)                         
table(outcome[,'msex'])
table(outcome[,'study'])
dim(outcome)

colnames(outcome)[ncol(outcome)]
n_peak=ncol(outcome)-7 # the first 7 columns in 'outcome' data are covariates
print('No. of peaks')
outcome[,'amylsqrt']=as.numeric(outcome[,'amylsqrt'])
outcome[,'age_death']=as.numeric(outcome[,'age_death'])
outcome[,'educ']=as.numeric(outcome[,'educ'])
outcome[,'pmi']=as.numeric(outcome[,'pmi'])
output=matrix(0,n_peak,4)
for (i in 8:(7+n_peak)){
  outcome[,'peak']=as.numeric(outcome[,i])
  fit=lm(amylsqrt~peak+age_death+msex+educ+pmi,outcome)
  ctable <- coef(summary(fit))
  output[j,]=ctable[rownames(ctable)=='peak',]
  j=j+1
}
dim(output)
output=cbind(colnames(outcome)[start:end],output)
colnames(output)=c('peak',"beta","se","t","p")
write.table(output,'lm_map_Abeta.csv',sep=',',row.names = F, col.names=T)
#write.table(output,'lm_ros_Abeta.csv',sep=',',row.names = F, col.names=T)







############################ Tau ############################
a=read.delim('~/5hmcData/Hydro_dataset_655_basic_10-11-2020.Neat1060.txt')
#colnames(a)
a=data.frame(a)
#head(a)
a=a[a[,'study']=='MAP',] # train: MAP; test: ROS
#a=a[a[,'study']=='ROS',]
### read 5hmC data 
hmC=readRDS('5hmC_tau.RDS') 
hmC[1:2,1:3]

sum(hmC[,'projid']%in%a[,'projid']) 
outcome=merge.data.frame(a[,c('projid','study','tangles','age_death','msex','pmi','educ')],hmC,by.x = "projid",by.y="projid")
dim(outcome)                         
table(outcome[,'msex'])
table(outcome[,'study'])
dim(outcome)

colnames(outcome)[ncol(outcome)]
n_peak=ncol(outcome)-7 # the first 7 columns in 'outcome' data are covariates
print('No. of peaks')
outcome[,'tangles']=as.numeric(outcome[,'tangles'])
outcome[,'age_death']=as.numeric(outcome[,'age_death'])
outcome[,'educ']=as.numeric(outcome[,'educ'])
outcome[,'pmi']=as.numeric(outcome[,'pmi'])
output=matrix(0,n_peak,4)
for (i in 8:(7+n_peak)){
  outcome[,'peak']=as.numeric(outcome[,i])
  fit=lm(tangles~peak+age_death+msex+educ+pmi,outcome)
  ctable <- coef(summary(fit))
  output[j,]=ctable[rownames(ctable)=='peak',]
  j=j+1
}
dim(output)
output=cbind(colnames(outcome)[start:end],output)
colnames(output)=c('peak',"beta","se","t","p")
write.table(output,'lm_map_Tau.csv',sep=',',row.names = F, col.names=T)
#write.table(output,'lm_ros_Tau.csv',sep=',',row.names = F, col.names=T)

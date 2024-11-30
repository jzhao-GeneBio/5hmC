####---------------------------------------------------------------####
####                  Demographics of Multi-omics Data             ####
####---------------------------------------------------------------####

cli.phe <- read.csv("H:/Documents/dataset_1009_basic_06-01-2023.csv",header = TRUE)
phe.go <- c("projid", "study","age_death","msex","ad_reagan", "amyloid","tangles")
cli.phe.sub <- cli.phe[,colnames(cli.phe)%in%phe.go]
mo.id <- read.csv("H:/Documents/ID_multi-omics.csv",header = TRUE)
hmc.phe <- cli.phe.sub[cli.phe.sub$projid%in%mo.id$id.5hmC,]
hmc.bar <- data.frame(Status=c("Case","Control","Case","Control","Case","Control"),
                      Cohort=c("MAP","MAP","ROS","ROS","All","All"),
                      Sample_size=c(392,208,288,162,680,370))

hmc.plot <- ggplot(data=hmc.bar, aes(x=Cohort, y=Sample_size, fill=Status)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=Sample_size, label=Sample_size), vjust=1.6, 
            color="darkblue", size=8)+
  scale_fill_brewer(palette=7)+labs(y="Sample size")+theme_classic()+
  theme(axis.text.y = element_text(face="bold",size=17),
        axis.text.x = element_text(face="bold",size=17),axis.title.y = element_text(face="bold",size=17),
        axis.title.x = element_text(face="bold",size=17),legend.text=element_text(face="bold",size=17),
        legend.title=element_text(face="bold",size=17))

###-----------------------###
gene.phe <- cli.phe.sub[cli.phe.sub$projid%in%mo.id$id.gene,]
gene.bar <- data.frame(Status=c("Case","Control","Case","Control","Case","Control"),
                       Cohort=c("MAP","MAP","ROS","ROS","All","All"),
                       Sample_size=c(316,175,200,128,516,303))

gene.plot <- ggplot(data=gene.bar, aes(x=Cohort, y=Sample_size, fill=Status)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=Sample_size, label=Sample_size), vjust=1.6, 
            color="darkblue", size=8)+
  scale_fill_brewer(palette=7)+labs(y="Sample size")+theme_classic()+
  theme(axis.text.y = element_text(face="bold",size=17),
        axis.text.x = element_text(face="bold",size=17),axis.title.y = element_text(face="bold",size=17),
        axis.title.x = element_text(face="bold",size=17),legend.text=element_text(face="bold",size=17),
        legend.title=element_text(face="bold",size=17))

###------------###
histone.phe <- cli.phe.sub[cli.phe.sub$projid%in%mo.id$id.histone,]
histone.bar <- data.frame(Status=c("Case","Control","Case","Control","Case","Control"),
                          Cohort=c("MAP","MAP","ROS","ROS","All","All"),
                          Sample_size=c(168,112,144,89,312,201))

histone.plot <- ggplot(data=histone.bar, aes(x=Cohort, y=Sample_size, fill=Status)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=Sample_size, label=Sample_size), vjust=1.6, 
            color="darkblue", size=8)+
  scale_fill_brewer(palette=7)+labs(y="Sample size")+theme_classic()+
  theme(axis.text.y = element_text(face="bold",size=17),
        axis.text.x = element_text(face="bold",size=17),axis.title.y = element_text(face="bold",size=17),
        axis.title.x = element_text(face="bold",size=17),legend.text=element_text(face="bold",size=17),
        legend.title=element_text(face="bold",size=17))

###--------------###
protein.phe <- cli.phe.sub[cli.phe.sub$projid%in%mo.id$id.protein,]
protein.bar <- data.frame(Status=c("Case","Control","Case","Control","Case","Control"),
                          Cohort=c("MAP","MAP","ROS","ROS","All","All"),
                          Sample_size=c(300,173,157,101,457,274))

protein.plot <- ggplot(data=protein.bar, aes(x=Cohort, y=Sample_size, fill=Status)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=Sample_size, label=Sample_size), vjust=1.6, 
            color="darkblue", size=8)+
  scale_fill_brewer(palette=7)+labs(y="Sample size")+theme_classic()+
  theme(axis.text.y = element_text(face="bold",size=17),
        axis.text.x = element_text(face="bold",size=17),axis.title.y = element_text(face="bold",size=17),
        axis.title.x = element_text(face="bold",size=17),legend.text=element_text(face="bold",size=17),
        legend.title=element_text(face="bold",size=17))

###-------------------###
hmc.phe2 <- hmc.phe[,2:4]
hmc.phe2[,1] <- "All"
hmc.phe.stack <- rbind(hmc.phe[,2:4],hmc.phe2)
colnames(hmc.phe.stack) <- c("Cohort","Age","Gender")
hmc.phe.stack[hmc.phe.stack[,3]==1,3] <- "Male"
hmc.phe.stack[hmc.phe.stack[,3]==0,3] <- "Female"

hmc.age.plot <- ggplot(hmc.phe.stack, aes(x=Cohort, y=Age, fill=Gender)) + 
  geom_boxplot()+
  theme(axis.text.y = element_text(face="bold",size=17),
        axis.text.x = element_text(face="bold",size=17),axis.title.y = element_text(face="bold",size=17),
        axis.title.x = element_text(face="bold",size=17),legend.text=element_text(face="bold",size=17),
        legend.title=element_text(face="bold",size=17))

###--------------------------------###
gene.phe2 <- gene.phe[,2:4]
gene.phe2[,1] <- "All"
gene.phe.stack <- rbind(gene.phe[,2:4],gene.phe2)
colnames(gene.phe.stack) <- c("Cohort","Age","Gender")
gene.phe.stack[gene.phe.stack[,3]==1,3] <- "Male"
gene.phe.stack[gene.phe.stack[,3]==0,3] <- "Female"

gene.age.plot <- ggplot(gene.phe.stack, aes(x=Cohort, y=Age, fill=Gender)) + 
  geom_boxplot()+
  theme(axis.text.y = element_text(face="bold",size=17),
        axis.text.x = element_text(face="bold",size=17),axis.title.y = element_text(face="bold",size=17),
        axis.title.x = element_text(face="bold",size=17),legend.text=element_text(face="bold",size=17),
        legend.title=element_text(face="bold",size=17))

###----------------------------###
histone.phe2 <- histone.phe[,2:4]
histone.phe2[,1] <- "All"
histone.phe.stack <- rbind(histone.phe[,2:4],histone.phe2)
colnames(histone.phe.stack) <- c("Cohort","Age","Gender")
histone.phe.stack[histone.phe.stack[,3]==1,3] <- "Male"
histone.phe.stack[histone.phe.stack[,3]==0,3] <- "Female"

histone.age.plot <- ggplot(histone.phe.stack, aes(x=Cohort, y=Age, fill=Gender)) + 
  geom_boxplot()+
  theme(axis.text.y = element_text(face="bold",size=17),
        axis.text.x = element_text(face="bold",size=17),axis.title.y = element_text(face="bold",size=17),
        axis.title.x = element_text(face="bold",size=17),legend.text=element_text(face="bold",size=17),
        legend.title=element_text(face="bold",size=17))

###-----------------------###
protein.phe2 <- protein.phe[,2:4]
protein.phe2[,1] <- "All"
protein.phe.stack <- rbind(protein.phe[,2:4],protein.phe2)
colnames(protein.phe.stack) <- c("Cohort","Age","Gender")
protein.phe.stack[protein.phe.stack[,3]==1,3] <- "Male"
protein.phe.stack[protein.phe.stack[,3]==0,3] <- "Female"

protein.age.plot <- ggplot(protein.phe.stack, aes(x=Cohort, y=Age, fill=Gender)) + 
  geom_boxplot()+
  theme(axis.text.y = element_text(face="bold",size=17),
        axis.text.x = element_text(face="bold",size=17),axis.title.y = element_text(face="bold",size=17),
        axis.title.x = element_text(face="bold",size=17),legend.text=element_text(face="bold",size=17),
        legend.title=element_text(face="bold",size=17))

####------------------------####
jpeg(file = "H:/Documents/hmc_omics_plot.jpeg",width = 20000, height = 12000,res=1050)
ggarrange(ggarrange(hmc.plot,gene.plot,histone.plot,protein.plot,
                    common.legend = TRUE, legend = "bottom",ncol=4, nrow=1),
          ggarrange(hmc.age.plot,gene.age.plot,histone.age.plot,protein.age.plot,
                    common.legend = TRUE, legend = "bottom",ncol=4, nrow=1), ncol=1, nrow=2)
dev.off() 


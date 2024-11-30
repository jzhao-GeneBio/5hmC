####---------------------------------------------------------------####
####                     colocalization analysis                   ####
####---------------------------------------------------------------####

folder <- c("As","En","Ex","In","Mi","Ol","Op","Pe")
library(coloc)
best.causal.snp <- c()
sink("/path/Celine_Extend500kb_GRCh38_Coloc_Res_4_18_2024.txt")
all.peak.gene <- read.table("/path/AD_Reagan_Amy_Tan_peak_anno.txt",header=TRUE)
for(i in 1:2821){
  peak <- all.peak.gene[i,4]
  gene <- all.peak.gene[i,5]
  if(is.na(gene)){
    next
  }else{
    celine <- read.table(paste0("/path/All_Peak_Celine_SummaryData_SNP155_GRCh38_Extend500kb_AD_Reagan_Amy_Tan_4_18_2024/",peak,"_Celine_SummaryData_SNP155_GRCh38_Extend500kb_AD_Reagan_Amy_Tan_4_18_2024.txt"),header = TRUE)
    celine.delete <- which(duplicated(celine[,3]) == TRUE)
    if(length(celine.delete) > 0){
      celine <- celine[-celine.delete,]}
    celine.list <- list()
    celine.list[[1]] <- celine$b
    celine.list[[2]] <- (celine$se)^2
    celine.list[[3]] <- celine$SNP
    celine.list[[4]] <- celine$pos
    celine.list[[5]] <- "quant"
    celine.list[[6]] <- celine$N
    celine.list[[7]] <- pmin(celine$freq, 1 - celine$freq)
    names(celine.list) <- c("beta","varbeta","snp","position","type","N","MAF")
    for(j in 1:8){
      if(!file.exists(paste0("/path/CT_eQTL_GRCh38_37464/",folder[j],"/",folder[j],"_",gene,"_eQTL_SumStats_GRCh38.txt"))){
        next
      }else if(file.size(paste0("/path/CT_eQTL_GRCh38_37464/",folder[j],"/",folder[j],"_",gene,"_eQTL_SumStats_GRCh38.txt")) == 58){
        next
      } else{
        file <- read.table(paste0("/path/CT_eQTL_GRCh38_37464/",folder[j],"/",folder[j],"_",gene,"_eQTL_SumStats_GRCh38.txt"),header = TRUE)
        file.delete1 <- which(duplicated(file[,3]) == TRUE)
        file.delete2 <- which((file[,6]<= 0)|(file[,6] >= 1))
        file.delete <- sort(unique(c(file.delete1,file.delete2)))
        if(length(file.delete) > 0){
          file <- file[-file.delete,]}
        file.list <- list()
        file.list[[1]] <- file$Beta
        file.list[[2]] <- (file$Se)^2
        file.list[[3]] <- file$SNP
        file.list[[4]] <- file$Pos
        file.list[[5]] <- "quant"
        file.list[[6]] <- 192
        file.list[[7]] <- file$MAF
        names(file.list) <- c("beta","varbeta","snp","position","type","N","MAF")
        possibleError <- tryCatch({res <- coloc.abf(dataset1=celine.list,dataset2=file.list)}, error=function(e) e)
        if(inherits(possibleError, "error")){ 
          next
        }else if(res$summary[6] < 0.8){
          next 
        }else{
          best <- subset(res$results,SNP.PP.H4 == max(res$results$SNP.PP.H4))[,c(1,2,12)]
          best.causal.snp <- rbind(best.causal.snp, c(peak, gene,folder[j], best, res$summary[1],res$summary[2],res$summary[3],res$summary[4],res$summary[5],res$summary[6]))
        }
      }
    }
  }
}
colnames(best.causal.snp) <- c("Peak","Gene","CT","SNP","Pos","SNP.PP.H4", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
write.table(best.causal.snp, file = "/path/Celine_Coloc_Res/bf_best_causal_snp_Celine_Coloc.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)




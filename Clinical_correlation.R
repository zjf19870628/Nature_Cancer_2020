ALL_mutations <- read.csv("../../ALL_Figures/data_coding/ALL_coding_mutation_table.csv",
                          colClasses = "character")
ALL_mutations <- ALL_mutations[ALL_mutations$Synonymous..0..vs..non.Synonymous..1.=="1" &
                                 ALL_mutations$Splicing.site.mutation..0..no..1..yes.=="0",]
ALL_mutations_D <- ALL_mutations[ALL_mutations$Diagnosis.variant.frequency!="0",]
ALL_mutations_R <- ALL_mutations[ALL_mutations$Relapse.variant.frequency!="0",]

CNV_Del_summary <- read.table("CNV_Del_summary.txt",header=TRUE,
                              sep="\t",fill=TRUE,colClasses = "character")

mut_Genes <- read.csv("../Recurrence_analysis/Recurrence_analysis_I3.csv",
                      colClasses = "character")[1:21,1]
#cnv_Genes <- c("CDKN2A", "BRF1", "IKZF1", "ETV6","G6PD", "PHF6","NAA10", "RPL10", "HPRT1", "MED12", "USP9X")
cnv_Genes <- c("CDKN2A", "IKZF1", "ETV6", "HPRT1")

N <- length(c(mut_Genes,cnv_Genes))

#gene fusion
library(stringr)
TARGET_RNA_samples <- intersect(sapply(read.csv("~/Dropbox/CUMC/Projects/ALL/TARGET/RNAseq/pair_samples_with_SRA_ID.csv",colClasses = "character")[,1],
                                       function(x) substring(x,11)),unique(ALL_mutations$Sample))
RNA_samples <- c("4","5","35","37","6787","11936","6206","6329","14840","8711","7300","8148","8173","5818","15342",
                 TARGET_RNA_samples)

fusion_table <- read.csv("~/Dropbox/CUMC/Projects/ALL/RNAseq/Primary/Primary_report_star-fusion.csv",colClasses = "character")
fusion_table$Case.ID <- str_extract(fusion_table$Sample.ID, "[0-9]+")
fusion_table$Sample <- str_extract(fusion_table$Sample.ID, "[A-Z]+")

TARGET_fusion_table <- read.csv("~/Dropbox/CUMC/Projects/ALL/TARGET/RNAseq/TARGET_report_star-fusion.csv",colClasses = "character")
TARGET_fusion_table[,"Case.ID"] <- sapply(TARGET_fusion_table[,"Case.ID"],function(x) unlist(strsplit(x,split="-"))[3])
TARGET_fusion_table[,"Sample"] <- sapply(TARGET_fusion_table[,"Sample"],function(x) substring(x,1,1))

fusion_list <- c("ETV6--RUNX1","PICALM--MLLT10","NUP214--ABL1")

fus_samples_D <- list()
for(fusion in fusion_list){
  fusion_rev <- paste0(rev(unlist(strsplit(fusion,split="--"))),collapse = "--")
  fus_samples_D[[fusion]] <- c(unique(fusion_table$Case.ID[fusion_table$X.FusionName %in% c(fusion,fusion_rev) & 
                                                           fusion_table$Sample=="D"]),
                               unique(TARGET_fusion_table$Case.ID[TARGET_fusion_table$X.FusionName %in% c(fusion,fusion_rev) & 
                                                                    TARGET_fusion_table$Sample=="D"])
  )
}

fus_samples_R <- list()
for(fusion in fusion_list){
  fusion_rev <- paste0(rev(unlist(strsplit(fusion,split="--"))),collapse = "--")
  fus_samples_R[[fusion]] <- c(unique(fusion_table$Case.ID[fusion_table$X.FusionName %in% c(fusion,fusion_rev) & 
                                                             fusion_table$Sample=="R"]),
                               unique(TARGET_fusion_table$Case.ID[TARGET_fusion_table$X.FusionName %in% c(fusion,fusion_rev) & 
                                                                    TARGET_fusion_table$Sample=="R"])
  )
}

####################################################
Sample_table <- read.csv("../ALL_Sample_Full_List.csv",header=TRUE,colClasses="character")
Sample_table <- Sample_table[!(grepl("Included",Sample_table$Annotation)),]
T_samples <- Sample_table$Sample.ID[Sample_table$Origin=="T"]
B_samples <- Sample_table$Sample.ID[Sample_table$Origin=="B"]

#Diagnosis

mut_samples_D <- list()
for(gene in mut_Genes){
  mut_samples_D[[paste0(gene,"_MUT")]] <- unique(ALL_mutations_D$Sample[ALL_mutations_D$Gene==gene])
}
for(gene in cnv_Genes){
  mut_samples_D[[paste0(gene,"_DEL")]] <- unlist(strsplit(CNV_Del_summary$Diagnosis[CNV_Del_summary$Gene==gene],split=","))
}

#Relapse

mut_samples_R <- list()
for(gene in mut_Genes){
  mut_samples_R[[paste0(gene,"_MUT")]] <- unique(ALL_mutations_R$Sample[ALL_mutations_R$Gene==gene])
}
for(gene in cnv_Genes){
  mut_samples_R[[paste0(gene,"_DEL")]] <- unlist(strsplit(CNV_Del_summary$Relapse[CNV_Del_summary$Gene==gene],split=","))
}


library(RColorBrewer)

mycol <- brewer.pal(9,"Set1")[c(2,1)]

pdf("Clinicla_correlation.pdf",width=6,height=8,useDingbats = FALSE)
par(mar=c(2,6,7,4))

plot(0,
     xlim=c(1,12),
     ylim=c(-N,-1),
     type="n",
     xaxt="n",
     yaxt="n",
     bty="n",
     xlab=NA,
     ylab=NA)

text(0.5,-c(1:N),
     labels=names(mut_samples_D),
     xpd=TRUE,
     cex=0.6,
     pos=2,
     font=3)
text(c(1:12),1,
     labels=c("Diagnosis","Relapse","T-ALL","B-ALL",fusion_list,"T-ALL","B-ALL",fusion_list),
     xpd=TRUE,
     cex=0.8,
     adj=0,
     font=3,
     srt=45)

rect(2.6,0,7.4,0.5,col=mycol[1],border=NA,xpd=TRUE)
rect(7.6,0,12.4,0.5,col=mycol[2],border=NA,xpd=TRUE)

for(i in 1:12){
  abline(v=i,lty=1,col="grey90")
}

for(i in 1:N){
  abline(h=-i,lty=1,col="grey90")
}

#D or R
for(i in 1:N){
  t <- fisher.test(matrix(c(length(mut_samples_D[[i]]),175-length(mut_samples_D[[i]]),length(mut_samples_R[[i]]),175-length(mut_samples_R[[i]])),byrow = TRUE,ncol=2))
  if(t$p.value < 0.1){
    if(t$estimate > 1){
      points(1,-i,pch=16,cex=-log10(max(t$p.value,0.001)),col=brewer.pal(9,"Set1")[5],xpd=TRUE)
    }
    else{
      points(2,-i,pch=16,cex=-log10(max(t$p.value,0.001)),col=brewer.pal(9,"Set1")[5],xpd=TRUE)
    }
    
  }
}

#T or B in D
for(i in 1:N){
  T_mut <- length(intersect(T_samples,mut_samples_D[[i]]))
  B_mut <- length(intersect(B_samples,mut_samples_D[[i]]))
  t <- fisher.test(matrix(c(T_mut,length(T_samples)-T_mut,B_mut,length(B_samples)-B_mut),byrow = TRUE,ncol=2))
  if(t$p.value < 0.1){
    if(t$estimate > 1){
      points(3,-i,pch=16,cex=-log10(max(t$p.value,0.001)),col=brewer.pal(9,"Set1")[5],xpd=TRUE)
    }
    else{
      points(4,-i,pch=16,cex=-log10(max(t$p.value,0.001)),col=brewer.pal(9,"Set1")[5],xpd=TRUE)
    }
  }
}

#FUS in D
for(fus in fusion_list){
  for(i in 1:N){
    t <- fisher.test(matrix(c(length(intersect(fus_samples_D[[fus]],mut_samples_D[[i]])),
                              length(intersect(setdiff(mut_samples_D[[i]],fus_samples_D[[fus]]),RNA_samples)),
                              length(setdiff(fus_samples_D[[fus]],mut_samples_D[[i]])),
                              length(setdiff(setdiff(RNA_samples,fus_samples_D[[fus]]),mut_samples_D[[i]]))),
                            byrow = TRUE,ncol=2))
    if(t$p.value < 0.1){
      points(4 + which(fusion_list == fus),-i,pch=16,cex=-log10(max(t$p.value,0.001)),col=brewer.pal(9,"Set1")[5],xpd=TRUE)
    }
  }
}

##############################################################

#T or B in R
for(i in 1:N){
  T_mut <- length(intersect(T_samples,mut_samples_R[[i]]))
  B_mut <- length(intersect(B_samples,mut_samples_R[[i]]))
  t <- fisher.test(matrix(c(T_mut,length(T_samples)-T_mut,B_mut,length(B_samples)-B_mut),byrow = TRUE,ncol=2))
  if(t$p.value < 0.1){
    if(t$estimate > 1){
      points(8,-i,pch=16,cex=-log10(max(t$p.value,0.001)),col=brewer.pal(9,"Set1")[5],xpd=TRUE)
    }
    else{
      points(9,-i,pch=16,cex=-log10(max(t$p.value,0.001)),col=brewer.pal(9,"Set1")[5],xpd=TRUE)
    }
  }
}

#FUS in R
for(fus in fusion_list){
  for(i in 1:N){
    t <- fisher.test(matrix(c(length(intersect(fus_samples_R[[fus]],mut_samples_R[[i]])),
                              length(intersect(setdiff(mut_samples_R[[i]],fus_samples_R[[fus]]),RNA_samples)),
                              length(setdiff(fus_samples_R[[fus]],mut_samples_R[[i]])),
                              length(setdiff(setdiff(RNA_samples,fus_samples_R[[fus]]),mut_samples_R[[i]]))),
                            byrow = TRUE,ncol=2))
    if(t$p.value < 0.1){
      points(9 + which(fusion_list == fus),-i,pch=16,cex=-log10(max(t$p.value,0.001)),col=brewer.pal(9,"Set1")[5],xpd=TRUE)
    }
  }
}


dev.off()

pdf("Clinicla_correlation_legend.pdf",width=6,height=8,useDingbats = FALSE)
par(mar=c(2,6,7,4))

plot(0,
     xlim=c(0,1),
     ylim=c(0,3),
     type="n",
     xaxt="n",
     yaxt="n",
     bty="n",
     xlab=NA,
     ylab=NA)

points(0,3,pch=16,cex=1,col=brewer.pal(9,"Set1")[5],xpd=TRUE)
points(0,2,pch=16,cex=2,col=brewer.pal(9,"Set1")[5],xpd=TRUE)
points(0,1,pch=16,cex=3,col=brewer.pal(9,"Set1")[5],xpd=TRUE)

dev.off()




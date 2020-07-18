library(ComplexHeatmap)
library(RColorBrewer)

mycols <- brewer.pal(7,"Set1")[c(3,2,1)]

###################################################################
#mutation landscape
Recurrent_genes <- read.csv("/Users/junfeizhao/Dropbox/CUMC/Projects/ALL/Figures/Recurrence_analysis/Recurrence_analysis_I3.csv",
                            header=TRUE,colClasses="character")
ALL_genes<- Recurrent_genes$Gene[as.numeric(Recurrent_genes$Number.of.Mutant.Cases)>=4]

df <- read.csv ("/Users/junfeizhao/Dropbox/CUMC/Projects/ALL/ALL_Figures/data_coding/ALL_coding_mutation_table.csv", 
                header = TRUE,colClasses="character")
Samples <- unique(df[,"Sample"])
Sample_table <- read.csv("../ALL_Sample_Full_List.csv",header=TRUE,colClasses="character")
Sample_table <- Sample_table[!grepl("Included",Sample_table$Annotation),]
rownames(Sample_table) <- Sample_table$Sample.ID
Origin <- Sample_table[Samples,"Origin",drop=FALSE]
Platform <- rep("Illumina",nrow(Origin))
names(Platform) <- rownames(Origin)
Platform[startsWith(names(Platform),"P")] <- "CG"
Adult_IDs <- c(scan("../../Figures/Evolution_tree/Adult_IDs.txt",what = "character"),
               scan("../../TARGET/WGS/Adult_IDs.txt",what = "character"))
Ped_or_Adult <- rep("Ped",nrow(Origin))
names(Ped_or_Adult) <- rownames(Origin)
Ped_or_Adult[names(Ped_or_Adult) %in% Adult_IDs] <- "Adult"
Origin <- cbind(Origin,Platform,Ped_or_Adult)
Origin_col <- brewer.pal(10,"Paired")[c(7,8)]
names(Origin_col) <- c("T","B")
Platform_col <- brewer.pal(10,"Paired")[c(9,10)]
names(Platform_col) <- c("Illumina","CG")
Ped_or_Adult_col <- brewer.pal(10,"Paired")[c(9,10)]
names(Ped_or_Adult_col) <- c("Ped","Adult")

mutation_burden <- matrix(0,nrow=3,ncol=length(Samples))
dimnames(mutation_burden) <- list(c("D_R","D","R"),Samples)
for(S in Samples){
  tmp <- df[df$Sample==S,]
  mutation_burden["D_R",S] <- sum(as.numeric(tmp$Diagnosis.variant.frequency)>0 & as.numeric(tmp$Relapse.variant.frequency)>0)
  mutation_burden["D",S] <- sum(as.numeric(tmp$Diagnosis.variant.frequency)>0 & as.numeric(tmp$Relapse.variant.frequency)==0)
  mutation_burden["R",S] <- sum(as.numeric(tmp$Diagnosis.variant.frequency)==0 & as.numeric(tmp$Relapse.variant.frequency)>0)
}

HM_samples <- colnames(mutation_burden)[colSums(mutation_burden)>200]
col_fontsize <- rep(0,ncol(mutation_burden))
names(col_fontsize) <- colnames(mutation_burden)
col_fontsize[HM_samples] <- 4

Hyper_Mut <- rep("No",nrow(Origin))
names(Hyper_Mut) <- rownames(Origin)
Hyper_Mut[HM_samples] <- "Yes"
Origin <- cbind(Origin,Hyper_Mut)
Hyper_Mut_col <- brewer.pal(10,"Paired")[c(1,2)]
names(Hyper_Mut_col) <- c("No","Yes")

driver_mutations <- df[df$Gene %in% ALL_genes,]
ALL_genes <- intersect(ALL_genes,driver_mutations$Gene)
mat <- matrix("",nrow=length(ALL_genes),ncol=length(Samples),
              dimnames=list(ALL_genes,Samples))
for(i in 1:nrow(driver_mutations)){
  ID <- driver_mutations[i,"Sample"]
  Gene <- driver_mutations[i,"Gene"]
  if(driver_mutations[i,"Diagnosis.variant.frequency"] != "0" & driver_mutations[i,"Relapse.variant.frequency"] != "0"){
    if(any(driver_mutations[i,"Reference.Sequence"]=="-", driver_mutations[i,"Variant.Sequence"]=="-") | any(nchar(driver_mutations[i,"Reference.Sequence"])>1, nchar(driver_mutations[i,"Variant.Sequence"])>1)){
      mat[Gene,ID] <- "D_R;Truncating"
    }
    else{
      if(grepl("\\*",driver_mutations[i,"Predicted.protein.product" ])){
        mat[Gene,ID] <- "D_R;Truncating"
      }
      else{
        if(driver_mutations[i,"Splicing.site.mutation..0..no..1..yes."]!="0"){
          mat[Gene,ID] <- "D_R"
        }
        else{
          mat[Gene,ID] <- "D_R"
        }
      }
    }
  }
  if(driver_mutations[i,"Diagnosis.variant.frequency"] != "0" & driver_mutations[i,"Relapse.variant.frequency"] == "0"){
    if(any(driver_mutations[i,"Reference.Sequence"]=="-", driver_mutations[i,"Variant.Sequence"]=="-") | any(nchar(driver_mutations[i,"Reference.Sequence"])>1, nchar(driver_mutations[i,"Variant.Sequence"])>1)){
      mat[Gene,ID] <- "D;Truncating"
    }
    else{
      if(grepl("\\*",driver_mutations[i,"Predicted.protein.product" ])){
        mat[Gene,ID] <- "D;Truncating"
      }
      else{
        if(driver_mutations[i,"Splicing.site.mutation..0..no..1..yes."]!="0"){
          mat[Gene,ID] <- "D"
        }
        else{
          mat[Gene,ID] <- "D"
        }
      }
    }
  } 
  if(driver_mutations[i,"Diagnosis.variant.frequency"] == "0" & driver_mutations[i,"Relapse.variant.frequency"] != "0"){
    if(any(driver_mutations[i,"Reference.Sequence"]=="-", driver_mutations[i,"Variant.Sequence"]=="-") | any(nchar(driver_mutations[i,"Reference.Sequence"])>1, nchar(driver_mutations[i,"Variant.Sequence"])>1)){
      mat[Gene,ID] <- "R;Truncating"
    }
    else{
      if(grepl("\\*",driver_mutations[i,"Predicted.protein.product" ])){
        mat[Gene,ID] <- "R;Truncating"
      }
      else{
        if(driver_mutations[i,"Splicing.site.mutation..0..no..1..yes."]!="0"){
          mat[Gene,ID] <- "R"
        }
        else{
          mat[Gene,ID] <- "R"
        }
      }
    }
  } 
}

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "grey98", col = NA))
  },
  D = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycols[2], col = NA))
  },
  R = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycols[3], col = NA))
  },
  D_R = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = mycols[1], col = NA))
  },
  Truncating = function(x, y, w, h) {
    grid.points(x,y,pch=4,size = unit(0.3, "char"))
  }
)
col = c("D" = mycols[2], "R" = mycols[3], "D_R" = mycols[1],"Truncating" = "black")

###################################################################
#copy number landscape
CNV_genes <- c("CDKN2A", "IKZF1", "ETV6", "HPRT1","RB1","TP53")
CNV_data <- data.matrix(read.table("../../WES_Combined_analysis/CNV/CNV_by_genes.txt",
                       header=TRUE,row.names=1,colClasses="character",check.names = FALSE))
CNV_mat <- matrix("",nrow=length(CNV_genes),ncol=length(Samples),
              dimnames=list(CNV_genes,Samples))
for(Gene in CNV_genes){
  for(ID in Samples){
    if(paste0(ID,"-D") %in% colnames(CNV_data) & paste0(ID,"-R") %in% colnames(CNV_data)){
      if(CNV_data[Gene,paste0(ID,"-D")] < 2 & CNV_data[Gene,paste0(ID,"-R")] < 2){
        CNV_mat[Gene,ID] <- "D_R"
      }
      if(CNV_data[Gene,paste0(ID,"-D")] < 2 & CNV_data[Gene,paste0(ID,"-R")] >= 2){
        CNV_mat[Gene,ID] <- "D"
      } 
      if(CNV_data[Gene,paste0(ID,"-D")] >= 2 & CNV_data[Gene,paste0(ID,"-R")] < 2){
        CNV_mat[Gene,ID] <- "R"
      } 
    }
  }
}
rownames(CNV_mat) <- c("CDKN2A/B_DEL", "IKZF1_DEL", "ETV6_DEL", "HPRT1_DEL","RB1_DEL","TP53_DEL")

#############################################################
#gene fusion
fusion_table <- read.csv("~/Dropbox/CUMC/Projects/ALL/RNAseq/Primary/Primary_report_star-fusion.csv",colClasses = "character")
TARGET_fusion_table <- read.csv("~/Dropbox/CUMC/Projects/ALL/TARGET/RNAseq/TARGET_report_star-fusion.csv",colClasses = "character")
TARGET_fusion_table[,"Case.ID"] <- sapply(TARGET_fusion_table[,"Case.ID"],function(x) unlist(strsplit(x,split="-"))[3])
TARGET_fusion_table[,"Sample"] <- sapply(TARGET_fusion_table[,"Sample"],function(x) substring(x,1,1))

fusion_list <- c("PICALM--MLLT10","NUP214--ABL1")
Fusion_mat <- matrix("",nrow=length(fusion_list),ncol=length(Samples),
                  dimnames=list(fusion_list,Samples))
for(fusion in fusion_list){
  fusion_ <- paste0(rev(unlist(strsplit(fusion,"--"))),collapse = "--")
  Sample_IDs <- unique(fusion_table$Sample.ID[fusion_table$X.FusionName %in% c(fusion,fusion_)])
  Patient_IDs <- unique(sapply(Sample_IDs,function(x) unlist(strsplit(x,split="-"))[2]))
  for(ID in Patient_IDs){
    Fusion_mat[fusion,ID] <- paste0(sapply(Sample_IDs[endsWith(Sample_IDs,ID)],function(x) unlist(strsplit(x,split="-"))[1]),collapse = "_")
  }
}
for(fusion in fusion_list){
  fusion_ <- paste0(rev(unlist(strsplit(fusion,"--"))),collapse = "--")
  tmp <- TARGET_fusion_table[TARGET_fusion_table$X.FusionName %in% c(fusion,fusion_),,drop=FALSE]
  if(nrow(tmp) > 0){
    Case_IDs <- unique(tmp$Case.ID)
    for(ID in Case_IDs){
      Fusion_mat[fusion,ID] <- paste0(unique(tmp$Sample[tmp$Case.ID==ID]),collapse = "_")
    }
  }
}
  
n_mutated_patients <- matrix(0,nrow=length(c(ALL_genes,rownames(CNV_mat),fusion_list)),ncol=3)
dimnames(n_mutated_patients) <- list(c(ALL_genes,rownames(CNV_mat),fusion_list),c("D_R","D","R"))
for(gene in ALL_genes){
  tmp <- sapply(mat[gene,],function(x) unlist(strsplit(x,split=";"))[1])
  n_mutated_patients[gene,"D_R"] <- sum(tmp=="D_R",na.rm = TRUE)
  n_mutated_patients[gene,"D"] <- sum(tmp=="D",na.rm = TRUE)
  n_mutated_patients[gene,"R"] <- sum(tmp=="R",na.rm = TRUE)
}
for(gene in rownames(CNV_mat)){
  tmp <- sapply(CNV_mat[gene,],function(x) unlist(strsplit(x,split=";"))[1])
  n_mutated_patients[gene,"D_R"] <- sum(tmp=="D_R",na.rm = TRUE)
  n_mutated_patients[gene,"D"] <- sum(tmp=="D",na.rm = TRUE)
  n_mutated_patients[gene,"R"] <- sum(tmp=="R",na.rm = TRUE)
}
for(fusion in rownames(Fusion_mat)){
  tmp <- Fusion_mat[fusion,]
  n_mutated_patients[fusion,"D_R"] <- sum(tmp=="D_R",na.rm = TRUE)
  n_mutated_patients[fusion,"D"] <- sum(tmp=="D",na.rm = TRUE)
  n_mutated_patients[fusion,"R"] <- sum(tmp=="R",na.rm = TRUE)
}

ht = oncoPrint(rbind(mat,CNV_mat,Fusion_mat), get_type = function(x) strsplit(x, ";")[[1]],
                        alter_fun = alter_fun, col = col, 
                        row_names_gp = gpar(fontsize = 8),
                        column_names_gp = gpar(fontsize = col_fontsize),
                        show_pct = FALSE,
                        show_column_names = FALSE,
                        top_annotation = HeatmapAnnotation(column_barplot = anno_barplot(t(mutation_burden), axis = TRUE, border=FALSE, ylim = c(0, 200), gp = gpar(fill = mycols, col = mycols),height=unit(4, "cm"))),
                        bottom_annotation = HeatmapAnnotation(df = Origin[,c(1,3,4),drop=FALSE],col=list(Origin=Origin_col,Ped_or_Adult=Ped_or_Adult_col,Hyper_Mut=Hyper_Mut_col)),
                        right_annotation = NULL,
                        row_order = NULL,
                        row_split = c(rep("A",length(ALL_genes)),rep("B",length(CNV_genes)),rep("C",length(fusion_list)))
               ) +
  rowAnnotation(row_barplot = anno_barplot(n_mutated_patients, ylim = c(0, 80), gp = gpar(fill = mycols, col = mycols),
                                               axis = TRUE, axis_param = list(side = "top"), border=FALSE, which = "row"), width = unit(2, "cm"))

#make the plot
pdf("mutational_landscape_test.pdf",width=16,height=8)
par(mar=c(2,5,2,5))
draw(ht)
dev.off()


# #######################
# pdf("legend.pdf",height=1.5,width=1.5)
# 
# par(mar=c(0,0,0,0))
# 
# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# legend("topleft", legend =c('Frameshift', 'Nonsense', 'Splicing'), 
#        pch=c(8,4,3), cex=1, bty='n',
#        col = c('black'))
# #mtext("Alteration", at=0.2, cex=2)
# 
# dev.off()



Gene_list <- read.csv("../Recurrence_analysis/Recurrence_analysis_I3.csv",
                      colClasses = "character")[1:21,1]
#Gene_list <- c("NT5C2","KRAS","NRAS","WT1","NOTCH1","PHF6","DNM2","JAK2","FBXW7","FLT3",
#               "CREBBP","KMT2D","USP9X","ZFHX3","CACNA1H","TP53","JAK3","SHROOM3","EPHA3","NR3C1")
mutation_table <- read.csv("../../ALL_WES_Figures/data_coding/ALL_coding_mutation_table.csv",
                             colClasses = "character")
mutation_table <- mutation_table[mutation_table$Synonymous..0..vs..non.Synonymous..1.=="1" &
                                     mutation_table$Splicing.site.mutation..0..no..1..yes.=="0",]

mut_num <- sapply(Gene_list,function(x) length(unique(mutation_table$Sample[mutation_table$Gene==x])))
mutation_table_D <- mutation_table[apply(mutation_table,1,function(x) as.numeric(x["Relapse.variant.frequency"])==0),]
D_num <- sapply(Gene_list,function(x) length(unique(mutation_table_D$Sample[mutation_table_D$Gene==x])))
mutation_table_R <- mutation_table[apply(mutation_table,1,function(x) as.numeric(x["Diagnosis.variant.frequency"])==0),]
R_num <- sapply(Gene_list,function(x) length(unique(mutation_table_R$Sample[mutation_table_R$Gene==x])))
mutation_table_D_R <- mutation_table[apply(mutation_table,1,function(x) as.numeric(x["Diagnosis.variant.frequency"])>0 & as.numeric(x["Relapse.variant.frequency"])>0),]
D_R_num <- sapply(Gene_list,function(x) length(unique(mutation_table_D_R$Sample[mutation_table_D_R$Gene==x])))

circos_table <- rbind(D_num,D_R_num,R_num)
dimnames(circos_table) <- list(c("Diagnosis","Diagnosis & Relapse","Relapse"),
                               Gene_list)
#write.table(circos_table,"circos_table.txt",
#            row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)

Gene_list_exp <- as.vector(sapply(Gene_list,function(x) c(paste0(x,"_D"),paste0(x,"_D_R"),paste0(x,"_R"))))
circos_table_exp <- matrix(0,nrow=3,ncol=length(Gene_list_exp))
for(i in 1:length(Gene_list)){
  circos_table_exp[1,1+(i-1)*3] <- circos_table[1,i]
  circos_table_exp[2,2+(i-1)*3] <- circos_table[2,i]
  circos_table_exp[3,3+(i-1)*3] <- circos_table[3,i]
}
dimnames(circos_table_exp) <- list(c("Diagnosis","Diagnosis & Relapse","Relapse"),
                               Gene_list_exp)
library(circlize)
library(RColorBrewer)

mycol <- brewer.pal(9,"Set1")[c(2,3,1)]

pdf("circos_plot.pdf",width=8,height=8)
par(mar=c(2,2,2,2))

grid.col <- rep(mycol,length(Gene_list)+1)
names(grid.col) <- c(c("Diagnosis","Diagnosis & Relapse","Relapse"),Gene_list_exp)
circos.par(gap.after = c(rep(3, nrow(circos_table)-1), 10, rep(c(0,0,3), length(Gene_list)-1)[c(-22,-46)],0,0,10))
chordDiagram(circos_table_exp[,c(-24,-48)],
             grid.col = grid.col,
             order = c(c("Diagnosis","Diagnosis & Relapse","Relapse"),as.vector(sapply(Gene_list,function(x) c(paste0(x,"_R"),paste0(x,"_D_R"),paste0(x,"_D"))))[c(-22,-46)]),
             annotationTrack = c("grid"),
             preAllocateTracks = list(track.height = uh(4, "mm"),
                                      track.margin = c(uh(4, "mm"),0)))

highlight.sector("Diagnosis", track.index = 1, col = "white", 
                 text = "Diagnosis", cex = 1.5, text.col = "black", niceFacing = TRUE)
highlight.sector("Diagnosis & Relapse", track.index = 1, col = "white", 
                 text = "Diagnosis & Relapse", cex = 1.5, text.col = "black", niceFacing = TRUE)
highlight.sector("Relapse", track.index = 1, col = "white", 
                 text = "Relapse", cex = 1.5, text.col = "black", niceFacing = TRUE)
for(gene in c("NRAS","KRAS","NT5C2","NOTCH1","CREBBP","TP53", "WT1","PHF6","DNM2","JAK3","PTPN11","FBXW7","MYC","FLT3")){
  sector_names <- c(paste0(gene,"_D"),paste0(gene,"_D_R"),paste0(gene,"_R"))
  index <- c(circos_table_exp[1,paste0(gene,"_D")],
             circos_table_exp[2,paste0(gene,"_D_R")],
             circos_table_exp[3,paste0(gene,"_R")]) > 0
  sector_names <- sector_names[index]
  highlight.sector(sector_names, track.index = 1, col = NA, 
                   text = gene, cex = 1, text.col = "black", niceFacing = TRUE)
}
for(gene in c("WHSC1","CACNA1H","JAK2","HTR3A","ABL1","SHROOM3","SETD2")){
  sector_names <- c(paste0(gene,"_D"),paste0(gene,"_D_R"),paste0(gene,"_R"))
  index <- c(circos_table_exp[1,paste0(gene,"_D")],
             circos_table_exp[2,paste0(gene,"_D_R")],
             circos_table_exp[3,paste0(gene,"_R")]) > 0
  sector_names <- sector_names[index]
  highlight.sector(sector_names, track.index = 1, col = NA, 
                   text = gene, cex = 1, text.col = "black", niceFacing = TRUE,text.vjust="8mm",xpd=TRUE)
}
# for(gene in c("NR3C1")){
#   sector_names <- c(paste0(gene,"_D"),paste0(gene,"_D_R"),paste0(gene,"_R"))
#   index <- c(circos_table_exp[1,paste0(gene,"_D")],
#              circos_table_exp[2,paste0(gene,"_D_R")],
#              circos_table_exp[3,paste0(gene,"_R")]) > 0
#   sector_names <- sector_names[index]
#   highlight.sector(sector_names, track.index = 1, col = NA, 
#                    text = gene, cex = 1, text.col = "black", niceFacing = TRUE,text.vjust="4mm",xpd=TRUE)
# }


dev.off()

circos.clear()


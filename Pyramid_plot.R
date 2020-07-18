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

sink("Pyramid_plot.txt")

pdf("Pyramid_plot.pdf",width=6,height=6,useDingbats = FALSE)
#par( mfrow = c( 1, 2 ) )
par( mar = c( 3,3,3,3 ) )

#Diagnosis

mut_samples <- list()
for(gene in mut_Genes){
  mut_samples[[paste0(gene,"_MUT")]] <- unique(ALL_mutations_D$Sample[ALL_mutations_D$Gene==gene])
}
for(gene in cnv_Genes){
  mut_samples[[paste0(gene,"_DEL")]] <- unlist(strsplit(CNV_Del_summary$Diagnosis[CNV_Del_summary$Gene==gene],split=","))
}

#pvalues <- c()
#for(i in 2:N){
#  for(j in 1:(i-1)){
#    gene_1_samples <- mut_samples[[i]]
#    gene_2_samples <- mut_samples[[j]]
#    if(length(gene_1_samples)>=3 & length(gene_2_samples)>=3){
#      pvalues <- c(pvalues,phyper(length(intersect(gene_1_samples,gene_2_samples)),length(gene_1_samples),175-length(gene_1_samples),length(gene_2_samples),lower.tail = TRUE))
#    }
#  }
#}

plot(0,
     xlim=c(1,N),
     ylim=c(0,N),
     type="n",
     xaxt="n",
     yaxt="n",
     bty="n",
     xlab=NA,
     ylab=NA)

for(i in 2:N){
  for(j in 1:(i-1)){
    gene_1_samples <- mut_samples[[i]]
    gene_2_samples <- mut_samples[[j]]
    if(length(gene_1_samples)>=5 & length(gene_2_samples)>=5){
      p1 <- phyper(length(intersect(gene_1_samples,gene_2_samples)),length(gene_1_samples),175-length(gene_1_samples),length(gene_2_samples),lower.tail = TRUE)
      p2 <- phyper(length(intersect(gene_1_samples,gene_2_samples)),length(gene_1_samples),175-length(gene_1_samples),length(gene_2_samples),lower.tail = FALSE)
      #print(paste(c(names(mut_samples)[i],names(mut_samples)[j],p1,p2)))
      if(any(c(p1,p2)<=0.05)){
        if(p1<=0.05){
          points(i,j,pch=16,cex=-log10(max(p1,0.000001))/2,col="red")
        }
        else{
          points(i,j,pch=16,cex=-log10(max(p2,0.000001))/2,col="green")
        }
      }
      else{
        points(i,j,pch=16,cex=0.5,col="grey")
      }
    }
    else{
      points(i,j,pch=16,cex=0.5,col="grey")
    }
  }
}
text(c(2:N)-0.5,c(1:(N-1))+0.5,labels=names(mut_samples)[2:N],xpd=TRUE,cex=0.5,pos=2,font=3,srt=-45)

#Relapse

mut_samples <- list()
for(gene in mut_Genes){
  mut_samples[[paste0(gene,"_MUT")]] <- unique(ALL_mutations_R$Sample[ALL_mutations_R$Gene==gene])
}
for(gene in cnv_Genes){
  mut_samples[[paste0(gene,"_DEL")]] <- unlist(strsplit(CNV_Del_summary$Relapse[CNV_Del_summary$Gene==gene],split=","))
}

plot(0,
     xlim=c(1,N),
     ylim=c(0,N),
     type="n",
     xaxt="n",
     yaxt="n",
     bty="n",
     xlab=NA,
     ylab=NA)

for(i in N:2){
  for(j in 1:(i-1)){
    gene_1_samples <- mut_samples[[i]]
    gene_2_samples <- mut_samples[[j]]
    if(length(gene_1_samples)>=5 & length(gene_2_samples)>=5){
      p1 <- phyper(length(intersect(gene_1_samples,gene_2_samples)),length(gene_1_samples),175-length(gene_1_samples),length(gene_2_samples),lower.tail = TRUE)
      p2 <- phyper(length(intersect(gene_1_samples,gene_2_samples)),length(gene_1_samples),175-length(gene_1_samples),length(gene_2_samples),lower.tail = FALSE)
      #print(paste(c(i,j,p1,p2)))
      if(any(c(p1,p2)<=0.05)){
        if(p1<=0.05){
          points(N-i+2,j,pch=16,cex=-log10(max(p1,0.000001))/2,col="red")
        }
        else{
          points(N-i+2,j,pch=16,cex=-log10(max(p2,0.000001))/2,col="green")
        }
      }
      else{
        points(N-i+2,j,pch=16,cex=0.5,col="grey")
      }
    }
    else{
      points(N-i+2,j,pch=16,cex=0.5,col="grey")
    }
  }
}
text(-1,c(1:(N-1)),labels=names(mut_samples)[1:(N-1)],xpd=TRUE,cex=0.5,adj=0.5,font=3)
text(c(N:2)+0.5,c(1:(N-1))+0.6,labels=names(mut_samples)[2:N],xpd=TRUE,cex=0.5,pos=4,font=3,srt=45)

dev.off()

sink()

###################
# pdf("legend.pdf",height=6,width=6,useDingbats = FALSE)
# 
# par(mar=c(0,0,0,0))
# 
# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# legend("topleft", legend =c('0.05', '0.01', '0.001'), 
#        pch=1, cex=-log10(0.05)/2, bty='n',
#        col = c('black'))
# #mtext("Alteration", at=0.2, cex=2)
# 
# dev.off()


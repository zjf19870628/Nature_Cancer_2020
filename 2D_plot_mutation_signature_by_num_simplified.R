# signature_contribution <- data.matrix(read.csv("Illumina_signature_contribution_combined.csv",
#                                                row.names=1,header=TRUE,check.names=FALSE))
# for(i in 1:ncol(signature_contribution)){
#   signature_contribution[,i] <- signature_contribution[,i]/sum(signature_contribution[,i]) * 100
# }
signature_contribution <- data.matrix(read.csv("Illumina_signature_contribution.csv",
                                               row.names=1,header=TRUE,check.names=FALSE))
signature_contribution <- signature_contribution[,sapply(colnames(signature_contribution),function(x) !(unlist(strsplit(x,split="_"))[3] %in% c("MH1574","MH3406","15342")))]
for(i in 1:ncol(signature_contribution)){
  signature_contribution[,i] <- signature_contribution[,i]/sum(signature_contribution[,i]) * 100
}

D_samples <- colnames(signature_contribution)[1:46]
R_samples <- colnames(signature_contribution)[47:91]

signature_contribution <- cbind(apply(signature_contribution[,1:46]>0,1,sum),
                                apply(signature_contribution[,47:91]>0,1,sum))

signature_qvalues <- data.matrix(read.csv("signature_qvalues.csv",row.names=1))

library(RColorBrewer)
mycol <- col2rgb(brewer.pal(9,"Set1")[c(2,1)])

plot_data <- list(D=signature_contribution[,1],
                  R=signature_contribution[,2])
plot_data <- as.data.frame(plot_data)
rownames(plot_data) <- rownames(signature_contribution)
for(i in 1:nrow(plot_data)){
  plot_data$color[i] <- rgb(t(as.matrix((mycol[,1] * plot_data$D[i]
                                         + mycol[,2] * plot_data$R[i])/sum(plot_data[i,1:2])/255,byrow=FALSE,ncol=3)))
  plot_data$cex[i] <- max(1,-log10(signature_qvalues[rownames(plot_data)[i],1]))
}

plot_data <- plot_data[signature_qvalues[,1] < 0.01,]

pdf("2D_plot_mutation_signature_by_num_simplified.pdf",width=5,height=5,useDingbats = FALSE)
par(mar=c(5,5,2,2))

plot(plot_data$D,
     plot_data$R,
     xlim=c(-0.1,45),
     ylim=c(-0.1,45),
     col="grey30",
     bg=plot_data$color,
     pch=21,
     bty="l",
     cex=plot_data$cex,
     xlab="Number of samples (diagnosis)",
     ylab="Number of samples (relapse)",
     cex.lab=1.5,
     xpd=TRUE,
     axes=FALSE)
axis(1,at=c(0,10,20,30,40),labels=c(0,10,20,30,40),col=brewer.pal(9,"Set1")[2],cex.axis=1.5)
axis(2,at=c(0,10,20,30,40),labels=c(0,10,20,30,40),col=brewer.pal(9,"Set1")[1],cex.axis=1.5)
lines(x = c(0,35), y = c(0,35),lty=2,col="grey")

dev.off()

###############################################

pdf("2D_plot_mutation_signature_labeled_by_num_simplified.pdf",width=5,height=5,useDingbats = FALSE)
par(mar=c(5,5,2,2))

plot(plot_data$D,
     plot_data$R,
     xlim=c(-0.1,45),
     ylim=c(-0.1,45),
     col="grey30",
     bg=plot_data$color,
     pch=21,
     bty="l",
     cex=plot_data$cex,
     xlab="Number of samples (diagnosis)",
     ylab="Number of samples (relapse)",
     cex.lab=1.5,
     xpd=TRUE,
     axes=FALSE)
axis(1,at=c(0,10,20,30,40),labels=c(0,10,20,30,40),col=brewer.pal(9,"Set1")[2],cex.axis=1.5)
axis(2,at=c(0,10,20,30,40),labels=c(0,10,20,30,40),col=brewer.pal(9,"Set1")[1],cex.axis=1.5)
lines(x = c(0,35), y = c(0,35),lty=2,col="grey")

text(plot_data$D,
     plot_data$R,
     labels = sapply(rownames(plot_data),function(x) unlist(strsplit(x,split="\\."))[2]),
     pos=3,xpd=TRUE)

dev.off()

#########################
# pdf("legend.pdf",height=5,width=5,useDingbats = FALSE)
# 
# par(mar=c(0,0,0,0))
# 
# plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
# legend("topleft", legend =c('0.01', '0.001', '0.0001','0.00001'), 
#        pch=1, cex=1, bty='n',
#        col = c('black'))
# #mtext("Alteration", at=0.2, cex=2)
# 
# dev.off()

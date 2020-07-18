mut_fre <- read.csv("../Recurrence_analysis/Recurrence_analysis_I3.csv")

library(RColorBrewer)
mycol <- col2rgb(brewer.pal(9,"Set1")[c(2,1,3)])

plot_data <- list(D=as.numeric(mut_fre$Number.of.diagnosis.specific.mutant.cases),
                  R=as.numeric(mut_fre$Number.of.relapse.specific.mutant.cases),
                  D_R=as.numeric(mut_fre$Number.of.Mutant.Cases) - (as.numeric(mut_fre$Number.of.diagnosis.specific.mutant.cases) + as.numeric(mut_fre$Number.of.relapse.specific.mutant.cases)))
plot_data <- as.data.frame(plot_data)
rownames(plot_data) <- mut_fre$Gene
for(i in 1:nrow(plot_data)){
  plot_data$color[i] <- rgb(t(as.matrix((mycol[,1] * plot_data$D[i]
                      + mycol[,2] * plot_data$R[i]
                      + mycol[,3] * plot_data$D_R[i])/sum(plot_data[i,1:3])/255,byrow=FALSE,ncol=3)))
}

# library("gg3D")
# 
# pdf("3D_plot.pdf",width=8,height=7)
# 
# ggplot(plot_data, aes(x=D, y=R, z=D_R,color=color)) + 
#   theme_void() +
#   axes_3D() +
#   stat_3D(size=10*as.numeric(mut_fre$Number.of.Mutant.Cases)/31) +
#   labs_3D(
#     labs=c("Diagnosis", "Relapse", "Diagnosis & Relapse"),
#     hjust=c(0,1,1), vjust=c(1, 1, -0.2), angle=c(0, 0, 0)) +
#   theme(legend.position = "none")
# 
# dev.off()
# 
# 
library("plot3D")
# library('rgl')
# library('car')
# library("scatterplot3d")
# 
# scatter3d(x = plot_data$D, 
#           y = plot_data$R , 
#           z = plot_data$D_R, 
#           point.col = "red",
#           surface= FALSE, 
#           xlab = 'Diagnosis',
#           ylab = 'Relapse',
#           zlab = 'Diagnosis & Relapse',
#           axis.col = c('firebrick3','black','goldenrod1'), 
#           axis.ticks = F,axis.scales = T,sphere.size =1.2,
#           theta=-40)
# 
# rgl.postscript("final_lung_expanded_mutation.pdf",fmt="pdf")

x0 <- c(0, 0, 0)
y0 <- c(10, 10, 10)
z0 <- c(0, 0, 0)
x1 <- c(0, 30, 0)
y1 <- c(0, 10, 10)
z1 <- c(0, 0, 15)
cols <- brewer.pal(9,"Set1")[c(2,1,3)]

pdf("3D_plot.pdf",width=7,height=7,useDingbats = FALSE)
par(mar=c(2,2,2,2))

arrows3D(x0, y0, z0, x1, y1, z1, col = cols,
         xlim=c(0,30),ylim=c(0,10),zlim=c(0,15),
         lwd = 2, d = 3, bty ="n",margin = c(0, 0))

for(i in 1:nrow(plot_data)){
  points3D(plot_data$R[i], 10-plot_data$D[i], plot_data$D_R[i], 
           add = TRUE, col=plot_data$color[i], border="black",
           colkey = FALSE, pch = 19, cex = max(1,sum(plot_data[i,1:3])/5))
}

text3D(plot_data["NRAS","R"] + 0.5, 10-plot_data["NRAS","D"] + 0.6, plot_data["NRAS","D_R"] + 0.5, "NRAS",
      col = plot_data["NRAS","color"], pos=3,add=TRUE, colkey = FALSE)

text3D(plot_data["KRAS","R"] + 0.5, 10-plot_data["KRAS","D"] - 1, plot_data["KRAS","D_R"] + 3, "KRAS",
       col = plot_data["KRAS","color"], pos=3,add=TRUE, colkey = FALSE)

text3D(plot_data["NT5C2","R"] + 0.5, 10-plot_data["NT5C2","D"] + 0.5, plot_data["NT5C2","D_R"] + 0.5, "NT5C2",
       col = plot_data["NT5C2","color"], pos=3,add=TRUE, colkey = FALSE)

text3D(plot_data["NOTCH1","R"] + 0.5, 10-plot_data["NOTCH1","D"] - 1, plot_data["NOTCH1","D_R"] + 2.5, "NOTCH1",
       col = plot_data["NOTCH1","color"], pos=3,add=TRUE, colkey = FALSE)

text3D(plot_data["TP53","R"] + 0.5, 10-plot_data["TP53","D"] + 0.5, plot_data["TP53","D_R"] + 0.5, "TP53",
       col = plot_data["TP53","color"], pos=3,add=TRUE, colkey = FALSE)

text3D(plot_data["CREBBP","R"] + 0.5, 10-plot_data["CREBBP","D"] + 0.5, plot_data["CREBBP","D_R"] + 0.5, "CREBBP",
       col = plot_data["CREBBP","color"], pos=3,add=TRUE, colkey = FALSE)

text3D(plot_data["WT1","R"], 10-plot_data["WT1","D"], plot_data["WT1","D_R"] + 1, "WT1",
       col = plot_data["WT1","color"], pos=3,add=TRUE, colkey = FALSE)

dev.off()


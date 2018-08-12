library(limma)
library(iChip)
library(geneplotter)
library(parallel)
options(digits=3)

set.seed(777)
targets <- readTargets("EA02QC.exp",sep="\t")
targets_lof<-targets[1:2,]        #Files of loss of function of WOX5
RGraw <- read.maimages(targets_lof, source="agilent", columns = list(G = "gMeanSignal", Gb = "gBGMedianSignal", R = "rMeanSignal", Rb = "rBGMedianSignal"), annotation = c("Start", "Sequence", "ProbeUID", "ControlType", "ProbeName", "GeneName", "SystematicName"))
#Set names to sample names
colnames(RGraw$G) <- targets_lof$GreenID
colnames(RGraw$Gb) <- targets_lof$GreenID
colnames(RGraw$R) <- targets_lof$RedID
colnames(RGraw$Rb) <- targets_lof$RedID
#Set Rownames 
#Plot diagnostics of the raw values.
nG <- ncol(RGraw$G)
nR <- ncol(RGraw$R)
allGR <- cbind(RGraw$G,RGraw$R)
allGRb <- cbind(RGraw$Gb,RGraw$Rb)
#Boxplots have options for common log2(y) limits, and vertical x-labels
boxplot(data.frame(log2(allGRb)),main="Raw Background Intensities",ylim=c(0,20),las=2, col=c(rep("green",nG), rep("red",nR)))
#Reds have a slightly higher (6 versus 5) bkgd.
savepng("BG_BN")
dev.off()
boxplot(data.frame(log2(allGR)),main="Raw Foreground Intensities",ylim=c(0,20),las=2, col=c(rep("green",nG), rep("red",nR)))
#Medians similar for all. Red slightly larger spread.
savepng("FG_BN")
dev.off()
#Cluster analysis of all samples against all.
#Generate correlation matrix for green and red
pcor <- cor(allGR, method="spearman")
#Generate coefficient of correlation
rsq <- pcor*pcor
#Save the coefficients of correlation
write.table(rsq, file="correlation_BN.txt", quote=FALSE, sep="\t", col.names=NA)
#Inputs cluster separately from the IPs.
plot(hclust(dist(1-pcor)), main="Spearman Correlation of Raw Intensities")
savepdf("dendro_BN")
#IP contrasts greatly from input.
heatmap(pcor, Rowv=NA, Colv=NA, symm=TRUE)
savepdf("heatmap_BN")
plotDensities(RGraw)
#Both red and green have peaks around 11 and 16. The green has a higher peak at 11, and the red has a higher peak at 16.
savepdf("density_BN")

#maplots
MAraw <- MA.RG(RGraw, bc.method="none")
#All plots have bananas that bend towards input.
plotMA3by2(MAraw, prefix="MA_BN")
#Not run, imageplots
#imageplot(MAraw$M[,1], layout=list(ngrid.r=1, ngrid.c=1, nspot.r=1068, nspot.c=912), zlim=c(-3,3))
#Unfortunately, 16574 spots are not accounted for in the FE text file, giving the above error.
#Get log intensities
RGlog <- RGraw
RGlog$G <- log2(RGraw$G)
RGlog$R <- log2(RGraw$R)
#Slight variance in the medians.
boxplot(MAraw$M, main="Raw Log-Ratios", ylim=c(-10,10), las=2)
savepng("Mbox_BN")
#Slight variance in the medians.
boxplot(MAraw$A, main="Raw Average Intensities", ylim=c(0,20), las=2)
savepng("Abox_BN")

#Normalization
RG <- RGraw
#Removed the control features before generating log ratios
i <- RG$genes$ControlType==0
#950368 probes remain
RGnc <- RG[i,]
#Try first the normexp BC method. Then check the next boxplot (FG_AN.png).
MA <- normalizeWithinArrays(RGnc, method="loess", bc.method="normexp", offset=32)
#Use the edwards BC method some arrays have compressed distribution of foreground.
#MA <- normalizeWithinArrays(RGnc, method="loess", bc.method="edwards", offset=32)
#Convert normalized MAlist back to RG intensities
RGnorm <- RG.MA(MA)
colnames(RGnorm$G) <- colnames(RGraw$G)
colnames(RGnorm$R) <- colnames(RGraw$R)

#Plot diagnostics of the normalized values.
allGRn <- cbind(RGnorm$G,RGnorm$R)
#Boxplots have options for common log2(y) limits, and vertical x-labels
boxplot(data.frame(log2(allGRn)),main="Normalized Foreground Intensities",ylim=c(0,20),las=2, col=c(rep("green",nG), rep("red",nR)))
#Slight variance in the medians.
savepng("FG_AN")

#Cluster analysis of all samples against all.
#Generate correlation matrix for green and red
pcorn <- cor(allGRn, method="spearman")
#Generate coefficient of correlation
rsqn <- pcorn*pcorn
#Save the coefficients of correlation
write.table(rsqn, file="correlation_AN.txt", quote=FALSE, sep="\t", col.names=NA)
#Same as before: inputs cluster separately from IPs.
plot(hclust(dist(1-pcorn)), main="Spearman Correlation of Normalized Intensities")
savepdf("dendro_AN")
#Inputs different from IPs.
heatmap(pcorn, Rowv=NA, Colv=NA, symm=TRUE)
savepdf("heatmap_AN")
#Fair red-green overlap at ca. 12 peak. Red has ca. 15 peak; green has ca. 16 peak.
plotDensities(RGnorm)
savepdf("density_AN")

#Arrays fairly straight, but seem to be slightly negative.
plotMA3by2(MA, prefix="MA_AN")
#Could not do imageplots, because I had taken out the control spots.
#Medians and spreads similar.
boxplot(MA$M, main="Normalized Log-Ratios", ylim=c(-10,10), las=2)
savepng("Mbox_AN")
#Slight variance in medians.
boxplot(MA$A, main="Normalized Average Intensities", ylim=c(0,20), las=2)
savepng("Abox_AN")

Systematic_names<-strsplit(RGnorm$genes$SystematicName,":")
chr<-do.call( rbind, Systematic_names)[,1]
start<-do.call(rbind,(strsplit((do.call(rbind, Systematic_names)[,2]),"-")))[,1]
Wox5<-cbind(chr,start,allGRn)
oct4<-data.frame(Wox5)
rownames(oct4)<-c()
oct4[3:6]<-lapply(oct4[3:6], as.numeric)
oct4 = oct4[order(oct4[,1],oct4[,2]),]

oct4lmt = lmtstat(oct4[,5:6],(oct4[,3:4]))
oct4Y = cbind(oct4[,1],oct4lmt)

oct4res2 = iChip2(Y=oct4Y,burnin=2000,sampling=10000,winsize=2,sdcut=2,beta=1.25,verbose=FALSE)

enrichreg(pos=oct4Y[,1:2],enrich=oct4lmt,pp=oct4res2$pp,cutoff=0.9,method="ppcut",maxgap=500)
###############################################For WOX5#################################
targets_wox5<-targets[3,]        #Files of loss of function of WOX5
RGraw_wox5 <- read.maimages(targets_wox5, source="agilent", columns = list(G = "gMeanSignal", Gb = "gBGMedianSignal", R = "rMeanSignal", Rb = "rBGMedianSignal"), annotation = c("Start", "Sequence", "ProbeUID", "ControlType", "ProbeName", "GeneName", "SystematicName"))
#Set names to sample names
colnames(RGraw_wox5$G) <- targets_wox5$GreenID
colnames(RGraw_wox5$Gb) <- targets_wox5$GreenID
colnames(RGraw_wox5$R) <- targets_wox5$RedID
colnames(RGraw_wox5$Rb) <- targets_wox5$RedID
#Set Rownames 
#Plot diagnostics of the raw values.
nG_wox5 <- ncol(RGraw_wox5$G)
nR_wox5 <- ncol(RGraw_wox5$R)
allGR_wox5 <- cbind(RGraw_wox5$G,RGraw_wox5$R)
allGRb_wox5 <- cbind(RGraw_wox5$Gb,RGraw_wox5$Rb)
#Boxplots have options for common log2(y) limits, and vertical x-labels
boxplot(data.frame(log2(allGRb_wox5)),main="Raw Background Intensities",ylim=c(0,20),las=2, col=c(rep("green",nG_wox5), rep("red",nR_wox5)))
#Reds have a slightly higher (6 versus 5) bkgd.
savepng("BG_BN")
dev.off()
boxplot(data.frame(log2(allGR_wox5)),main="Raw Foreground Intensities",ylim=c(0,20),las=2, col=c(rep("green",nG_wox5), rep("red",nR_wox5)))
#Medians similar for all. Red slightly larger spread.
savepng("FG_BN")
dev.off()
#Cluster analysis of all samples against all.
#Generate correlation matrix for green and red
pcor_wox5 <- cor(allGR_wox5, method="spearman")
#Generate coefficient of correlation
rsq_wox5 <- pcor_wox5*pcor_wox5
#Save the coefficients of correlation
write.table(rsq_wox5, file="correlation_BN_wox5x5.txt", quote=FALSE, sep="\t", col.names=NA)
#Inputs cluster separately from the IPs.
plot(hclust(dist(1-pcor_wox5)), main="Spearman Correlation of Raw Intensities")
savepdf("dendro_BN")
#IP contrasts greatly from input.
heatmap(pcor_wox5, Rowv=NA, Colv=NA, symm=TRUE)
savepdf("heatmap_BN")
plotDensities(RGraw_wox5)
#Both red and green have peaks around 11 and 16. The green has a higher peak at 11, and the red has a higher peak at 16.
savepdf("density_BN")

#maplots
MAraw_wox5 <- MA.RG(RGraw_wox5, bc.method="none")
#All plots have bananas that bend towards input.
plotMA3by2(MAraw_wox5, prefix="MA_BN_wox5")
#Not run, imageplots
#imageplot(MAraw$M[,1], layout=list(ngrid.r=1, ngrid.c=1, nspot.r=1068, nspot.c=912), zlim=c(-3,3))
#Unfortunately, 16574 spots are not accounted for in the FE text file, giving the above error.
#Get log intensities
RGlog_wox5 <- RGraw_wox5
RGlog_wox5$G <- log2(RGraw_wox5$G)
RGlog_wox5$R <- log2(RGraw_wox5$R)
#Slight variance in the medians.
boxplot(MAraw_wox5$M, main="Raw Log-Ratios", ylim=c(-10,10), las=2)
savepng("Mbox_BN")
#Slight variance in the medians.
boxplot(MAraw_wox5$A, main="Raw Average Intensities", ylim=c(0,20), las=2)
savepng("Abox_BN")

#Normalization
RG_wox5 <- RGraw_wox5
#Removed the control features before generating log ratios
i_wox5 <- RG_wox5$genes$ControlType==0
#950368 probes remain
RGnc_wox5 <- RG_wox5[i_wox5,]
#Try first the normexp BC method. Then check the next boxplot (FG_AN.png).
MA_wox5 <- normalizeWithinArrays(RGnc_wox5, method="loess", bc.method="normexp", offset=32)
#Use the edwards BC method some arrays have compressed distribution of foreground.
#MA <- normalizeWithinArrays(RGnc, method="loess", bc.method="edwards", offset=32)
#Convert normalized MAlist back to RG intensities
RGnorm_wox5 <- RG.MA(MA_wox5)
colnames(RGnorm_wox5$G) <- colnames(RGraw_wox5$G)
colnames(RGnorm_wox5$R) <- colnames(RGraw_wox5$R)

#Plot diagnostics of the normalized values.
allGRn_wox5 <- cbind(RGnorm_wox5$G,RGnorm_wox5$R)
#Boxplots have options for common log2(y) limits, and vertical x-labels
boxplot(data.frame(log2(allGRn_wox5)),main="Normalized Foreground Intensities",ylim=c(0,20),las=2, col=c(rep("green",nG_wox5), rep("red",nR_wox5)))
#Slight variance in the medians.
savepng("FG_AN")

#Cluster analysis of all samples against all.
#Generate correlation matrix for green and red
pcorn_wox5 <- cor(allGRn_wox5, method="spearman")
#Generate coefficient of correlation
rsqn_wox5 <- pcorn_wox5*pcorn_wox5
#Save the coefficients of correlation
write.table(rsqn_wox5, file="correlation_AN.txt", quote=FALSE, sep="\t", col.names=NA)
#Same as before: inputs cluster separately from IPs.
plot(hclust(dist(1-pcorn_wox5)), main="Spearman Correlation of Normalized Intensities")
savepdf("dendro_AN")
#Inputs different from IPs.
heatmap(pcorn_wox5, Rowv=NA, Colv=NA, symm=TRUE)
savepdf("heatmap_AN")
#Fair red-green overlap at ca. 12 peak. Red has ca. 15 peak; green has ca. 16 peak.
plotDensities(RGnorm_wox5)
savepdf("density_AN")

#Arrays fairly straight, but seem to be slightly negative.
plotMA3by2(MA_wox5, prefix="MA_AN")
#Could not do imageplots, because I had taken out the control spots.
#Medians and spreads similar.
boxplot(MA_wox5$M, main="Normalized Log-Ratios", ylim=c(-10,10), las=2)
savepng("Mbox_AN_wox5")
#Slight variance in medians.
boxplot(MA_wox5$A, main="Normalized Average Intensities", ylim=c(0,20), las=2)
savepng("Abox_AN")

Systematic_names_wox5<-strsplit(RGnorm_wox5$genes$SystematicName,":")
chr_wox5<-do.call( rbind, Systematic_names_wox5)[,1]
start_wox5<-do.call(rbind,(strsplit((do.call(rbind, Systematic_names_wox5)[,2]),"-")))[,1]
Wox5_wox5<-cbind(chr_wox5,start_wox5,allGRn_wox5)



oct4_wox5<-data.frame(Wox5_wox5)
rownames(oct4_wox5)<-c()
oct4_wox5[3:4]<-lapply(oct4_wox5[3:4], as.numeric)
oct4_wox5 = oct4_wox5[order(oct4_wox5[,1],oct4_wox5[,2]),]

oct4lmt_wox5 = lmtstat((oct4_wox5[,4]),(oct4_wox5[,3]))
oct4Y_wox5 = cbind(oct4_wox5[,1],oct4lmt_wox5)

oct4Y_wox5 = cbind(chr_wox5,rownames(allGRn_wox5),oct4lmt_wox5)

oct4res2_wox5 = iChip2(Y=oct4Y_wox5,burnin=2000,sampling=10000,winsize=2,sdcut=2,beta=1.25,verbose=FALSE)

enrichreg(pos=oct4Y_wox5[,1:2],enrich=oct4lmt_wox5,pp=oct4res2_wox5$pp,cutoff=0.9,method="ppcut",maxgap=500)

library(limma)
library(dplyr)
require(Biobase) 
library(geneplotter)
library(survcomp)
library(snow)
library(fabia)
library(biclust)
options(digits=3)
targets_single_channel<-readTargets("expressionList.exp")   #reading single channel replicates
targets_two_channels <- readTargets("two_channel.exp")     #reading Two channels replicates
#Read The raw data
targets2 <- targetsA2C(targets_two_channels)                #Convert a two-color targets dataframe with one row per array to one with one row per channel
targets_double_channel<- subset(targets2,targets2$Target =="WT")  #Select WT replicates
targets_two_channels$RedID<-NULL
targets_all=bind_rows(targets_single_channel,targets_double_channel)   #Bind all replicates from single channel and two channels in one file
targets_all$Cy3[13:18]<-"QC"                                           #Replace name WT to QC             
targets_all<-targets_all[,-6]                                          #Select certain columns 
targets_all$Array[15:18]<-c(15,16,17,18)                               #Assign Array numbers
#####################################Analysis_of_QC VS CRC#######################################################################################################################################
eset_target<- read.maimages(targets_all,source="agilent", columns = list(E = "gMeanSignal", Eb = "gBGMedianSignal"), annotation = c("ControlType","Probe ID", "ProbeName","SystematicName")) #Read the selected replicates with annotations 
colnames(eset_target)<-eset_target$targets$GreenID
colnames(eset_target)
rownames(eset_target)<- eset_target$genes$SystematicName
boxplot(data.frame(log2(eset_target$E)),main="Green foreground",ylim=c(0,20),las=2)
#savepdf("boxplot_A1")
savepdf("All_replicates/boxplot_A1")
eset_nc <-eset_target[eset_target$genes$ControlType ==0,]
boxplot(data.frame(log2(eset_nc$E)),main="Green foreground After Quality Control",ylim=c(0,20),las=2)
#savepdf("boxplot_Quality_Control_A1")
savepdf("All_replicates/boxplot_A1_After_Quality_Control")
#Choose edwards  methods because its after applying this method its distribution is more similar to real distribution############
EG_edwards <- backgroundCorrect(eset_nc,method="edwards", offset=32)
boxplot(data.frame(log2(EG_edwards$E)),main="After edwards Background",ylim=c(0,20),las=2)
savepdf("All_replicates/boxplot_edwards_A1")
EG_quantile<- normalizeBetweenArrays(EG_edwards, method="quantile")
boxplot(data.frame(EG_quantile$E),main="After quantile normalization",ylim=c(0,20),las=2)
savepdf("All_replicates/boxplot_quantile_A1")
EG_cyclicloess<- normalizeBetweenArrays(EG_edwards, method="cyclicloess")
boxplot(data.frame(EG_cyclicloess$E),main="After cyclicloess Normalization",ylim=c(0,20),las=2)
savepdf("All_replicates/boxplot_cyclicloess_A1")
###################################edwards################################
#EG_normexp <- backgroundCorrect(eset_nc,method="normexp", offset=32)
#boxplot(data.frame(log2(EG_normexp$E)),main="After normexp Background",ylim=c(0,20),las=2)
#savepdf("boxplot_normexp_A1")
#savepdf(paste(var_foldername,"boxplot_normexp_A1",sep="/"))
###################################Hierarchical clustering##############################
EG_avg<-avereps(EG_cyclicloess)           #Microarray data object so that values for within-array replicate probes are replaced with their average.
boxplot(data.frame(EG_avg$E),main="Average Over Irregular Replicate Probes",ylim=c(0,20),las=2)
savepdf("All_replicates/Average_Irregular_Replicate")

pcor <- cor(eset_target$E, method="spearman")
heatmap(pcor, Rowv=NA, Colv=NA, symm=TRUE)
plot(hclust(dist(1-pcor)), main="Spearman Correlation of Raw Intensities")
pcor2 <- cor(EG_avg$E, method="spearman")
heatmap(pcor2, Rowv=NA, Colv=NA, symm=TRUE)
plot(hclust(dist(1-pcor2)), main="Spearman Correlation of Raw Intensities")
##########################Design###########################
f <- factor(EG_avg$targets$Cy3, levels=c("QC","CSC","CRC"))
design <- model.matrix(~0+f)
colnames(design)<-c("QC","CSC","CRC")
fit <- lmFit(EG_avg, design)
contrast.matrix <- makeContrasts("QCvsCSC"= QC-CSC,"QCvsCRC"=QC-CRC,"CSCvsQC"=CSC-QC,"CSCvsCRC"=CSC-CRC,"CRCvsQC"=CRC-QC,"CRCvsCSC"=CRC-CSC, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2  <- eBayes(fit2)
lQCvsCSC<-topTable(fit2, coef="QCvsCSC", adjust="BH",n=Inf) #Result for QC differs from QC
lQCvsCRC<-topTable(fit2, coef="QCvsCRC", adjust="BH",n=Inf) #Result QC differ from CRC
lCSCvsQC<-topTable(fit2, coef="CSCvsQC", adjust="BH",n=Inf) #Result CSC differs from CRC
lCSCvsCRC<-topTable(fit2, coef="CSCvsCRC", adjust="BH",n=Inf) #Result CRC differs from QC
lCRCvsQC<-topTable(fit2, coef="CRCvsQC", adjust="BH",n=Inf)  # Result CRC differs from CSC
lCRCvsCSC<-topTable(fit2, coef="CRCvsCSC", adjust="BH",n=Inf) #Result CRC differs from QC

write.table(lQCvsCSC[,c("SystematicName","logFC","adj.P.Val")],"complete_gene_lists/QCvsCSC.txt",row.names=FALSE,sep = "\t")
write.table(lQCvsCRC[,c("SystematicName","logFC","adj.P.Val")],"complete_gene_lists/QCvsCRC.txt",row.names=FALSE,sep = "\t")
write.table(lCSCvsQC[,c("SystematicName","logFC","adj.P.Val")],"complete_gene_lists/CSCvsQC.txt",row.names=FALSE,sep = "\t")
write.table(lCSCvsCRC[,c("SystematicName","logFC","adj.P.Val")],"complete_gene_lists/CSCvsCRC.txt",row.names=FALSE,sep = "\t")
write.table(lCRCvsQC[,c("SystematicName","logFC","adj.P.Val")],"complete_gene_lists/CRCvsQC.txt",row.names=FALSE,sep = "\t")
write.table(lCRCvsCSC[,c("SystematicName","logFC","adj.P.Val")],"complete_gene_lists/CRCvsCSC.txt",row.names=FALSE,sep = "\t")
##################################Gredient Lists##################################################################

QCvsCSC_lfcPos1<-lQCvsCSC[(lQCvsCSC$adj.P.Val<0.05)&(lQCvsCSC$logFC>1),] # cutoff  of lQCvsCSC into QCvsCSC_lfcPos1
CSCvsCRC_lfcPos1<-lCSCvsCRC[(lCSCvsCRC$adj.P.Val<0.05)&(lCSCvsCRC$logFC>1),]# cutoff of CSCvsCRC into CSCvsCRC_lfcPos1
QCvsCSC_lfcNeg1<-lQCvsCSC[(lQCvsCSC$adj.P.Val<0.05)&(lQCvsCSC$logFC< -1),] #cutoff of lQCvsCRC into QCvsCSC_lfcNeg1
CSCvsCRC_lfcNeg1<-lCSCvsCRC[(lCSCvsCRC$adj.P.Val<0.05)&(lCSCvsCRC$logFC< -1),] #cutoff of CSCvsCRC into CSCvsCRC_lfcNeg1


#write.table(QCvsCSC_lfcPos1[,c("SystematicName","logFC","adj.P.Val")],"Genes_lists_feb_22/QCvsCSC_lfcPos1.txt",row.names=FALSE,sep = "\t")
#write.table(CSCvsCRC_lfcPos1[,c("SystematicName","logFC","adj.P.Val")],"Genes_lists_feb_22/CSCvsCRC_lfcPos1.txt",row.names=FALSE,sep = "\t")
#write.table(QCvsCSC_lfcNeg1[,c("SystematicName","logFC","adj.P.Val")],"Genes_lists_feb_22/QCvsCSC_lfcNeg1.txt",row.names=FALSE,sep = "\t")
#write.table(CSCvsCRC_lfcNeg1[,c("SystematicName","logFC","adj.P.Val")],"Genes_lists_feb_22/QCvsCSC_lfcNeg1.txt",row.names=FALSE,sep = "\t")
##########################0.01 cutoff#################################
QCvsCSC_lfcPos1_0.01<-lQCvsCSC[(lQCvsCSC$adj.P.Val<0.01&lQCvsCSC$logFC>1),]   # cutoff  of lQCvsCSC into QCvsCSC_lfcPos1
CSCvsCRC_lfcPos1_0.01<-lCSCvsCRC[(lCSCvsCRC$adj.P.Val<0.01&lCSCvsCRC$logFC>1),]# cutoff of CSCvsCRC into CSCvsCRC_lfcPos1
QCvsCSC_lfcNeg1_0.01<-lQCvsCSC[(lQCvsCSC$adj.P.Val<0.01&lQCvsCSC$logFC< -1),]#cutoff of lCSCvsCRC into QCvsCSC_lfcNeg1
CSCvsCRC_lfcNeg1_0.01<-lCSCvsCRC[(lCSCvsCRC$adj.P.Val<0.01&lCSCvsCRC$logFC< -1),] #cutoff of CSCvsCRC into CSCvsCRC_lfcNeg1

#write.table(QCvsCSC_lfcPos1_0.01[,c("SystematicName","logFC","adj.P.Val")],"Genes_lists_feb_22/QCvsCSC_lfcPos1_0.01.txt",row.names=FALSE,sep = "\t")
#write.table(CSCvsCRC_lfcPos1_0.01[,c("SystematicName","logFC","adj.P.Val")],"Genes_lists_feb_22/CSCvsCRC_lfcPos1_0.01.txt",row.names=FALSE,sep = "\t")
#write.table(QCvsCSC_lfcNeg1_0.01[,c("SystematicName","logFC","adj.P.Val")],"Genes_lists_feb_22/QCvsCSC_lfcNeg1_0.01.txt",row.names=FALSE,sep = "\t")
#write.table(QCvsCSC_lfcNeg1_0.01[,c("SystematicName","logFC","adj.P.Val")],"Genes_lists_feb_22/QCvsCSC_lfcNeg1_0.01.txt",row.names=FALSE,sep = "\t")

####################gradient ASC###########################


incAGInum<-intersect(QCvsCSC_lfcPos1$SystematicName,CSCvsCRC_lfcPos1$SystematicName) #intersection of systematic names of QCvsCSC_lfcPos1 and CSCvsCRC_lfcPos1
QCvsCSC_lfcPos1_indices<-match(incAGInum,QCvsCSC_lfcPos1$SystematicName)  #get Indices of QCvsCSC_lfcPos1 of common systematic name of both lists
CSCvsCRC_lfcPos1_indices<-match(incAGInum,CSCvsCRC_lfcPos1$SystematicName) #get indices of CSCvsCRC_lfcPos1 of common systematic names of both lists

QCvsCSC_lfcPos1_data.frame<-QCvsCSC_lfcPos1[QCvsCSC_lfcPos1_indices,] #get the rows using indices 
QCvsCRC_lfcPos1_data.frame<-CSCvsCRC_lfcPos1[CSCvsCRC_lfcPos1_indices,] #get the rows using the indices
incL<-merge(QCvsCSC_lfcPos1_data.frame,QCvsCRC_lfcPos1_data.frame,by='SystematicName',all = TRUE) #merge the both list by systematic names
incL<-incL[,c("SystematicName","logFC.x","logFC.y")] #Select columns of Systematic name ,logFC of QCvsCSC_lfcPos1_data_frame and QCvsCRC_lfcPos1
colnames(incL)[2]<-"QCvsCSC_lfcPos1_LogFC" #Change column name of logFC.x
colnames(incL)[3]<-"QCvsCRC_lfcPos1_LogFC" #Change column name of logFC.y
incL$log2FCs<-incL$QCvsCSC_lfcPos1_LogFC+incL$QCvsCRC_lfcPos1_LogFC     #Add LogFc from both lists
incL<-incL[order(incL$log2FCs,decreasing = TRUE),]  #Sort the sum of two log2FCS by decreasing order
write.table(incL,"gradient/increasing_logFC.txt",row.names=FALSE,sep = "\t")
##################Gradient DSC#########################
decAGInum<-intersect(QCvsCSC_lfcNeg1$SystematicName,CSCvsCRC_lfcNeg1$SystematicName) #interaction of systematic names 
QCvsCSC_lfcNeg1_indices<-match(decAGInum,QCvsCSC_lfcNeg1$SystematicName) #get Indices of QCvsCSC_lfcNeg1 of common systematic name of both lists
CSCvsCRC_lfcNeg1_indices<-match(decAGInum,CSCvsCRC_lfcNeg1$SystematicName) #get Indices of CSCvsCRC_lfcNeg1 of common systematic name of both lists
QCvsCSC_lfcNeg1_data.frame<-QCvsCSC_lfcNeg1[QCvsCSC_lfcNeg1_indices,] #get the rows using the indices
QCvsCRC_lfcNeg1_data.frame<-CSCvsCRC_lfcNeg1[CSCvsCRC_lfcNeg1_indices,] #get the rows using the indices
dincL<-merge(QCvsCSC_lfcNeg1_data.frame,QCvsCRC_lfcNeg1_data.frame,by='SystematicName',all = TRUE) #merge the both list by systematic names
dincL<-dincL[,c("SystematicName","logFC.x","logFC.y")]#Select columns of Systematic name ,logFC of QCvsCSC_lfcNeg1_data_frame and QCvsCRC_lfcNeg1_data_frame
colnames(dincL)[2]<-"QCvsCSC_lfcNeg1_LogFC" #Change column name of logFC.x
colnames(dincL)[3]<-"CSCvsCRC_lfcNeg1_LogFC" #Change column name of logFC.y
dincL$log2FCs<-dincL$QCvsCSC_lfcNeg1_LogFC+dincL$CSCvsCRC_lfcNeg1_LogFC #Sum LogFcs of both lists 
dincL<-dincL[order(dincL$log2FCs,decreasing = FALSE),]  #Sort the sum of two log2FCS by increasing order
write.table(dincL,"gradient/decreasing_logFC.txt",row.names=FALSE,sep = "\t")



#################Function for cell-Type specific Lists#########################
cell_cutoff<-function(data.frame1,data.frame2,cutoff){
  data.frame1<-data.frame1[order(data.frame1$ProbeName),] #sort the dataframe by increasin order
  data.frame2<-data.frame2[order(data.frame2$ProbeName),] #sort the dataframe by increasing order
  cell_data.frame<-NA                                     #cell_data.frame
  cell_data.frame<-merge(data.frame1,data.frame2,by='ProbeName',all = TRUE)  #merge the both frame by Probe names 
  cell_data.frame$logFC.x<-(cell_data.frame$logFC.x-mean(cell_data.frame$logFC.x))/sd(cell_data.frame$logFC.x)
  cell_data.frame$logFC.y<-(cell_data.frame$logFC.y-mean(cell_data.frame$logFC.y))/sd(cell_data.frame$logFC.y)
  cell_data.frame$mPVal<-NA                                                  #initiate the new column for combine P vlaues
  cell_data.frame$mZ<-NA                                                     #Initiate the new column for combine zscore
  for (x in seq_along(cell_data.frame$ProbeName)) { #iterate oveer merger datafram by proble names and combine P value and Zscore by respectively Staufer method and fisher method
    cell_data.frame$mZ[x]<-(cell_data.frame$logFC.x[x]+cell_data.frame$logFC.y[x])/sqrt(2)
    cell_data.frame$mPVal[x]<-combine.test(p=c(cell_data.frame$adj.P.Val.x[x],cell_data.frame$adj.P.Val.y[x]),method ="fisher", na.rm=TRUE)
  }
  cell_data.frame<-cell_data.frame[(cell_data.frame$adj.P.Val.x<cutoff)&(cell_data.frame$adj.P.Val.y<cutoff),] #cutoff P value of both of indivisual dataset
  cell_data.frame<-cell_data.frame[cell_data.frame$mPVal<cutoff,]        #cutoff combined P value                                      
  cell_data.frame<-cell_data.frame[,c("SystematicName.x","mPVal","mZ")] #Select three columns from dataframe
  colnames(cell_data.frame)[1]<-"SystematicName"                          #Change the name of SytematicName.x to Systematic Name
  return(cell_data.frame)                                                #return the cell_data.frame
}

##################QC#########################
l05QC<-cell_cutoff(lQCvsCSC,lQCvsCRC,0.05)                #use the function for cutoff value 0.05 for QC
l01QC<-cell_cutoff(lQCvsCSC,lQCvsCRC,0.01)                #use the function for cutoff value 0.01 for QC
##################CSC######################## 
l05CSC<-cell_cutoff(lCSCvsQC,lCSCvsCRC,0.05)              #use the function for cutoff value 0.01 for CRC
l01CSC<-cell_cutoff(lCSCvsQC,lCSCvsCRC,0.01)              #use the function for cutoff value 0.01 for CRC
#################CRC#########################
l05CRC<-cell_cutoff(lCRCvsQC,lCRCvsCSC,0.05)              #use the function for cutoff value 0.01 for CSC
l01CRC<-cell_cutoff(lCRCvsQC,lCRCvsCSC,0.01)              #use the function for cutoff value 0.01 for CSC
##################Save Files#################
write.table(l05QC,"cell-type-specific_Lists/QC_005_cutoff.txt",row.names=FALSE,sep = "\t")
write.table(l01QC,"cell-type-specific_Lists/QC_001_cutoff.txt",row.names=FALSE,sep = "\t")
write.table(l05CSC,"cell-type-specific_Lists/CSC_005_cutoff.txt",row.names=FALSE,sep = "\t")
write.table(l01CSC,"cell-type-specific_Lists/CSC_001_cutoff.txt",row.names=FALSE,sep = "\t")
write.table(l05CRC,"cell-type-specific_Lists/CRC_005_cutoff.txt",row.names=FALSE,sep = "\t")
write.table(l01CRC,"cell-type-specific_Lists/CRC_001_cutoff.txt",row.names=FALSE,sep = "\t")
########################files of two lists############
l05QC_mlog2FC<-l05QC[(l05QC$mZ< -2)|(l05QC$mZ>2),]    #cutoff list of genes by mZ cutoff by -2 and 2
l05CSC_mlog2FC<-l05CSC[(l05CSC$mZ< -2)|(l05CSC$mZ>2),] #cutoff list of genes by mZ cutoff by -2 and 2
l05CRC_mlog2FC<-l05CRC[(l05CRC$mZ< -2)|(l05CRC$mZ>2),]#cutoff list of genes by mZ cutoff by -2 and 2
#################Save The Files########################################
write.table(l05QC_mlog2FC,"cell_specific_logFC_cutoff/QC_logFC_cutoff.txt",row.names=FALSE,sep = "\t")
write.table(l05CSC_mlog2FC,"cell_specific_logFC_cutoff/CSC_logFC_cutoff.txt",row.names=FALSE,sep = "\t")
write.table(l05CRC_mlog2FC,"cell_specific_logFC_cutoff/CRC_logFC_cutoff.txt",row.names=FALSE,sep = "\t")
################################Marker Genes by Ernst########################
QC_genes<-read.table("confirmed_genes/QC.txt",stringsAsFactors = FALSE)    #Read files of QC Maker Genes
CRC_genes<-read.table("confirmed_genes/CRC.txt",stringsAsFactors = FALSE)  #Read files of CRC Maker Genes
CSC_genes<-read.table("confirmed_genes/CSC.txt",stringsAsFactors = FALSE)  #Read files of CSC Maker Genes
QC_genes$V1<-substr(QC_genes$V1,1,9)                                       #Select first 9 character of QC Marker Genes
CRC_genes$V1<-substr(CRC_genes$V1,1,9)                                     #Select first 9 character of CRC Marker Genes
CSC_genes$V1<-substr(CSC_genes$V1,1,9)                                     #Select first 9 character of CSC Marker Genes

####################QC Adjusted P value with cutoff of 0.5#################
QC_genes$V1 %in% substr(l05QC$SystematicName,1,9)                         #Query QC marker Genes in QC Genes list with cutoff of P value of 0.5
CSC_genes$V1 %in% substr(l05QC$SystematicName,1,9)                        #Query CSC marker Genes in QC Genes list with cutoff of P value of 0.5
CRC_genes$V1 %in% substr(l05QC$SystematicName,1,9)                        #Query CRC marker Genes in QC Genes list with cutoff of P value of 0.5
####################CSC Adjusted P value with cutoff of 0.5#################
QC_genes$V1 %in% substr(l05CSC$SystematicName,1,9)                        #Query QC marker Genes in CSC Genes list with cutoff of P value of 0.5
CSC_genes$V1 %in% substr(l05CSC$SystematicName,1,9)                       #Query CSC marker Genes in CSC Genes list with cutoff of P value of 0.5
CRC_genes$V1 %in% substr(l05CSC$SystematicName,1,9)                       #Query CRC marker Genes in CSC Genes list with cutoff of P value of 0.5
####################CRC Adjusted P value with cutoff of 0.5#################
QC_genes$V1 %in% substr(l05CRC$SystematicName,1,9)                        #Query QC marker Genes in CRC Genes list with cutoff of P value of 0.5
CSC_genes$V1 %in% substr(l05CRC$SystematicName,1,9)                       #Query CSC marker Genes in CRC Genes list with cutoff of P value of 0.5
CRC_genes$V1 %in% substr(l05CRC$SystematicName,1,9)                       #Query CRC marker Genes in CRC Genes list with cutoff of P value of 0.5
####################QC Adjusted P value with cutoff of 0.1#################
QC_genes$V1 %in% substr(l01QC$SystematicName,1,9)                         #Query QC marker Genes in QC Genes list with cutoff of P value of 0.1
CSC_genes$V1 %in% substr(l01QC$SystematicName,1,9)                        #Query CSC marker Genes in QC Genes list with cutoff of P value of 0.1
CRC_genes$V1 %in% substr(l01QC$SystematicName,1,9)                        #Query CRC marker Genes in QC Genes list with cutoff of P value of 0.1
####################CSC Adjusted P value with cutoff of 0.1#################
QC_genes$V1 %in% substr(l01CSC$SystematicName,1,9)                        #Query QC marker Genes in CSC Genes list with cutoff of P value of 0.1
CSC_genes$V1 %in% substr(l01CSC$SystematicName,1,9)                       #Query CSC marker Genes in CSC Genes list with cutoff of P value of 0.1
CRC_genes$V1 %in% substr(l01CSC$SystematicName,1,9)                       #Query CRC marker Genes in CSC Genes list with cutoff of P value of 0.1  
####################CRC Adjusted P value with cutoff of 0.1#################
QC_genes$V1 %in% substr(l01CRC$SystematicName,1,9)                         #Query QC marker Genes in CRC Genes list with cutoff of P value of 0.1
CSC_genes$V1 %in% substr(l01CRC$SystematicName,1,9)                        #Query CSC marker Genes in CRC Genes list with cutoff of P value of 0.1
CRC_genes$V1 %in% substr(l01CRC$SystematicName,1,9)                       #Query CRC marker Genes in CRC Genes list with cutoff of P value of 0.1
##################ernst_list############################################
ernst05QC<-intersect(QC_genes$V1,substr(l05QC$SystematicName,1,9))         #Validation of  QC Ernst List and assign the marker Genes which are in the list with cutoff 0.05               
ernst05CSC<-intersect(CSC_genes$V1,substr(l05CSC$SystematicName,1,9))      #Validation of  CSC Ernst List and assign the marker Genes which are in the list with cutoff 0.05                 
ernst05CRC<-intersect(CRC_genes$V1,substr(l05CRC$SystematicName,1,9))     #Validation of  CRC Ernst List and assign the marker Genes which are in the list  with cutoff 0.05                 
ernst01QC<-intersect(QC_genes$V1,substr(l01QC$SystematicName,1,9))        #Validation of  QC Ernst List and assign the marker Genes which are in the list  with cutoff 0.01                
ernst01CSC<-intersect(CSC_genes$V1,substr(l01CSC$SystematicName,1,9))     #Validation of  CSC Ernst List and assign the marker Genes which are in the list with cutoff 0.01
ernst01CRC<-intersect(CRC_genes$V1,substr(l01CRC$SystematicName,1,9))     #Validation of  CRC Ernst List and assign the marker Genes which are in the list with cutoff 0.01
###################################
QC_genes$V1<-ernst05QC
CSC_genes$V1<-ernst05CSC                       
CRC_genes$V1<-ernst05CRC

####################################Save the Validated Marker Genes Lists#####################
write.table(ernst05QC,"ernstlists/ernst05QC.txt",row.names=FALSE,sep = "\t")
write.table(ernst05CSC,"ernstlists/ernst05CSC.txt",row.names=FALSE,sep = "\t")
write.table(ernst05CRC,"ernstlists/ernst05CRC.txt",row.names=FALSE,sep = "\t")
write.table(ernst01QC,"ernstlists/ernst01QC.txt",row.names=FALSE,sep = "\t")
write.table(ernst01CSC,"ernstlists/ernst01CSC.txt",row.names=FALSE,sep = "\t")
write.table(ernst01CRC,"ernstlists/ernst01CRC.txt",row.names=FALSE,sep = "\t")
####################Analysis##############################
CRC_cluster<-read.table("Data_Analysis/CRC_cluster.txt",stringsAsFactors=FALSE)  
QC_cluster<-read.table("Data_Analysis/QC_cluster.txt",stringsAsFactors=FALSE)
CSC_cluster<-read.table("Data_Analysis/CSC_cluster.txt",stringsAsFactors=FALSE)

####################################################Combining lists of 0.05########################
all_systematicNames<-unique(c(l05QC$SystematicName,l05CRC$SystematicName,l05CSC$SystematicName)) #combining the genes lists cutoff by 0.05 of three cells for cluster Analysis
indices_of_systematicNames<-match(all_systematicNames,rownames(EG_avg))    #Indices of Systematic Names of Average Expression List
data_for_clustering<-EG_avg[indices_of_systematicNames,]                   #Selection of rows on basis of Average Expression
####################################################Combining lists of 0.01########################
all_systematicNames<-unique(c(l01QC$SystematicName,l01CRC$SystematicName,l01CSC$SystematicName)) #combining the genes lists cutoff by 0.01 of three cells for cluster Analysis
indices_of_systematicNames<-match(all_systematicNames,rownames(EG_avg))    #Indices of Systematic Names of Average Expression List
data_for_clustering_l01<-EG_avg[indices_of_systematicNames,]               #Selection of rows on basis of Average Expression
########################################Biclustering#####################################################################
#Biclustering technique used here for clustering genes and conditions highly related in the sub problem data.
#########Algorithm for finding the optimal parameter values for plaid algorithm with 0.05 cutoff############################################
y<-0      
continue <- TRUE
while(y <100){       #intiate loop untill number of iteration
  value1 <- sample(1:500, 1, replace=F)  #Generate Random Numbers for the iteration layer 
  print(value1)   
  y<-y+1                          #Increase y on every iteration
  res1 <- biclust(data_for_clustering$E, method=BCPlaid(), iter.startup=1,iter.layer=value1,verbose=FALSE)  #Performs plaid model biclustering
  myrows <- attributes(res1)$RowxNumber; 
  mycols <- attributes(res1)$NumberxCol; 
  myclusters <- lapply(1:length(myrows[1,]), function(x) list(rownames(data_for_clustering$E[myrows[,x], ]), colnames(data_for_clustering$E[, mycols[x,]]))) #Extract the clusteres
  names(myclusters) <- paste("CL", 1:length(myclusters), sep="_") #Name each Cluster by number
  number_length<-FALSE
  number_length_2<-FALSE
  number_length_3<-FALSE
  
  
  print("########################################################")
  for (x in 1:length(myclusters)) {     #iterate the all clusters and save cluster number which is most compitable with QC marker Genes 
    Number_of_true<-QC_genes$V1 %in% substr(myclusters[[x]][[1]],1,9)
    if(all(Number_of_true)){
      number_length<-TRUE
      a<-x
    }
  }
  for (y in 1:length(myclusters)) {    #iterate the all clusters and save cluster number which is most compitable with CRC marker Genes 
    Number_of_true_2<-CRC_genes$V1 %in% substr(myclusters[[y]][[1]],1,9)
    if(all(Number_of_true_2)){
      number_length_2<-TRUE
      b<-y
    }
  }
  for (z in 1:length(myclusters)) {    #iterate the all clusters and save cluster number which is most compitable with CSC marker Genes 
    Number_of_true_3<-CSC_genes$V1 %in% substr(myclusters[[z]][[1]],1,9)
    if(all(Number_of_true_3)){
      number_length_3<-TRUE
      c<-z
    }
  }
  
  
  if((number_length)&(number_length_2)&(number_length_3)){         #if three of the cluster numbers are true then print cluster information
    print(number_length)
    print(number_length_2)
    print(number_length_3)
    print(value1)
    print(QC_genes$V1 %in% substr(myclusters[[a]][[1]],1,9))
    print(CRC_genes$V1 %in% substr(myclusters[[b]][[1]],1,9))
    print(CSC_genes$V1 %in% substr(myclusters[[c]][[1]],1,9))
  }
}
#Show how much most compatible cluster is compatible with marker genes and sheow the cluster number of most compatible cluster.
total_number <- number_length+number_length_2+number_length_3
total_length<-(length(QC_genes$V1)+length(CRC_genes$V1)+length(CSC_genes$V1))  
total=total_number/total_length 
print(total)
print(paste0("Cluster number: ", i))
########################0.01 Marker Genes#############################################
QC_genes$V1<-ernst01QC
CSC_genes$V1<-ernst01CSC                       
CRC_genes$V1<-ernst01CRC
###################################Plaid Biclustering for 0.05 cutoff##############################
res1 <- biclust(data_for_clustering$E, method=BCPlaid(),iter.startup=66,iter.layer=1,verbose=TRUE) #Apply the Plaid algorithm on the Data set for clustering
# Performs plaid model biclustering as described in Turner et al, 2003. This algorithm models data matrices to a sum of
# layers, the model is fitted to data through minimization of error.
summary(res1) # The summary() functions returns the sizes of all clusters found
show(res1)    #the show() function provides an overview of the method applied.
names(attributes(res1)) 
# Converts res1 object of class Biclust into a list and returns the names of its components (slots). The results for the
# row clustering are stored in the "RowxNumber" slot and the results for the column cluster
myrows <- attributes(res1)$RowxNumber;    
mycols <- attributes(res1)$NumberxCol; 
myclusters <- lapply(1:length(myrows[1,]), function(x) list(rownames(data_for_clustering$E[myrows[,x], ]), colnames(data_for_clustering$E[, mycols[x,]])))
names(myclusters) <- paste("CL", 1:length(myclusters), sep="_")
myclusters

write.table(file="clusters_0.05/Plaid/clust_plaid_01(QC_genes_validated_marker_genes).txt", substr(myclusters$CL_1[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.05/Plaid/clust_plaid_02(CRC_genes_validated_marker_genes).txt", substr(myclusters$CL_2[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.05/Plaid/clust_plaid_03.txt", substr(myclusters$CL_3[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.05/Plaid/clust_plaid_04(CSC_genes_validated_marker_genes).txt", substr(myclusters$CL_4[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.05/Plaid/clust_plaid_05.txt", substr(myclusters$CL_5[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.05/Plaid/clust_plaid_06.txt", substr(myclusters$CL_6[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.05/Plaid/clust_plaid_07.txt", substr(myclusters$CL_7[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.05/Plaid/clust_plaid_08.txt", substr(myclusters$CL_8[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.05/Plaid/clust_plaid_09.txt", substr(myclusters$CL_9[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.05/Plaid/clust_plaid_10.txt", substr(myclusters$CL_10[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)

############################Fabia_0.05######################################
number<-0
total_number<-0
y<-0
while(y <100){ 
  i <- sample(1:18,1, replace=F) #Generate Random Numbers for the total number of cluster
  value2<-sample(1:500, 1, replace=F) #Generate Random Numbers for the cyc
  number_length<-0
  number_length_2<-0
  number_length_3<-0
  print("########################################################")
  
  
  res_fabia <- fabia(X=data_for_clustering$E,p=i,alpha=0.01,cyc=value2,center=1) #Performs fabia model biclustering for data of 0.05 cutoff
  rb <- extractBic(res_fabia)  #extract clusters of fabia
  for (x in 1:i) {            #iterate the all clusters and save cluster number which is most compitable with QC marker Genes 
    Number_of_true<-length(which(QC_genes$V1 %in% substr(rb$bic[x,]$bixn,1,9)))
    if(Number_of_true>number_length){
      number_length<-Number_of_true
      a<-x
    }
  }
  for (y in 1:i) {      #iterate the all clusters and save cluster number which is most compitable with CRC marker Genes 
    Number_of_true_2<-length(which(CRC_genes$V1 %in% substr(rb$bic[y,]$bixn,1,9)))
    if(Number_of_true_2>number_length_2){
      number_length_2<-Number_of_true_2
      b<-y
    }
  }
  for (z in 1:i) {  #iterate the all clusters and save cluster number which is most compitable with CSC marker Genes 
    Number_of_true_3<-length(which(CSC_genes$V1 %in% substr(rb$bic[z,]$bixn,1,9)))
    if(Number_of_true_3>number_length_3){
      number_length_3<-Number_of_true_3
      c<-z
    }
  }
  total_length<-(length(QC_genes$V1)+length(CRC_genes$V1)+length(CSC_genes$V1))
  total_number <- number_length+number_length_2+number_length_3
  print(total_number)
  total=total_number/total_length    #Percentage of compatibilty with Marker Genes
  print(total)
  print(paste0("Cluster number: ", i))
  if(total==1&a!=b&b!=c&a!=c){       #if cluster is 100 percent compatible with Marker Genese then show the information of clusters
    if(!all(CRC_genes$V1 %in% substr(rb$bic[a,]$bixn,1,9))){    #if all CRC marker genes is not QC specific cluster, show information
      print(paste0("cluster QC: ",a))                 #print the cluster number of QC           
      print(paste0("Cluster CRC: ",b))                #print the cluster number of CRC
      print(paste0("cluster CSC: ",c))                #print the cluster number of CSC
      print(paste0("QC: ", (number_length/length(QC_genes$V1))))  #print how much QC specific cluster is compatible with QC Marker Genes
      print(paste0("CRC: ", (number_length_2/length(CRC_genes$V1))))   #print how much CRC specific cluster is compatible with CRC Marker Genes
      print(paste0("CSC: ", (number_length_3/length(CSC_genes$V1))))    #print how much CSC specific cluster is compatible with CSC Marker Genes
      print(paste0("Total: ", (total_number/total_length)))             #print Total compatiblity with Total Marker genes 
      print("###########################################################")
      QC_genes$V1 %in% substr(rb$bic[a,]$bixn,1,9)    #Query QC marker Genes in QC specific cluster
      QC_genes$V1 %in% substr(rb$bic[b,]$bixn,1,9)    #Query QC marker Genes in CRC specific cluster
      QC_genes$V1 %in% substr(rb$bic[c,]$bixn,1,9)    #Query QC marker Genes in CSC specific cluster
      print("############################################################")
      CRC_genes$V1 %in% substr(rb$bic[a,]$bixn,1,9)    #Query CRC marker Genes in QC specific cluster
      CRC_genes$V1 %in% substr(rb$bic[b,]$bixn,1,9)    #Query CRC marker Genes in CSC specific cluster
      CRC_genes$V1 %in% substr(rb$bic[c,]$bixn,1,9)    #Query CRC marker Genes in CRC specific cluster
      print("#############################################################")
      CSC_genes$V1 %in% substr(rb$bic[a,]$bixn,1,9)   #Query CSC marker Genes in QC specific cluster
      CSC_genes$V1 %in% substr(rb$bic[b,]$bixn,1,9)   #Query CSC marker Genes in CSC specific cluster
      CSC_genes$V1 %in% substr(rb$bic[c,]$bixn,1,9)    #Query CSC marker Genes in CRC specific cluster
      print("##########################################################")
      clusters_numbers<-c(a,b,c)                   #save clutster numbers of QC ,CSC and CRC specific genes in cluster numbers
      #Save QC,CSC,CRC specific cluster in the txt file
      write.table(substr(rb$bic[a,]$bixn,1,9),file = "QC_cluster_genes_cluster_1_15.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
      write.table(substr(rb$bic[b,]$bixn,1,9),file = "CRC_cluster_genes_cluster_2_15.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
      write.table(substr(rb$bic[c,]$bixn,1,9),file = "CSC_cluster_genes_cluster_3_15.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
      #save clusters number in txt file
      write.table(clusters_numbers,file = "cluster_number.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)   
      #Save all clusters
      for (s in 1:i ) {
        write.table(substr(rb$bic[s,]$bixn,1,9),file = paste0("cluster_genes_05_",s),quote = FALSE,row.names = FALSE,col.names = FALSE)
        
      }
    }
  }
}



########Algorithm for finding the optimal parameter values for plaid algorithm with 0.01 cutoff############################################
y<-0      
continue <- TRUE
while(y <100){       #intiate loop untill number of iteration
  value1 <- sample(1:500, 1, replace=F)  #Generate Random Numbers for the iteration layer 
  print(value1)   
  y<-y+1                          #Increase y on every iteration
  res1 <- biclust(data_for_clustering_l01$E, method=BCPlaid(), iter.startup=1,iter.layer=value1,verbose=FALSE)  #Performs plaid model biclustering
  myrows <- attributes(res1)$RowxNumber; 
  mycols <- attributes(res1)$NumberxCol; 
  myclusters <- lapply(1:length(myrows[1,]), function(x) list(rownames(data_for_clustering_l01$E[myrows[,x], ]), colnames(data_for_clustering_l01$E[, mycols[x,]]))) #Extract the clusteres
  names(myclusters) <- paste("CL", 1:length(myclusters), sep="_") #Name each Cluster by number
  number_length<-FALSE
  number_length_2<-FALSE
  number_length_3<-FALSE
  
  
  print("########################################################")
  for (x in 1:length(myclusters)) {     #iterate the all clusters and save cluster number which is most compitable with QC marker Genes 
    Number_of_true<-QC_genes$V1 %in% substr(myclusters[[x]][[1]],1,9)
    if(all(Number_of_true)){
      number_length<-TRUE
      a<-x
    }
  }
  for (y in 1:length(myclusters)) {    #iterate the all clusters and save cluster number which is most compitable with CRC marker Genes 
    Number_of_true_2<-CRC_genes$V1 %in% substr(myclusters[[y]][[1]],1,9)
    if(all(Number_of_true_2)){
      number_length_2<-TRUE
      b<-y
    }
  }
  for (z in 1:length(myclusters)) {    #iterate the all clusters and save cluster number which is most compitable with CSC marker Genes 
    Number_of_true_3<-CSC_genes$V1 %in% substr(myclusters[[z]][[1]],1,9)
    if(all(Number_of_true_3)){
      number_length_3<-TRUE
      c<-z
    }
  }
  
  
  if((number_length)&(number_length_2)&(number_length_3)){         #if three of the cluster numbers are true then print cluster information
    print(number_length)
    print(number_length_2)
    print(number_length_3)
    print(value1)
    print(QC_genes$V1 %in% substr(myclusters[[a]][[1]],1,9))
    print(CRC_genes$V1 %in% substr(myclusters[[b]][[1]],1,9))
    print(CSC_genes$V1 %in% substr(myclusters[[c]][[1]],1,9))
  }
}
#Show how much most compatible cluster is compatible with marker genes and sheow the cluster number of most compatible cluster.
total_number <- number_length+number_length_2+number_length_3
total_length<-(length(QC_genes$V1)+length(CRC_genes$V1)+length(CSC_genes$V1))  
total=total_number/total_length 
print(total)
print(paste0("Cluster number: ", i))
############################Fabia_0.01######################################
number<-0
total_number<-0
y<-0
while(y <100){  #intiate loop untill number of iteration
  i <- sample(1:18,1, replace=F)   #Generate Random Numbers for the total number of cluster
  value2<-sample(1:500, 1, replace=F) #Generate Random Numbers for the cyc
  number_length<-0
  number_length_2<-0
  number_length_3<-0
  print("########################################################")
  
  
  res_fabia <- fabia(X=data_for_clustering_l01$E,p=i,alpha=0.01,cyc=value2,center=1) #Performs fabia model biclustering for data of 0.01 cutoff
  rb <- extractBic(res_fabia) #extract clusters of fabia
  for (x in 1:i) {            #iterate the all clusters and save cluster number which is most compitable with QC marker Genes 
    Number_of_true<-length(which(QC_genes$V1 %in% substr(rb$bic[x,]$bixn,1,9)))
    if(Number_of_true>number_length){
      number_length<-Number_of_true
      a<-x
    }
  }
  for (y in 1:i) {      #iterate the all clusters and save cluster number which is most compitable with CRC marker Genes 
    Number_of_true_2<-length(which(CRC_genes$V1 %in% substr(rb$bic[y,]$bixn,1,9)))
    if(Number_of_true_2>number_length_2){
      number_length_2<-Number_of_true_2
      b<-y
    }
  }
  for (z in 1:i) {    #iterate the all clusters and save cluster number which is most compitable with CSC marker Genes 
    Number_of_true_3<-length(which(CSC_genes$V1 %in% substr(rb$bic[z,]$bixn,1,9)))
    if(Number_of_true_3>number_length_3){
      number_length_3<-Number_of_true_3
      c<-z
    }
  }#Show how much most compatible cluster is compatible with marker genes and sheow the cluster number of most compatible cluster.
  total_length<-(length(QC_genes$V1)+length(CRC_genes$V1)+length(CSC_genes$V1))  
  total_number <- number_length+number_length_2+number_length_3
  print(total_number)
  total=total_number/total_length           #Percentage of compatibilty with Marker Genes
  print(total)
  print(paste0("Cluster number: ", i))      #Print Cluster numbers
  if(total==1&a!=b&b!=c&a!=c){              #if cluster is 100 percent compatible with Marker Genese then show the information of clusters
    if(!all(CRC_genes$V1 %in% substr(rb$bic[a,]$bixn,1,9))){ #if all CRC marker genes is not QC specific cluster, show information
      print(paste0("cluster QC: ",a))       #print the cluster number of QC 
      print(paste0("Cluster CRC: ",b))      #print the cluster number of CSC
      print(paste0("cluster CSC: ",c))      #print the cluster number of CRC 
      print(paste0("QC: ", (number_length/length(QC_genes$V1)))) #print how much QC specific cluster is compatible with QC Marker Genes
      print(paste0("CRC: ", (number_length_2/length(CRC_genes$V1))))  #print how much CRC specific cluster is compatible with CRC Marker Genes
      print(paste0("CSC: ", (number_length_3/length(CSC_genes$V1))))  #print how much CSC specific cluster is compatible with CSC Marker Genes
      print(paste0("Total: ", (total_number/total_length)))           #print Total compatiblity with Total Marker genes 
      print("###########################################################")
      QC_genes$V1 %in% substr(rb$bic[a,]$bixn,1,9)                    #Query QC marker Genes in QC specific cluster
      QC_genes$V1 %in% substr(rb$bic[b,]$bixn,1,9)                    #Query QC marker Genes in CRC specific cluster
      QC_genes$V1 %in% substr(rb$bic[c,]$bixn,1,9)                    #Query QC marker Genes in CSC specific cluster
      print("############################################################")
      CRC_genes$V1 %in% substr(rb$bic[a,]$bixn,1,9)                   #Query CRC marker Genes in QC specific cluster
      CRC_genes$V1 %in% substr(rb$bic[b,]$bixn,1,9)                   #Query CRC marker Genes in CSC specific cluster
      CRC_genes$V1 %in% substr(rb$bic[c,]$bixn,1,9)                   #Query CRC marker Genes in CRC specific cluster
      print("#############################################################")
      CSC_genes$V1 %in% substr(rb$bic[a,]$bixn,1,9)                   #Query CSC marker Genes in QC specific cluster
      CSC_genes$V1 %in% substr(rb$bic[b,]$bixn,1,9)                   #Query CSC marker Genes in CSC specific cluster
      CSC_genes$V1 %in% substr(rb$bic[c,]$bixn,1,9)                   #Query CSC marker Genes in CRC specific cluster
      print("##########################################################")
      clusters_numbers<-c(a,b,c)                                       #save clutster numbers of QC ,CSC and CRC specific genes in cluster numbers
      #Save QC,CSC,CRC specific cluster in the txt file
      write.table(substr(rb$bic[a,]$bixn,1,9),file = "QC_cluster_genes_cluster_1_11.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
      write.table(substr(rb$bic[b,]$bixn,1,9),file = "CRC_cluster_genes_cluster_2_11.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
      write.table(substr(rb$bic[c,]$bixn,1,9),file = "CSC_cluster_genes_cluster_3_11.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
      #save clusters number in txt file
      write.table(clusters_numbers,file = "cluster_number.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
      #Save all clusters
      for (s in 1:i ) {
        write.table(substr(rb$bic[s,]$bixn,1,9),file = paste0("QC_cluster_genes_3_11_",s),quote = FALSE,row.names = FALSE,col.names = FALSE)
        
      }
    }
  }
}

###################################Plaid Biclustering for 0.01 cutoff##############################
res2 <- biclust(data_for_clustering_l01$E, method=BCPlaid(),iter.startup=66,iter.layer=1,verbose=TRUE)
summary(res2); 
show(res2) 
names(attributes(res2)) 
myrows_l01 <- attributes(res2)$RowxNumber; 
mycols_l01 <- attributes(res2)$NumberxCol; 
myclusters_l01 <- lapply(1:length(myrows_l01[1,]), function(x) list(rownames(data_for_clustering_l01$E[myrows_l01[,x], ]), colnames(data_for_clustering_l01$E[, mycols_l01[x,]])))
names(myclusters_l01) <- paste("CL", 1:length(myclusters), sep="_")

write.table(file="clusters_0.01/Plaid/clust_plaid_01(CRC_genes_validated_marker_genes).txt", substr(myclusters_l01$CL_1[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.01/Plaid/clust_plaid_02(QC_genes_validated_marker_genes).txt", substr(myclusters_l01$CL_2[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.01/Plaid/clust_plaid_03(CSC_genes_validated_marker_genes).txt", substr(myclusters_l01$CL_3[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.01/Plaid/clust_plaid_04.txt", substr(myclusters_l01$CL_4[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.01/Plaid/clust_plaid_05.txt", substr(myclusters_l01$CL_5[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.01/Plaid/clust_plaid_06.txt", substr(myclusters_l01$CL_6[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.01/Plaid/clust_plaid_07.txt", substr(myclusters_l01$CL_7[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.01/Plaid/clust_plaid_08.txt", substr(myclusters_l01$CL_8[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.01/Plaid/clust_plaid_09.txt", substr(myclusters_l01$CL_9[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)
write.table(file="clusters_0.01/Plaid/clust_plaid_10.txt", substr(myclusters_l01$CL_10[[1]],1,9),quote=FALSE,row.names=FALSE,col.names = FALSE)

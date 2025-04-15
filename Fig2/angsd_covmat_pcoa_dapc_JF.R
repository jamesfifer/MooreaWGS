getwd()

setwd("C:/Users/james/OneDrive - UC San Diego/Documents/BOSTON/Davies/Strader/code/PCA")
#read in metadata
library(stringr)
library(tidyverse)
#library(adegenet)
library(vegan)
library(dplyr)
library(ggplot2)

bams=data.frame(read.table("../bamscl", header = FALSE)) # list of bam files
#For MayJuv PCA
bams=data.frame(read.table("MayJuvbamscl", header = FALSE)) # list of bam files


#for May only all sites
bams=data.frame(read.table("Maybamscl", header = FALSE)) # list of bam files
#for dapc same sites MayOct Juv
bams=data.frame(read.table("MayOctJuvbamscl", header = FALSE)) # list of bam files
#redoing pca without clone/outliers same sites MayOct Juv
bams=data.frame(read.table("redoMayOctJuvbamscl", header = FALSE)) # list of bam files



#for dapc Juv Assignment
bams=data.frame(read.table("MayOctJuv_2bamscl", header = FALSE)) # list of bam files
#for dapc MayOct overlap
bams=data.frame(read.table("MayOctbamscl4vcf",header=F))
bams=ifelse(length(grep("Oct",bams$V1))>0,data.frame(lapply(bams, function(x) {
  gsub("Oct.", "", x)
})),bams)
bams=as.data.frame(do.call(cbind, bams))
#for dapc MayOct2FR and MayOctFR outlier 
bams=data.frame(read.table("MayOctJuvbamscl", header = FALSE)) # list of bam files


colnames(bams)<- "Adapter"
#removing directory string
bams <- data.frame(lapply(bams, function(x) {
  gsub(".*/", "", x)
}))

#extracting sample name from .bam file
samplenames=word(bams$Adapter,1,sep = "\\_S")
#removing _'s from some of the samples
samplenames <- lapply(samplenames, function(x) {
  gsub("_", "", x)
})
bams$Sample=unlist(samplenames)
#rename two samples to match metadata
bams$Sample=str_replace(bams$Sample, "M1B1$", "M1B01")
bams$Sample=str_replace(bams$Sample, "N3F18$", "N3F18J")

#read in meta data
samples=readRDS(file="../forsubbamscl.RData")

#Making supplemental data file
# what=lapply(samples, function(x) {
#   as.vector(unlist(as.character(x)))
# })
# combined_df <- as.data.frame(do.call(cbind, what))
# combined_df2=combined_df %>%
#   select(Sample, CloneGrp,Adapter,Site, Habitat, Month_sampled,Depth_m,BleachingScoreMay,area)
# write.csv(combined_df2, file="../Allsamplesmeta.csv")


#merge bams and metadata
new<-merge(bams,samples, by ="Adapter")
new = new[!duplicated(new$Sample.y),] #some samples were duplicated unnecessarily for some reason
#remove the sample that was backreef in May and forereef in Oct (removed Oct only)
new=subset(new, Sample.y!="O5F35" & Sample.y!="N2F32J" & Sample.x!="N2F32J"
             )
  
#new = new[order(new[,'Adapter.x']),]
#new = new[!duplicated(new$Adapter),] #some samples were duplicated unnecessarily for some reason
#new<-new[ order(match(new$Adapter.x, bams$Adapter)), ] #MAKE SURE ORDER FOR THIS FILE MATCHES ORDER FOR BAM FILE USED FOR IBS CREATION OTHERWISE DENDO WILL BE WRONG
new$Site=str_replace(new$Site, "/Pile", "")
new$Habitat=str_replace(new$Habitat, "Unknown", "forereef")
new$Site=str_replace(new$Site, "Unknown", "LTER5")

new$Month_sampled[is.na(new$Month_sampled) ] <- "Juv"
new$Site.Habitat.Month=paste0(new$Site,new$Habitat,new$Month_sampled)
new$Site.Habitat=paste0(new$Site,new$Habitat)



#--------------------
#run above code 'loading invidividual.." first
# covariance / PCA 

library(vegan)
# choose either of the following covarince matrices:

#All loci all sampes
co = as.matrix(read.table("bayescan.covMat")) # covariance based on single-read resampling

#May all loci 
co = as.matrix(read.table("MayGLM.GWAS.covMat")) # covariance based on single-read resampling
#for dapc baypass loci
co = as.matrix(read.table("MayOctBaypassMayOctJuv.covMat")) # covariance based on single-read resampling
#for dapc all loci (Juv May Oct overlap)
co = as.matrix(read.table("MayOctJuvAllLoci.covMat")) # covariance based on single-read resampling
#for dapc all loci (Juv May Oct overlap redos)
co = as.matrix(read.table("redoMayOctJuvAllLoci.covMat")) # covariance based on single-read resampling


#for dapc baypass loci (May Oct overlap)
co = as.matrix(read.table("MayOctBaypassMayOctJuv_2.covMat")) # covariance based on single-read resampling
#for dapc all loci (May Oct overlap) for juv assignment
co = as.matrix(read.table("MayOctJuvAllLoci_2.covMat")) # covariance based on single-read resampling
#for dapc/pca all loci May Oct only
co = as.matrix(read.table("MayOctAllLoci.covMat")) # covariance based on single-read resampling

#for dapc/pca gwas intermediate check
co = as.matrix(read.table("MayOctJuvGWAS.MayOct.plink.all.PrePost.glm.logistic.0.01.sites.covMat")) # covariance based on single-read resampling
co = as.matrix(read.table("MayOctJuvGWAS.MayOct.plink.all.PrePost.glm.logistic.0.001.sites.covMat")) # covariance based on single-read resampling
co = as.matrix(read.table("MayOctJuvGWAS.MayOct.plink.all.PrePost.glm.logistic.0.0001.sites.covMat")) # covariance based on single-read resampling

#for pca MayJuv Gwas covariate
co = as.matrix(read.table("MayJuvAllLoci.covMat")) # covariance based on single-read resampling

#for dapc/pca gwas intermediate check FR
co = as.matrix(read.table("MayOctJuvGWAS.MayOct2FR.plink.all.PrePost.glm.logistic.0.1.sites.covMat")) # covariance based on single-read resampling




#co = as.matrix(read.table("ok.covar")) # covariance by ngsCovar
dimnames(co)=list(bams$Adapter,bams$Adapter)
#removing the NA

 co=co[!rownames(co) %in% "Lane3_10HB_O_S84_L003_minq20_sorted.bam", !colnames(co) %in% "Lane3_10HB_O_S84_L003_minq20_sorted.bam"]
 #removing the weird juvenile N2F32_J which is a clone at LTER3forereef and LTER2
 co=co[!rownames(co) %in% "N2F32J_S346.dedup.overlapclipped.bam", !colnames(co) %in% "N2F32J_S346.dedup.overlapclipped.bam"]
 #Only remvoe this for outlier analysis, was a weirdo/clone 
 co=co[!rownames(co) %in% "12N_O_S37.dedup.overlapclipped.bam", !colnames(co) %in% "12N_O_S37.dedup.overlapclipped.bam"]
 
 #^Just re-run the analysis without these weirdos
 #write.table(row.names(co),file="redoMayOctJuvbamscl", quote=F,row.names = F, col.names = F)

 
 
 
 # #removing all reps
# #create these using the Clonecheck.R script
# keep=noClones$Adapter
# newco=co[rownames(co) %in% keep, colnames(co) %in% keep]

 


# PCoA and CAP (constranied analysis of proximities)  
#pp0=capscale(as.dist(1-cov2cor(newco))~1) # PCoA


pp0=capscale(as.dist(1-cov2cor(co))~1) # PCoA
#pp=capscale(as.dist(1-cov2cor(co))~site,conds) # CAP

# significance of by-site divergence, based on 1-correlation as distance
#adonis(as.dist(1-cov2cor(co))~site,conds)

# eigenvectors: how many are interesting?
# Assuming you have already performed PCA and stored the eigenvalues in the variable "eigenvalues"
# eigenvalues should be a vector of eigenvalues in descending order

# Create a scree plot
library(ggplot2)


plot(pp0$CA$eig) 

kmeans.axes <- sum(diff(summary(eigenvals(pp0))[2,])*-1 > 0.01) #based on eigen plot
kmeans.axes
#0
# This function finds the total number of eigenvectors in which the difference between last eigenvalue and the one preceding is greater than 0.01
# Standardized method of selecting eigenvectors that contribute the most explanatory power in clustering

axes2plot=c(1,2)  
#quartz()
cc=pp0

cmd=pp0
cmd$CA$u

library(dplyr)
library(ggplot2)

pca_s <- as.data.frame(cmd$CA$u[,axes2plot])
#

row.names(pca_s)=str_replace(row.names(pca_s),"X","")
pca_s$Adapter=row.names(pca_s)
##Can go to DAPC section from here 

#saveRDS(pca_s, file="MayJuvPCA")

new2<-merge(pca_s,new, by ="Adapter")
#new2=new2[!duplicated(new2$Coral.ID),]

#levels(pca_s$site) <- list(NC_TriangleWrecks  = "NC", NC_RadioIsland = "NC_Radio",Virginia="VA",Massachusetts="WH",Rhode_Island="RI",
#                           TX_PortA="Texas_PortA",TX_Packery="Texas_Packery", NC_RadioIsland_Oculina="Oculina_NC",FL_Panacea="Panacea")

colors=c("#2b8cbe","#d40000ff","#ff8080ff")
#pca_s$site=factor(pca_s$site,levels=c("NC_RadioIsland","NC_TriangleWrecks","Virginia","Rhode_Island","Massachusetts","NC_RadioIsland_Oculina","TX_Packery","TX_PortA","FL_Panacea"))
names(colors)=c("October","May","Juv")
eigenvals=pp0$CA$eig
#calc variance explained
varexplained=eigenvals/sum(eigenvals)


e <- eigen(co)
eigenvals=e$values
  varexplained=eigenvals/sum(eigenvals)

plot(e$vectors[,1:2],lwd=2,ylab="PC 2",xlab="PC 2",main="Principal components",col=rep(1:3,each=10),pch=16)
varexplained=eigenvals/sum(eigenvals)


#All
#ggplot(pca_s, aes(MDS1, MDS2, label=pca_s$label)) +
ggplot(new2, aes(MDS1, MDS2)) +
  geom_point(data=new2, size=3.5,aes(colour=as.factor(Site.Habitat), shape=as.factor(Month_sampled)))  +
  theme_classic(base_size = 20) +
  #geom_text(hjust=0, vjust=0)+
  # geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
  stat_ellipse(geom="polygon",level=.75, alpha = 1/3, aes(fill = as.factor(new2$Month_sampled)))+
  #scale_fill_manual(values=c('#f0f9e8','#bae4bc','#7bccc4','#43a2ca','#0868ac'))+
  xlab(paste0("MDS1 (",formatC((varexplained[1]*100), digits = 2, format = "f"),"%)")) +
  ylab(paste0("MDS2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
  #scale_shape_manual(values = rep(10:19),name="")+
  #scale_fill_manual(values=colors, name="", labels = c("October", "May", "Juv", guide=FALSE)) #
  #scale_colour_manual(values=colors, name="", guide=FALSE) #
  #scale_fill_manual(values=colors, name="") #
xlim(c(-0.1,0.6))+ylim(c(-0.4,0.2))+guides(fill=FALSE) 
#ADONIS 
pca_table=cbind(pca_s$Adapter,cmd$CA$u)
colnames(pca_table)[1]="Adapter"
new2<-merge(pca_table,new, by ="Adapter")
#new2=new2[!duplicated(new2$Coral.ID),]
new2<-new2[ order(match(new2$Adapter, bams$Adapter)), ] #MAKE SURE ORDER FOR THIS FILE MATCHES ORDER FOR BAM FILE USED FOR IBS CREATION OTHERWISE DENDO WILL BE WRONG
conds=new2$Site.Habitat

AD=adonis2(as.dist(1-cov2cor(co))~as.factor(Site.Habitat), new2)
summary(AD)
p_value <- AD$aov.tab$"Pr(>F)"["new2$Site.Habitat"]
p_value
AD=adonis2(as.dist(1-cov2cor(co))~as.factor(Month_sampled), new2)
summary(AD)
p_value <- AD$aov.tab$"Pr(>F)"["new2$Site.Habitat"]
p_value


#May Only
new2$Month_sampled=as.factor(new2$Month_sampled)
new2=dplyr::filter(new2, Site.Habitat !="Unknown")
new2 <- dplyr::filter(new2, !(Month_sampled == "October"))

new2$Month_sampled <- droplevels(new2$Month_sampled)
levels(new2$Month_sampled)
levels(as.factor(new2$Site.Habitat))

shapes=rep(c(13,3:4,6,21:23,11,24))
shapes
oranges=rep("orange",9)
shapes
shape_names <- c("LTER1backreef", "LTER1crest", "LTER2forereef", "LTER2crest", "LTER3crest",
                 "LTER1forereef", "LTER3forereef", "LTER5backreef", "LTER5forereef")
names(shapes)=shape_names
shapes
levels(as.factor(new2$Site.Habitat))
Mayplot=ggplot(new2, aes(MDS1, MDS2)) +
  geom_point(data=new2, size=3.5,aes(fill=as.factor(Site.Habitat), shape=as.factor(Site.Habitat)))  +
  theme_classic(base_size = 26) +
  scale_fill_manual(values = oranges,labels=names(shapes),name="Group",drop=TRUE,limits = force) +
  #stat_ellipse(geom="polygon",level=.75, alpha = 1/3, aes(fill = as.factor(new2$Month_sampled)))+
  #scale_fill_manual(values=c('#f0f9e8','#bae4bc','#7bccc4','#43a2ca','#0868ac'))+
  xlab(paste0("PCo1 (",formatC((varexplained[1]*100), digits = 2, format = "f"),"%)")) +
  ylab(paste0("PCo2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
  scale_shape_manual(values=shapes, labels=names(shapes),name="Group") #+
 # xlim(c(-0.1,0.6))+ylim(c(-0.4,0.2))+guides(color=FALSE)+theme(aspect.ratio = 1) #+guides(color=FALSE) 
ggplotly(Mayplot)
Mayplot
saveRDS(Mayplot, file="Mayplot.RDS")
save.image(file="Mayplot.RData")
#ADONIS 
pca_table=cbind(pca_s$Adapter,cmd$CA$u)
colnames(pca_table)[1]="Adapter"
new2<-subset(merge(pca_table,new, by ="Adapter"),Month_sampled=="May")
new2<-new2[ order(match(new2$Adapter, bams$Adapter)), ] #MAKE SURE ORDER FOR THIS FILE MATCHES ORDER FOR BAM FILE USED FOR IBS CREATION OTHERWISE DENDO WILL BE WRONG
conds=new2$Site.Habitat

AD=adonis2(as.dist(1-cov2cor(co))~as.factor(Site.Habitat), new2)
summary(AD)
p_value <- AD$aov.tab$"Pr(>F)"["new2$Site.Habitat"]
p_value
#saveRDS(new2[1:3], file="MayPCA")


#Overlapping sites between May Oct only using MayOct
new2=subset(new2, Site.Habitat !="LTER5backreef")
ggplot(new2, aes(MDS1, MDS2)) +
  geom_point(data=new2, size=3.5,aes(fill=as.factor(Month_sampled),color=as.factor(Month_sampled), shape=as.factor(Site.Habitat)))  +
  theme_classic(base_size = 20) +
  #geom_text(hjust=0, vjust=0)+
  # geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
  scale_fill_manual(values = c("orange", "darkblue"),name="") +
  scale_color_manual(values = c("orange", "darkblue"),name="") +
  stat_ellipse(geom="polygon",level=.75, alpha = 1/3, aes(fill = as.factor(new2$Month_sampled)))+
  #scale_fill_manual(values=c('#f0f9e8','#bae4bc','#7bccc4','#43a2ca','#0868ac'))+
  xlab(paste0("MDS1 (",formatC((varexplained[1]*100), digits = 2, format = "f"),"%)")) +
  ylab(paste0("MDS2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
  scale_shape_manual(values = rep(c(18,21:26)),name="")+
    #scale_fill_manual(values=colors, name="", labels = c("October", "May", "Juv", guide=FALSE)) #
  #scale_colour_manual(values=colors, name="", guide=FALSE) #
  #scale_fill_manual(values=colors, name="") #
  xlim(c(-0.1,0.6))+ylim(c(-0.4,0.2))+guides(fill=FALSE) +guides(color=FALSE) 
saveRDS(new2[1:3], file="MayOctPCA")



#Overlapping sites only all 3

new2=subset(new2, Site.Habitat !="LTER5backreef")


p2=ggplot(new2, aes(MDS1, MDS2)) +
  geom_point(data=new2, size=3.5,aes(fill=as.factor(Month_sampled), shape=as.factor(Site.Habitat)))  +
  theme_classic(base_size = 20) +
  #geom_text(hjust=0, vjust=0, label=new2$Adapter)+
  # geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
  scale_fill_manual(values = c("pink", "orange", "darkblue"),name="") +
  #stat_ellipse(geom="polygon",level=.75, alpha = 1/3, aes(fill = as.factor(new2$Month_sampled)))+
  xlab(paste0("PCo1 (",formatC((varexplained[1]*100), digits = 2, format = "f"),"%)")) +
  ylab(paste0("PCo2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
  scale_shape_manual(values = rep(21:26),name="")+
 # scale_fill_manual(values=colors, name="", labels = c("October", "May", "Juv", guide=FALSE)) #
  #scale_colour_manual(values=colors, name="", guide=FALSE) #
  #scale_fill_manual(values=colors, name="") #
  xlim(c(-0.25,0.65))+ylim(c(-0.25,0.65))+
  #guides(colour=FALSE)+
guides(fill = guide_legend(override.aes = list(alpha = 1))) 
ggplotly(p2)
#N1F8J/N2F36J and N5F29J/N5F26J possible sibling pairs 
p2

saveRDS(p2, file="JuvOctMayPCA.RDS")  
save.image(file="JuvOctMayPCA.RData")
##Run for different months
for(i in seq_along(levels(as.factor(new$Month_sampled)))){
  i=2
  co = as.matrix(read.table("bayescan.covMat")) # covariance based on single-read resampling
  co=co[!rownames(co) %in% "Lane3_10HB_O_S84_L003_minq20_sorted.bam", !colnames(co) %in% "Lane3_10HB_O_S84_L003_minq20_sorted.bam"]
  #removing the weird juvenile N2F32_J which is a clone at LTER3forereef and LTER2
  co=co[!rownames(co) %in% "N2F32J_S346.dedup.overlapclipped.bam", !colnames(co) %in% "N2F32J_S346.dedup.overlapclipped.bam"]
  #co = as.matrix(read.table("ok.covar")) # covariance by ngsCovar
  dimnames(co)=list(bams$Adapter,bams$Adapter)
  SUB=levels(as.factor(new$Month_sampled))[i]
  SUBDF=subset(new, Month_sampled==SUB)
  keep=SUBDF$Adapter
  newco=co[rownames(co) %in% keep, colnames(co) %in% keep]
  pp0=capscale(as.dist(1-cov2cor(newco))~1) # PCoA
  # eigenvectors: how many are interesting?
  #plot(pp0$CA$eig) 
  kmeans.axes <- sum(diff(summary(eigenvals(pp0))[2,])*-1 > 0.01) #based on eigen plot
  kmeans.axes
  #2
  # This function finds the total number of eigenvectors in which the difference between last eigenvalue and the one preceding is greater than 0.01
  # Standardized method of selecting eigenvectors that contribute the most explanatory power in clustering
  axes2plot=c(1,2)  
  #quartz()
  cc=pp0
  cmd=pp0
  cmd$CA$u
  library(dplyr)
  library(ggplot2)
  pca_s <- as.data.frame(cmd$CA$u[,axes2plot])
  #
  pca_s$Adapter=row.names(pca_s)
  new2<-result <- merge(pca_s, new, by = "Adapter") %>%
    filter(Month_sampled == SUB)
  saveRDS(new2[1:3], file=paste0(SUB,"PCA"))
  eigenvals=pp0$CA$eig
  #calc variance explained
  varexplained=eigenvals/sum(eigenvals)
  
  p=ggplot(new2, aes(MDS1, MDS2)) +
    geom_point(data=new2, size=3.5,aes(colour=as.factor(Month_sampled), shape=as.factor(Site.Habitat)))  +
    theme_classic(base_size = 20) +
    #geom_text(hjust=0, vjust=0)+
    # geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
    #stat_ellipse(geom="polygon",level=.75, alpha = 1/3, aes(fill = as.factor(pca_s$site)),show.legend = NULL)+
    #scale_fill_manual(values=c('#f0f9e8','#bae4bc','#7bccc4','#43a2ca','#0868ac'))+
    xlab(paste0("MDS1 (",formatC((varexplained[1]*100), digits = 2, format = "f"),"%)")) +
    ylab(paste0("MDS2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
    #scale_shape_manual(values = rep(10:19),name="")+
    #scale_fill_manual(values=colors, name="", labels = c("October", "May", "Juv", guide=FALSE)) #
    #scale_colour_manual(values=colors, name="", guide=FALSE) #
    #scale_fill_manual(values=colors, name="") #
    xlim(c(-0.1,0.6))+ylim(c(-0.4,0.2))+guides(fill=FALSE) +
    ggtitle(SUB)
  print(p)
}

  



#May only looking at the sites sampled accross three timepoints only
##Run for different combos of site.habitat
for(i in seq_along(levels(as.factor(new$Site.Habitat)))){
co = as.matrix(read.table("bayescan.covMat")) # covariance based on single-read resampling

#co = as.matrix(read.table("ok.covar")) # covariance by ngsCovar
dimnames(co)=list(bams$Adapter,bams$Adapter)
co=co[!rownames(co) %in% "Lane3_10HB_O_S84_L003_minq20_sorted.bam", !colnames(co) %in% "Lane3_10HB_O_S84_L003_minq20_sorted.bam"]
#removing the weird juvenile N2F32_J which is a clone at LTER3forereef and LTER2
co=co[!rownames(co) %in% "N2F32J_S346.dedup.overlapclipped.bam", !colnames(co) %in% "N2F32J_S346.dedup.overlapclipped.bam"]


SUB=levels(as.factor(new$Site.Habitat))[i]
SUBDF=subset(new, Site.Habitat==SUB)
keep=SUBDF$Adapter
newco=co[rownames(co) %in% keep, colnames(co) %in% keep]
pp0=capscale(as.dist(1-cov2cor(newco))~1) # PCoA
# eigenvectors: how many are interesting?
#plot(pp0$CA$eig) 
kmeans.axes <- sum(diff(summary(eigenvals(pp0))[2,])*-1 > 0.01) #based on eigen plot
kmeans.axes
#2
# This function finds the total number of eigenvectors in which the difference between last eigenvalue and the one preceding is greater than 0.01
# Standardized method of selecting eigenvectors that contribute the most explanatory power in clustering
axes2plot=c(1,2)  
#quartz()
cc=pp0
cmd=pp0
cmd$CA$u
library(dplyr)
library(ggplot2)
pca_s <- as.data.frame(cmd$CA$u[,axes2plot])
#
pca_s$Adapter=row.names(pca_s)
new2<-merge(pca_s,new, by ="Adapter")
eigenvals=pp0$CA$eig
#calc variance explained
varexplained=eigenvals/sum(eigenvals)

p=ggplot(new2, aes(MDS1, MDS2)) +
  geom_point(data=new2, size=3.5,aes(colour=as.factor(Habitat), shape=as.factor(Month_sampled)))  +
  theme_classic(base_size = 20) +
  #geom_text(hjust=0, vjust=0)+
  # geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
  #stat_ellipse(geom="polygon",level=.75, alpha = 1/3, aes(fill = as.factor(pca_s$site)),show.legend = NULL)+
  #scale_fill_manual(values=c('#f0f9e8','#bae4bc','#7bccc4','#43a2ca','#0868ac'))+
  xlab(paste0("MDS1 (",formatC((varexplained[1]*100), digits = 2, format = "f"),"%)")) +
  ylab(paste0("MDS2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
  #scale_shape_manual(values = rep(10:19),name="")+
  #scale_fill_manual(values=colors, name="", labels = c("October", "May", "Juv", guide=FALSE)) #
  #scale_colour_manual(values=colors, name="", guide=FALSE) #
  #scale_fill_manual(values=colors, name="") #
  xlim(c(-0.1,0.6))+ylim(c(-0.4,0.2))+guides(fill=FALSE) +
  ggtitle(SUB)
print(p)
}

##Run for different combos of LTERS
for(i in seq_along(levels(as.factor(new$Site)))){
  co = as.matrix(read.table("bayescan.covMat")) # covariance based on single-read resampling
  #co = as.matrix(read.table("ok.covar")) # covariance by ngsCovar
  dimnames(co)=list(bams$Adapter,bams$Adapter)
  SUB=levels(as.factor(new$Site))[i]
  SUBDF=subset(new, Site==SUB)
  keep=SUBDF$Adapter
  newco=co[rownames(co) %in% keep, colnames(co) %in% keep]
  pp0=capscale(as.dist(1-cov2cor(newco))~1) # PCoA
  # eigenvectors: how many are interesting?
  #plot(pp0$CA$eig) 
  kmeans.axes <- sum(diff(summary(eigenvals(pp0))[2,])*-1 > 0.01) #based on eigen plot
  kmeans.axes
  #2
  # This function finds the total number of eigenvectors in which the difference between last eigenvalue and the one preceding is greater than 0.01
  # Standardized method of selecting eigenvectors that contribute the most explanatory power in clustering
  axes2plot=c(1,2)  
  #quartz()
  cc=pp0
  cmd=pp0
  cmd$CA$u
  library(dplyr)
  library(ggplot2)
  pca_s <- as.data.frame(cmd$CA$u[,axes2plot])
  #
  pca_s$Adapter=row.names(pca_s)
  new2<-merge(pca_s,new, by ="Adapter")
  eigenvals=pp0$CA$eig
  #calc variance explained
  varexplained=eigenvals/sum(eigenvals)
  
  p=ggplot(new2, aes(MDS1, MDS2)) +
    geom_point(data=new2, size=3.5,aes(colour=as.factor(Habitat), shape=as.factor(Month_sampled)))  +
    theme_classic(base_size = 20) +
    #geom_text(hjust=0, vjust=0)+
    # geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
    #stat_ellipse(geom="polygon",level=.75, alpha = 1/3, aes(fill = as.factor(pca_s$site)),show.legend = NULL)+
    #scale_fill_manual(values=c('#f0f9e8','#bae4bc','#7bccc4','#43a2ca','#0868ac'))+
    xlab(paste0("MDS1 (",formatC((varexplained[1]*100), digits = 2, format = "f"),"%)")) +
    ylab(paste0("MDS2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
    #scale_shape_manual(values = rep(10:19),name="")+
    #scale_fill_manual(values=colors, name="", labels = c("October", "May", "Juv", guide=FALSE)) #
    #scale_colour_manual(values=colors, name="", guide=FALSE) #
    #scale_fill_manual(values=colors, name="") #
    xlim(c(-0.1,0.6))+ylim(c(-0.4,0.2))+guides(fill=FALSE) +
    ggtitle(SUB)
  print(p)
}

#####DAPC

#
## PCA #############
#Can skip this step if using the PCA output created above (this just redoes the PCA)

  ## This function takes a covariance matrix and performs PCA. 
  # cov_matrix: a square covariance matrix generated by most pca softwares
  # ind_label: a vector in the same order and length as cov_matrix; it contains the individual labels of the individuals represented in the covariance matrix
  # pop_label: a vector in the same order and length as cov_matrix; it contains the population labels of the individuals represented in the covariance matrix
  # x_axis: an integer that determines which principal component to plot on the x axis
  # y_axis: an integer that determines which principal component to plot on the y axis
  # show.point: whether to show individual points
  # show.label: whether to show population labels
  # show.ellipse: whether to show population-specific ellipses
  # show.line: whether to show lines connecting population means with each individual point
  # alpha: the transparency of ellipses
  # index_exclude: the indices of individuals to exclude from the analysis
  # new2<-merge(pca_s,new, by ="Adapter")
  # new2_unique <- distinct(new2, Adapter, .keep_all = TRUE)
  # 
  # site_label=new2_unique$Site.Habitat
  # ind_label=new2_unique$Adapter
  # pop_label=new2_unique$Month_sampled
  # index_include <- setdiff(seq_along(ind_label), index_exclude)
  # m <- as.matrix(read.table("bayescan.covMat"))
  # m[is.na(m)]<- median(m, na.rm = T)
  # m<-m[index_include, index_include] ## Remove 4SJH, four 3Ps individuals, and contaminated ones
  # e <- eigen(m)
  # e_value<-e$values
  # x_axis=1
  # y_axis=2
  # x_variance<-e_value[x_axis]/sum(e_value)*100
  # y_variance<-e_value[y_axis]/sum(e_value)*100
  # e <- as.data.frame(e$vectors)
  # e <- cbind(ind_label[index_include], pop_label[index_include],site_label[index_include], e) ## with the above individuals removed
  # #colnames(e)[3:331]<-paste0("PC",1:329)
  # colnames(e)[3:(dim(e)[1])]<-paste0("PC",1:(dim(e)[1]-2)) ## with the above individuals removed
  # colnames(e)[1:3]<-c("individual", "population","site")
  # assign("pca_table", e, .GlobalEnv)
  # 

## DAPC ########

library(MASS)
library(cowplot)
library(dplyr)
library(vegan)


  ## This function should follow immediately after a PCA function. It takes the PCA output and performs a linear discriminant analysis
  # n: the number of principal components to use for discriminant analysis
  # x_axis: an integer that determines which linear discriminant to plot on the x axis
  # y_axis: an integer that determines which linear discriminant to plot on the y axis
  # show.point: whether to show individual points
  # show.label: whether to show population labels
  # show.ellipse: whether to show population-specific ellipses
  # show.line: whether to show lines connecting population means with each individual point
  # alpha: the transparency of ellipses
 
 #The DAPC evaluates shifts at all haplotypes at a locus, making this analysis potentially more sensitive to detecting selection than FST-based outlier analyses, which do not permit haplotype-specific comparisons. Pre- and post-SSWD groups were defined a priori. After running the a-score spline interpolation, which identified the optimal number of principal components to retain to be 124, we retained 89 to follow the <N/3 rule to avoid overfitting (68). The proportion of conserved variance was 47.7%, and one discriminant function (n groups b 1) was used.

  pca_table=cbind(pca_s$Adapter,cmd$CA$u)
  colnames(pca_table)[1]="Adapter"
  new2<-merge(pca_table,new, by ="Adapter") #this actually added a sample (M5B9, but that is correct dont know why it was missing from bam files for either May or Oct)

  pca_table1=cbind(new2$Adapter,new2$Site.Habitat,new2$Month_sampled, new2[c(2:(ncol(pca_table)))])
  colnames(pca_table1)[1:3]<-c("individual", "site","population")
pca_table=pca_table1
  
  # Convert columns 4 to ncol(pca_table) to numeric
  pca_table <- pca_table %>%
    mutate_at(vars(4:ncol(pca_table)), as.numeric)
  #will need to calculate A-score before doing n=50
  
  # Data frame for scree plot
  scree_data <- data.frame(Component = 1:length(pp0$CA$eig), Eigenvalue = pp0$CA$eig)
  # Plot the scree plot
  ggplot(scree_data, aes(x = Component, y = Eigenvalue)) +
    geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
    labs(x = "Component", y = "Eigenvalue", title = "Scree Plot") +
    theme_minimal()  
  #no obvious elbow for baypass so will just do number of components/3
  n=length(pp0$CA$eig)/3
  #obvious elbow at ~23 for all loci so will do 23 for all loci (also tried 15)
  n=23
#assumes that population is in column 3
  fit <- lda(population ~ ., data=pca_table[,3:(n+1)], na.action="na.omit", CV=F, output = "Scatterplot")

  plda <- predict(object = fit,
                  newdata = pca_table[,3:(n+1)])
  prop.lda = fit$svd^2/sum(fit$svd^2)
  dataset <- data.frame(group = pca_table[, 3], lda = plda$x)
  dataset$site=pca_table$site
  
#DAPC1 
# Calculate Euclidean distances to group centroids

  #
    newpca_table=subset(pca_table, population !="Juv")
    juvpca_table=subset(pca_table, population =="Juv")
    ?lda()
    fit <- lda(site ~ ., data=newpca_table[,c(2,4:(n+4))], na.action="na.omit", CV=F, output = "Scatterplot")
    # Empty vector to store cross-validation accuracy results
    accuracy_results <- c()
    # Iterate over different numbers of PCs
    for (n_pcs in 1:length(pp0$CA$eig)/3) {  # Adjust 10 to a reasonable number of PCs for your data
      fit <- lda(site ~ ., data=newpca_table[,c(2,4:(n_pcs+4))], na.action="na.omit", CV=TRUE)
      predictions <- fit$class
      predictions
      true_labels <- newpca_table$site
      true_labels
      valid_indices <- !is.na(predictions)
      filtered_predictions <- predictions[valid_indices]
      filtered_true_labels <- true_labels[valid_indices]
      # Compute accuracy only on valid predictions
      accuracy <- mean(filtered_predictions == filtered_true_labels)
      accuracy_results[n_pcs] <- accuracy
    }
    accuracy_results
    #for assignments:
    plda <- predict(object = fit,
                    newdata = juvpca_table[,c(2,4:(n+4))])
    observed_table <- table(plda$class)
    plda$class
    predicted_classes <- data.frame(plda$class)
    TEST=data.frame(ActualSite=juvpca_table$site,Pred=predicted_classes$plda.class)
    Output4Strader=TEST
    row.names(TEST)=juvpca_table$individual
    write.table(TEST,file="JuvPredictionsMoorea.txt",quote=F,row.names = T,col.names = T)
    
    new_df <- TEST %>%
      group_by(ActualSite, Pred) %>%
      summarise(count = n())
    # Print the counts of each class
    print(observed_table)
   
    library(ggalluvial)
    
    colors <- c("#FF5733", "#338DFF", "#33FF57", "#7E33FF", "#FFFF33", "#FFA233", "#FF33B5", "#33FFF9", "#9E6B4D")
    
    new_df <- data.frame(lapply(new_df, function(x) {
      gsub("forereef", " DFR", x)
    }))
    
    new_df <- data.frame(lapply(new_df, function(x) {
      gsub("crest", " SFR", x)
    }))
    
    new_df <- data.frame(lapply(new_df, function(x) {
      gsub("backreef", " BR", x)
    }))
    
    sourceplot=ggplot(data = new_df,
           aes(axis1 = as.factor(Pred), axis2 = as.factor(ActualSite), y = as.numeric(count))) +
      geom_alluvium(aes(fill = Pred)) +
      geom_stratum() +
      geom_text(stat = "stratum",
                aes(label = after_stat(stratum))) +
      scale_x_discrete(limits = c("ActualSite", "Pred"),
                       expand = c(0.15, 0.05)) +
      scale_fill_manual(values=colors, guide="none")+
      theme_void()
    #
    saveRDS(sourceplot,file="JuvSourcePlot.RDS")
    library(dplyr)
    library(reshape2)
    
    # Perform a t-test to compare the counts of the site of interest with the other sites
    lter3=subset(new_df, Pred=="LTER3 DFR")
    nonlter3=subset(new_df, Pred!="LTER3 DFR")
    observed_table
    
    counts_df <- data.frame(
      Site = c("LTER1backreef", "LTER1forereef", "LTER2crest", "LTER2forereef", "LTER3forereef", "LTER5backreef", "LTER5forereef"),
      Counts = c(0, 9, 0, 16, 49, 0, 0)
    )
    
    # Perform the Kruskal-Wallis test
    result <- kruskal.test(Counts ~ Site, data = counts_df)
    result
    # Create a data frame to store the results
    results <- data.frame(Level = character(0), ChiSquared = numeric(0), PValue = numeric(0))
    
    # Get unique levels in ActualSite
    levels_actual_site <- unique(new_df$ActualSite)
    
    # Loop through each level in ActualSite
    for (level in levels_actual_site) {
      # Create subsets of the data for the specific level
      matching_subset <- new_df %>% filter(ActualSite == level, Pred == level)
      non_matching_subset <- new_df %>% filter(ActualSite == level, Pred != level)
      
      # Calculate the counts for each subset
      matching_count <- sum(matching_subset$count)
      non_matching_count <- sum(non_matching_subset$count)
      
      # Create a 2x2 contingency table
      # Create a 1x2 contingency table
      contingency_table <- matrix(c(matching_count, non_matching_count), nrow = 1)
      
      # Perform the chi-squared test
      chi_squared_test <- chisq.test(contingency_table)
      chistat=chi_squared_test$statistic
      if (matching_count < non_matching_count) {
        chistat=chistat*-1
      }
      
      # Store the results in the data frame
      results <- rbind(results, data.frame(Level = level,
                                           ChiSquared = chistat,
                                           PValue = chi_squared_test$p.value))
    }
    
    
    
    # Chi-squared test
    chi_squared_test <- chisq.test(observed_table)
    # Print the test results
    print(chi_squared_test)
    
    
#Comparing juveniles outlier loci means 
    #for lda values
    newpca_table=subset(pca_table, population !="Juv")
    
    #creating the LTER3 plot jumpout
    # exclude=subset(newpca_table, population =="October" & site =="LTER3forereef")
    # newpca_table=subset(newpca_table, !(individual %in% exclude$individual))
    # exclude=subset(newpca_table, population =="May" & site =="LTER3forereef")
    # newpca_table=subset(newpca_table, !(individual %in% exclude$individual))
    # 
    juvpca_table=subset(pca_table, population =="Juv")
    fit <- lda(population ~ ., data=newpca_table[,c(3,4:(n+1))], na.action="na.omit", CV=F, output = "Scatterplot")
    
     plda <- predict(object = fit,
                    newdata = pca_table[,3:(n+1)])
    
    
    prop.lda = fit$svd^2/sum(fit$svd^2)
    dataset <- data.frame(group = pca_table[, 3], lda = plda$x)
    dataset$site=pca_table$site
    
    group_data_October<- subset(dataset, group == "October")
    group_data_May<- subset(dataset, group == "May")
    group_data_Juv<- subset(dataset, group == "Juv")
      Mdistances <- group_data_Juv$LD1 - mean(group_data_May$LD1)
      Odistances <- group_data_Juv$LD1 - mean(group_data_October$LD1)
      mean(group_data_May$LD1)
      #-2.516414
      mean(group_data_October$LD1)
      #3.127543
    df=cbind(Mdistances,Odistances)
    #if df LD1 values are <0 for May or >0 for Oct or greater than Oct then that juvenile falls outside the mean range
    #only one individual does. 
    
    #Putting together predictions for Juvs for outlier DAPC
    
    predicted_classes <- data.frame(plda$class)
    TEST=data.frame(ActualPop=pca_table$population,Pred=predicted_classes$plda.class, site=pca_table$site)
    
    new_df <- TEST %>% subset(ActualPop=="Juv") %>% 
      group_by(ActualPop, Pred, site) %>%
      summarise(count = n())
    # Print the counts of each class
    print(new_df)
  
  # DAPC_plot<- ggplot(dataset, aes(dataset[,1+x_axis], dataset[,1+y_axis], color= group, label=group, shape=group)) + 
  #   geom_enterotype(show.point=show.point, show.label=show.label, show.ellipse=show.ellipse, show.line=show.line, alpha=alpha) +
  #   scale_shape_manual(values = c(rep(c(15,16,17,18),7), 15, 16)) +
  #   theme_cowplot() 
   # labs(x = paste0("LD", x_axis," (", (prop.lda[x_axis]), ")", sep=""),
    #     y = paste0("LD", y_axis, " (", (prop.lda[y_axis]), ")", sep=""))


# Load required libraries
library(ggplot2)
library(dplyr)
library(ggridges)
# Calculate density values for each group
dataset=subset(dataset, group !="Unknown")

dataset$newgrp=ifelse(dataset$group=="October" & dataset$site=="LTER3forereef","Post-mortality LTER3",
                      ifelse(dataset$group=="May" & dataset$site=="LTER3forereef","Pre-mortality LTER3",dataset$group))

# Create a scatter plot for DAPC with density on the y-axis with LTER3 popped out
 ggplot(dataset, aes(x = dataset$LD1, y = 1)) +
  geom_density_ridges(aes(fill = newgrp), alpha = 0.6,) +
   scale_fill_manual(values = c("pink", "orange", "darkblue","#00FFFF","red")) +
   #scale_fill_manual(values = c("pink", "orange", "darkblue"),name=NULL) +
  labs(x = "Discriminant Function", y = "Density") +  theme_classic(base_size = 30) +
  theme()+theme(aspect.ratio = 1) #+ facet_wrap(~site, ncol=4)

 ggplot(dataset, aes(x = dataset$LD1, y = 1)) +
   geom_density_ridges(aes(fill = group), alpha = 0.6,) +
   scale_fill_manual(values = c("pink", "orange", "darkblue")) +
   #scale_fill_manual(values = c("pink", "orange", "darkblue"),name=NULL) +
   labs(x = "Discriminant Function", y = "Density") +  theme_classic(base_size = 30) +
   theme()+theme(aspect.ratio = 1) #+ facet_wrap(~site, ncol=4)
 

#Same plot but facet by site
 dataset$site=pca_table$site
 dataset=subset(dataset, site !="LTER5backreef")
 ggplot(dataset, aes(x = dataset$LD1, y = 1)) +
   geom_density_ridges(aes(fill = group), alpha = 0.6,) +
   scale_fill_manual(values = c("pink", "orange", "darkblue"),name=NULL) +
   labs(x = "Discriminant Function", y = "Density") +  theme_classic(base_size = 30) +
   theme()+theme(aspect.ratio = 1) + facet_wrap(~site, ncol=4)
 
 #Same plot but facet by site
 dataset$site=pca_table$site
 dataset=subset(dataset, site !="LTER5backreef")
 ggplot(dataset, aes(x = dataset$LD1, y = 1)) +
   geom_density_ridges(aes(fill = newgrp), alpha = 0.6,) +
   scale_fill_manual(values = c("pink", "orange", "darkblue","green")) +
   labs(x = "Discriminant Function", y = "Density") +  theme_classic(base_size = 30) +
   theme()+theme(aspect.ratio = 1) + facet_wrap(~site, ncol=4)


 
 #############SIB STUFF################
 
 files=list.files(pattern="bamscl.*sibs$")
 files
 PCoAs=list()
for (bam in 1:length(files)){
  bam=3
bams=data.frame(read.table(files[bam], header = FALSE)) # list of bam files
 colnames(bams)<- "Adapter"
 #removing directory string
 bams <- data.frame(lapply(bams, function(x) {
   gsub(".*/", "", x)
 }))
 
 #extracting sample name from .bam file
 samplenames=word(bams$Adapter,1,sep = "\\_S")
 #removing _'s from some of the samples
 samplenames <- lapply(samplenames, function(x) {
   gsub("_", "", x)
 })
 bams$Sample=unlist(samplenames)
 #rename two samples to match metadata
 bams$Sample=str_replace(bams$Sample, "M1B1$", "M1B01")
 bams$Sample=str_replace(bams$Sample, "N3F18$", "N3F18J")
 
 #read in meta data
 samples=readRDS(file="../forsubbamscl.RData")

 #merge bams and metadata
 new<-merge(bams,samples, by ="Adapter")
 new = new[!duplicated(new$Sample.y),] #some samples were duplicated unnecessarily for some reason
 #new = new[order(new[,'Adapter.x']),]
 #new = new[!duplicated(new$Adapter),] #some samples were duplicated unnecessarily for some reason
 #new<-new[ order(match(new$Adapter.x, bams$Adapter)), ] #MAKE SURE ORDER FOR THIS FILE MATCHES ORDER FOR BAM FILE USED FOR IBS CREATION OTHERWISE DENDO WILL BE WRONG
 new$Site=str_replace(new$Site, "/Pile", "")
 new$Habitat=str_replace(new$Habitat, "Unknown", "forereef")
 new$Site=str_replace(new$Site, "Unknown", "LTER5")
 new$Month_sampled[is.na(new$Month_sampled) ] <- "Juv"
 new$Site.Habitat.Month=paste0(new$Site,new$Habitat,new$Month_sampled)
 new$Site.Habitat=paste0(new$Site,new$Habitat)
 
 #All loci all sampes
 co = as.matrix(read.table(paste0(files[bam],".covMat"))) # covariance based on single-read resampling

 #co = as.matrix(read.table("ok.covar")) # covariance by ngsCovar
 dimnames(co)=list(bams$Adapter,bams$Adapter)
 
 pp0=capscale(as.dist(1-cov2cor(co))~1) # PCoA

 plot(pp0$CA$eig) 
 
 kmeans.axes <- sum(diff(summary(eigenvals(pp0))[2,])*-1 > 0.01) #based on eigen plot
 kmeans.axes
 #0
 # This function finds the total number of eigenvectors in which the difference between last eigenvalue and the one preceding is greater than 0.01
 # Standardized method of selecting eigenvectors that contribute the most explanatory power in clustering
 
 axes2plot=c(1,2)  
 #quartz()
 cc=pp0
 
 cmd=pp0
 cmd$CA$u
 
 
 pca_s <- as.data.frame(cmd$CA$u[,axes2plot])
 #
 
 row.names(pca_s)=str_replace(row.names(pca_s),"X","")
 pca_s$Adapter=row.names(pca_s)
 
 new2<-merge(pca_s,new, by ="Adapter")

 eigenvals=pp0$CA$eig
 #calc variance explained
 varexplained=eigenvals/sum(eigenvals)
 
 e <- eigen(co)
 eigenvals=e$values
 varexplained=eigenvals/sum(eigenvals)
 
 plot(e$vectors[,1:2],lwd=2,ylab="PC 2",xlab="PC 2",main="Principal components",col=rep(1:3,each=10),pch=16)
 varexplained=eigenvals/sum(eigenvals)
 bam_name <- files[bam]  # or extract from filename if better
 
 PCoAs[bam_name] <- list(new2)
 PCoAs[bam] <- list(new2)
 
}
 
 
 names(PCoAs)
 #May Only
 new2=PCoAs[["maybamscl_nohalfsibs_nosibs"]]
 
 new2$Month_sampled=as.factor(new2$Month_sampled)
 new2=dplyr::filter(new2, Site.Habitat !="Unknown")
 new2 <- dplyr::filter(new2, !(Month_sampled == "October"))
 
 new2$Month_sampled <- droplevels(new2$Month_sampled)
 levels(new2$Month_sampled)
 levels(as.factor(new2$Site.Habitat))
 
 shapes=rep(c(13,3:4,6,21:23,11,24))
 shapes
 oranges=rep("orange",9)
 shapes
 shape_names <- c("LTER1backreef", "LTER1crest", "LTER2forereef", "LTER2crest", "LTER3crest",
                  "LTER1forereef", "LTER3forereef", "LTER5backreef", "LTER5forereef")
 names(shapes)=shape_names
 shapes
 levels(as.factor(new2$Site.Habitat))
 Mayplot=ggplot(new2, aes(MDS1, MDS2)) +
   geom_point(data=new2, size=3.5,aes(fill=as.factor(Site.Habitat), shape=as.factor(Site.Habitat)))  +
   theme_classic(base_size = 26) +
   scale_fill_manual(values = oranges,labels=names(shapes),name="Group",drop=TRUE,limits = force) +
   #stat_ellipse(geom="polygon",level=.75, alpha = 1/3, aes(fill = as.factor(new2$Month_sampled)))+
   #scale_fill_manual(values=c('#f0f9e8','#bae4bc','#7bccc4','#43a2ca','#0868ac'))+
   xlab(paste0("PCo1 (",formatC((varexplained[1]*100), digits = 2, format = "f"),"%)")) +
   ylab(paste0("PCo2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
   scale_shape_manual(values=shapes, labels=names(shapes),name="Group") +
  xlim(c(-0.55,0.2))+ylim(c(-0.2,0.4))+guides(color=FALSE)+theme(aspect.ratio = 1) #+guides(color=FALSE) 
 ggplotly(Mayplot)
 Mayplot
 ggsave("MayOnly_nohalfsibs_nosibs_PCAplot_v1.pdf", plot = Mayplot, width = 7.5, height = 5.5)  # Adjust width/height as needed
 
##May Oct Juv
 
 new2=PCoAs[["MayOctJuvOnly_nohalfsibsbamscl_nosibs"]]
 
 p2=ggplot(new2, aes(MDS1, MDS2)) +
   geom_point(data=new2, size=3.5,aes(fill=as.factor(Month_sampled), shape=as.factor(Site.Habitat)))  +
   theme_classic(base_size = 26) +
   #geom_text(hjust=0, vjust=0, label=new2$Adapter)+
   # geom_polygon(data = hull_cyl, alpha = 0.5, color='black')+
   scale_fill_manual(values = c("pink", "orange", "darkblue"),name="") +
   #stat_ellipse(geom="polygon",level=.75, alpha = 1/3, aes(fill = as.factor(new2$Month_sampled)))+
   xlab(paste0("PCo1 (",formatC((varexplained[1]*100), digits = 2, format = "f"),"%)")) +
   ylab(paste0("PCo2 (",formatC((varexplained[2]*100), digits = 2, format = "f"),"%)")) +
   scale_shape_manual(values = rep(21:26),name="")+
   # scale_fill_manual(values=colors, name="", labels = c("October", "May", "Juv", guide=FALSE)) #
   #scale_colour_manual(values=colors, name="", guide=FALSE) #
   #scale_fill_manual(values=colors, name="") #
   xlim(c(-0.2,0.25))+ylim(c(-0.3,0.15))+theme(aspect.ratio = 1)+
   #xlim(c(-0.25,0.65))+ylim(c(-0.25,0.65))+
   #guides(colour=FALSE)+
   guides(fill = guide_legend(override.aes = list(alpha = 1))) 
 ggplotly(p2)
 #N1F8J/N2F36J and N5F29J/N5F26J possible sibling pairs 
 p2
 ggsave("MayOctJuvOnly_nohalfsibsbamscl_nosibss_PCAplot_v1.pdf", plot = p2, width = 7.5, height = 5.5)  # Adjust width/height as needed
 

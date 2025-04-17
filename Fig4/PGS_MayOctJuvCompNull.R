setwd("/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/")

library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(R.utils)
library(tidyr)
library(foreach)
library(doParallel)
library(ggplot2)
library(stats)

setDTthreads(0)

#Generating PGS from all samples to be used for May/Oct
Betafile=fread(file="MayJuv.plink.all.MayJuv.glm.logistic")
header=fread(file="MayJuv.plink.all.MayJuv.glm.logistic") #
HEAD=colnames(header)
colnames(Betafile)=HEAD
samples=readRDS(file="./forsubbamscl.RData")
samples <- data.frame(sapply(samples, function(x) {
  gsub("\\/Pile", "", x)
}))
samples$site.hab.month=paste(samples$Site,samples$Habitat,samples$Month_sampled,sep=".")
#read in bams
bams=data.frame(read.table("MayJuvbamscl4vcf", header = FALSE)) # list of bam files
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

drop=ncol(fread(file="MayJuvAllLoci.geno.gz"))
Genofile=fread(file="MayJuvAllLoci.geno.gz",drop = drop)
adaptnames=as.character(bams$Adapter)
colnames(Genofile)[c(3:ncol(Genofile))]=adaptnames


directory_path <- "/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/"  # Replace with the actual path
# Get the list of file names with the specified pattern
#Just the juveniles
file.names<-basename(list.files(path = directory_path,pattern = glob2rx("*NAbamscl.txt"))) #just taking one lane to get proportion
file.names

#randomly sample loci to create null distribution
allloci=data.frame(Betafile$ID)
MayOctBeta=grep("^MayOctBeta=",commandArgs()) 
if(length(MayOctBeta)>0) { MayOctBeta=as.character(sub("^MayOctBeta=","", commandArgs()[MayOctBeta])) } else { stop ("specify beta file to obtain number of loci to randomly sample") } 

#MayOctBeta="MayOct.plink.all.PrePost.glm.logistic.0.01"
pvalue <- str_extract(MayOctBeta, "\\d+\\.\\d+")
print(pvalue)

num_samples <- 500
sample_size=nrow(read.table(file=MayOctBeta))

# Create a list to store the 500 data frames
sampled_dfs <- vector("list", num_samples)

# Set seed for reproducibility
set.seed(42)

# Loop to create 100 random samples
n=nrow(Betafile)
for (i in 1:num_samples) {
  sample_indices <- sample(1:n, sample_size, replace = FALSE)
  sampled_dfs[[i]] <- Betafile[sample_indices, ]
}


randompgscomp=list()
for (ransam in 1:length(sampled_dfs)){
  #ransam=2
pgscomp=list()
#foreach(i = 1:length(file.names), .packages = c("data.table", "dplyr", "tidyr", "stringr")) %dopar% {
for(i in 1:length(file.names)){
  #i=2
  df <- fread(file.path(directory_path, file.names[i]), header = F)
  colnames(df)[1]="Adapter"
  df_name <- gsub("bamscl.txt", "", file.names[i]) 
  print(df_name)
  df <- data.frame(lapply(df, function(x) {
    gsub(".*/", "", x)
  }))
  newdf<-subset(merge(df,samples, by ="Adapter"),site.hab.month==df_name)
  adapter_columns <- df$Adapter
  tryCatch({
    newgeno <- cbind(Genofile$V1, Genofile$V2, Genofile[, ..adapter_columns])
    # Continue with the rest of your code
  }, error = function(e) {
    # Print the error message
    print(paste("Error:", e))
    # Continue with the loop
    return()
  })
  print(nrow(newgeno))
  
  listofdfs<- list()
  #listofdfs <- foreach(SAMPLE = 3:ncol(newgeno), .combine = c) %dopar% {
  for(SAMPLE in 3:ncol(newgeno)){
    #SAMPLE=3
    Adapter=colnames(newgeno)[SAMPLE]
    #Sample=newgeno %>% select("V1","V2",Adapter)
    Sample = newgeno[, c("V1", "V2", Adapter), with = FALSE]
    Sample.sep = Sample %>% separate(Adapter, c("A", "B"), sep = cumsum(c(1, 1)))
    Sample.sep$V3=paste(Sample.sep$V1,Sample.sep$V2,sep="_")
    #print(Adapters(Sample.sep))
    # Create logical vectors for conditions
    beta_file=Betafile
    Sample.sep <- Sample.sep[Sample.sep$V3 %in% beta_file$ID, ]
    beta_file <- beta_file[beta_file$ID %in% Sample.sep$V3, ]
    sigloci=data.frame(sampled_dfs[[ransam]])
    #colnames(sigloci)="ID"
    beta_file <- beta_file[beta_file$ID %in% sigloci$ID, ]
    Sample.sep <- Sample.sep[Sample.sep$V3 %in% beta_file$ID, ]
    cond1 <- Sample.sep$V3 == beta_file$ID
    cond2_A <- Sample.sep$A == beta_file$ALT
    cond2_B <- Sample.sep$B == beta_file$ALT
    
    # Compute Amult and Bmult using vectorized operations
    Sample.sep$Amult <- as.numeric(cond1 & cond2_A)
    Sample.sep$Bmult <- as.numeric(cond1 & cond2_B)
    
    # Compute mult and addbeta using vectorized operations
    Sample.sep$mult <- Sample.sep$Amult + Sample.sep$Bmult
    Sample.sep$addbeta <- Sample.sep$mult * as.numeric(beta_file$T_STAT)
    sumbeta=sum(Sample.sep$addbeta, na.rm = TRUE)
    #sumbetadf=cbind(Adapter,sumbeta)
    print("finished summing betas total:")
    print(sumbeta)
    listofdfs[[Adapter]] <- data.frame(Adapter = Adapter, sumbeta = sumbeta)
  }
  combined_df <- as.data.table(do.call(rbind, listofdfs))
  colnames(combined_df)=c("Adapter","sumbeta")
  #dflist[[df_name]]=newdf
  pgsdf_real=as.data.table(cbind(newdf$Adapter,newdf$BleachingScoreMay))
  colnames(pgsdf_real)=c("Adapter","BleachingScoreMay")
  combined_df$Adapter <- unlist(combined_df$Adapter)
  pgsdf_real$Adapter <- unlist(pgsdf_real$Adapter)
  pgsdf <- pgsdf_real[combined_df, on = "Adapter"]
  pgscomp[[df_name]]=pgsdf
  print("finished running:")
  print(df_name)
}

combined_list <- list()
# Iterate through each named data frame in your list
for (i in seq_along(pgscomp)) {
  df_name <- names(pgscomp)[i]
  df <- pgscomp[[i]]
  
  # Add a new column with the corresponding df name
  df$site.hab.month <- df_name
  
  # Append the modified df to the combined_list
  combined_list[[i]] <- df
}

# Combine the data frames in the combined_list using rbind
combined_df <- do.call(rbind, combined_list)
combined_df$site.hab <- sub("^([^\\.]+\\.[^\\.]+)\\..*$", "\\1", combined_df$site.hab.month)
combined_df$month_sampled <- gsub("^[^.]*\\.[^.]*\\.(.*)$", "\\1", combined_df$site.hab.month)
combined_df$month_sampled <- gsub("NA", "Juv", combined_df$month_sampled)
combined_df$hab <- sub(".*\\.(.*)\\..*", "\\1", combined_df$site.hab.month)


randompgscomp[[ransam]]=combined_df
}

outfile=paste0("pgscomp_ran",MayOctBeta,".RData")
saveRDS(randompgscomp,file=outfile)
#stopCluster()
quit(save="no")

BIGpgscomp=readRDS(file="pgscomp_ranMayOct2FR.plink.all.PrePost.glm.logistic.0.1.RData")


# Install and load necessary packages
library(dplyr)
library(purrr)
# Generate 500 sets of 10,000 random loci
# Calculate PGS for each set of random loci and store in a list
pgs_random_sets <- lapply(BIGpgscomp, function(df) df$sumbeta)

#


# Flatten the list of PGS values into a single vector
all_pgs <- abs(na.omit(unlist(pgs_random_sets)))
all_pgs

# Calculate 95th percentile for all PGS values across all random sets
percentile_95 <- quantile(all_pgs, 0.95)
percentile_95


# Check if your PGS falls within the 95th percentile
pgscomp=readRDS(file="pgscomp_MayOctJuv.RData")

pgs_real=lapply(pgscomp, function(df) df$sumbeta)
pgs_real_df=data.frame(abs(na.omit(unlist(pgs_real))))
colnames(pgs_real_df)="PGS"

is_within_95th_percentile <- any(pgs_your_individuals >= percentile_95)

#plotting:
df=data.frame(all_pgs)

ggplot(df,aes(x=all_pgs))+
  geom_density(aes(y = after_stat(density / sum(density))),color="darkblue", fill="lightblue")+
  #xlim(-.000000038,0.0000006)+
  # Assuming df2 is another dataframe with the data for the second density plot
  geom_density(data=pgs_real_df, aes(x=PGS, y = after_stat(density / sum(density))), color="red", fill="pink",alpha=.5) + # Add another density plot
  ylab("Proportion")+xlab("PGS")+
  geom_vline(aes(xintercept=percentile_95),color="black", linetype="dashed", size=.5)+
  theme_classic()

#Alternatively take the mean of every set
mean_list <- lapply(pgs_random_sets, mean)
all_pgs <- abs(na.omit(unlist(mean_list)))
percentile_95 <- quantile(all_pgs, 0.95)


df=data.frame(all_pgs)

ggplot(df,aes(x=all_pgs))+
  geom_density(aes(y = after_stat(density / sum(density))),color="darkblue", fill="lightblue")+
  xlim(-100,2000)+
  # Assuming df2 is another dataframe with the data for the second density plot
  geom_density(data=pgs_real_df, aes(x=PGS, y = after_stat(density / sum(density))), color="black", fill="#F7CBD3",alpha=.5) + # Add another density plot
  ylab("Proportion")+xlab("PGS")+
  geom_vline(aes(xintercept=percentile_95),color="black", linetype="dashed", size=.5)+
  theme_classic()







  
###Predictability

#linear modelhttps://scc-ondemand2.bu.edu/rnode/scc-c11/32414/graphics/97645966-ed3c-4db9-ac53-02f7cb7b4ab5.png
fitPCs <- lm(as.numeric(dfcovar$trait) ~ as.numeric(dfcovar$PC1)+as.numeric(dfcovar$PC2))
fitPCsEnv <- lm(as.numeric(dfcovar$trait) ~ as.numeric(dfcovar$PC1)+as.numeric(dfcovar$PC2)+as.numeric(dfcovar$depth)+as.numeric(dfcovar$SA))
fitPCsEnvSym<- lm(as.numeric(dfcovar$trait) ~  as.numeric(dfcovar$PropA)+ as.numeric(dfcovar$PropC)+as.numeric(dfcovar$depth)+as.numeric(dfcovar$PC1)+as.numeric(dfcovar$PC2)+as.numeric(dfcovar$SA))
fitPCsEnvBeta<- lm(as.numeric(dfcovar$trait) ~  as.numeric(dfcovar$sumbeta)+as.numeric(dfcovar$depth)+as.numeric(dfcovar$PC1)+as.numeric(dfcovar$PC2)+as.numeric(dfcovar$SA))

fitfull <- lm(as.numeric(dfcovar$trait) ~ as.numeric(dfcovar$sumbeta) + as.numeric(dfcovar$PropA)+ as.numeric(dfcovar$PropC)+as.numeric(dfcovar$depth)+as.numeric(dfcovar$PC1)+as.numeric(dfcovar$PC2)+as.numeric(dfcovar$SA))









#
combined_df_meta=merge(combined_df,samples, by="Adapter")
filtered_df <- combined_df_meta[!combined_df_meta[, duplicated(.SD, by = c("Adapter", "site.hab.month.x"))]]


df2$month_sampled=factor(df2$month_sampled,levels=c("Pre-mortality 2019","Post-mortality 2019","Juveniles 2021"))

ggplot(df2, aes(x = month_sampled, y = sumbeta)) +
  geom_point(aes(fill=month_sampled),shape = 21, size = 4) +
  geom_errorbar(aes(ymin = sumbeta - sd, ymax = sumbeta + sd), width = 0.2) +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~site.hab,scales = "free_x",labeller = labeller(site.hab = custom_labels))+
  theme_classic() +
  theme(legend.position="bottom",legend.title = element_blank())  +
    scale_fill_manual(values=c("orange", "darkblue", "pink"))+
xlab("")+ylab("Polygenic Score")+ theme(axis.title.x=element_blank(),
                                        axis.text.x=element_blank(),
                                        axis.ticks.x=element_blank())

df2 <- data_summary(combined_df, varname="sumbeta", 
                    groupnames=c("hab", "month_sampled"))

df2$month_sampled=factor(df2$month_sampled,levels=c("May","October","Juv"))

ggplot(df2, aes(x = month_sampled, y = sumbeta)) +
  geom_point(aes(fill=month_sampled),shape = 21, size = 4, show.legend = F) +
  geom_errorbar(aes(ymin = sumbeta - sd, ymax = sumbeta + sd), width = 0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~hab, scales = "free_x")+theme_classic() +scale_fill_manual(values=c("orange", "darkblue", "pink"))+
  xlab("")+ylab("Polygenic Score")


#testing significance


res=aov(sumbeta~month_sampled*site.hab ,data=combined_df)
summary(res)

res=aov(sumbeta~month_sampled ,data=subset(combined_df,hab=="forereef"))
summary(res)



res=aov(sumbeta~month_sampled ,data=subset(combined_df,site.hab=="LTER1.forereef"))
summary(res)
TukeyHSD(res, which = "month_sampled")

res=aov(sumbeta~month_sampled ,data=subset(combined_df,site.hab=="LTER2.forereef"))
summary(res)
TukeyHSD(res, which = "month_sampled")

res=aov(sumbeta~month_sampled ,data=subset(combined_df,site.hab=="LTER3.forereef"))
summary(res)
TukeyHSD(res, which = "month_sampled")

res=aov(sumbeta~month_sampled ,data=subset(combined_df,site.hab=="LTER5.forereef"))
summary(res)
TukeyHSD(res, which = "month_sampled")



Mayonly=subset(combined_df, month_sampled=="May")
p <- ggplot(Mayonly, aes(y = as.numeric(BleachingScoreMay), x=sumbeta) ,group=site.hab) + 
  geom_point(aes(fill=site.hab, group=site.hab),position = position_dodge(0.25), pch = 21, cex = 2,show.legend = FALSE) + theme_classic() +
  #scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#dd0000"), name = "", labels = c("Symbiodinium", "Cladocopoum", "Durusdinium")) +
  #xlab("Bleaching Score") + ylab("Proportion of Reads") +
  #theme(text = element_text(size = 15), legend.text = element_text(face = "italic")) +
  geom_smooth(aes(group=site.hab, fill=site.hab),method = 'lm', color = "black", show.legend = FALSE)+
  xlab("Polygenic Score") +ylab("Bleaching Score")
p +facet_wrap(~ site.hab) 


# Plot
########


rsq=data.frame(rsquared_values_model)
sd(rsq$rsquared_values)
mean(rsq$rsquared_values)
sum=summary(fit)
sum$r.squared
mean_rsquared <- mean(rsquared_values)
mean_rsquared
#1st column cov
#[1] 0.05116702
#7 covs, 0.000001 p val
# [1] 0.03093542
#7 covs, 0.0001 p val
#0.05609546

sd_rsquared <- sd(rsquared_values)
sd_rsquared
#0.06499134
#7 covs, 0.000001 p val
#[1] 0.06072632
#7 covs, 0.0001 p val
#[1] 0.0840278




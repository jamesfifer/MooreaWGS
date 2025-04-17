setwd("/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/")
library(ggpubr)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(R.utils)
library(tidyr)
library(foreach)
library(doParallel)
library(ggplot2)

setDTthreads(0)




##My own GWAS
directory_path <- "/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/MayOct2PGS_holdout/"  # Replace with the actual path
# Get the list of file names with the specified pattern
file.names<-basename(list.files(path = directory_path,pattern = glob2rx("*holdout*.psam"))) #just taking one lane to get proportion
file.names


samples=subset(readRDS(file="./forsubbamscl.RData"))

#read in bams
bams=data.frame(read.table("MayOctbamscl4vcf", header = FALSE)) # list of bam files
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

new<-merge(bams,samples, by ="Adapter")
new=subset(new, Sample.y!="M5B9")

new$NewAdapter=new$Adapter
to_modify <- new$Month_sampled == "October" & duplicated(new$Adapter, fromLast = TRUE)| duplicated(new$Adapter) & new$Month_sampled != "May"
# Apply the transformation to the identified rows
new$NewAdapter[to_modify] <- paste("Oct", new$Adapter[to_modify], sep=".")
new<-new[ order(match(new$NewAdapter, bams$Adapter)), ]


Genofile=fread(file="MayOctAllLoci.geno.gz",drop = c("V234")) #last column always produces NA for some reason so dropping
adaptnames=as.character(bams$Adapter)
colnames(Genofile)[c(3:ncol(Genofile))]=adaptnames

#num_cores <- 18
#registerDoParallel(cores = num_cores)
#pvalue 
pval=grep("^pval=",commandArgs()) 
if(length(pval)>0) { pval=as.numeric(sub("^pval=","", commandArgs()[pval])) } else { stop ("specify pval according to file format (e.g. pval=0.0001)") } 
#pval=0.2
pval <- format(pval, scientific = FALSE)#was weirding converting to scientific notation 
print(pval)


pgscomp=list()
#foreach(i = 1:length(file.names), .packages = c("data.table", "dplyr", "tidyr", "stringr")) %dopar% {
for(i in 1:length(file.names)){
  #i=1
  df <- fread(file.path(directory_path, file.names[i]), header = TRUE)
  colnames(df)[2]="Adapter"
  df_name <- gsub("MayOct2.holdout.|\\.psam", "", file.names[i]) 
  beta_filename <- paste0(directory_path, "MayOct2.train.", df_name, ".",pval,".beta.clumped")
  if (!file.exists(beta_filename)) {
    print(paste("Beta file not found for:", df_name))
    next  # Skip to the next iteration if beta file is not found
  }
  # Check if the file is empty
  if (file.info(beta_filename)$size == 0) {
    print(paste("Empty beta file found for:", df_name))
    next  # Skip to the next iteration if beta file is empty
  }
  beta_file <- fread(file = beta_filename)
  print(df_name)
  newdf<-merge(df,new, by.x="Adapter",by.y ="NewAdapter")
  adapter_columns <- df$Adapter
  newgeno <- cbind(Genofile$V1,Genofile$V2, Genofile[, ..adapter_columns])
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
    Sample.sep <- Sample.sep[Sample.sep$V3 %in% beta_file$ID, ]
    beta_file <- beta_file[beta_file$ID %in% Sample.sep$V3, ]
    
    cond1 <- Sample.sep$V3 == beta_file$ID
    cond2_A <- Sample.sep$A == beta_file$ALT
    cond2_B <- Sample.sep$B == beta_file$ALT
    
    # Compute Amult and Bmult using vectorized operations
    Sample.sep$Amult <- as.numeric(cond1 & cond2_A)
    Sample.sep$Bmult <- as.numeric(cond1 & cond2_B)
    
    # Compute mult and addbeta using vectorized operations
    Sample.sep$mult <- Sample.sep$Amult + Sample.sep$Bmult
    Sample.sep$addbeta <- Sample.sep$mult * as.numeric(beta_file$OR)
    
    sumbeta=sum(Sample.sep$addbeta)
    #sumbetadf=cbind(Adapter,sumbeta)
    print("finished summing betas total:")
    print(sumbeta)
    listofdfs[[Adapter]] <- data.frame(Adapter = Adapter, sumbeta = sumbeta)
  }
  combined_df <- as.data.table(do.call(rbind, listofdfs))
  colnames(combined_df)=c("Adapter","sumbeta")
  #dflist[[df_name]]=newdf
  newdf$trait=ifelse(newdf$Month_sampled=='May', "1", "0") #need to add for May/Oct comp
  pgsdf_real=as.data.table(cbind(newdf$Adapter,newdf$trait))
  colnames(pgsdf_real)=c("Adapter","trait")
  combined_df$Adapter <- unlist(combined_df$Adapter)
  pgsdf_real$Adapter <- unlist(pgsdf_real$Adapter)
  pgsdf <- pgsdf_real[combined_df, on = "Adapter"]
  pgscomp[[df_name]]=pgsdf
  print("finished running:")
  print(df_name)
}

outfile=paste("MayOct2.plink.log.pgscomp_allcovs.",pval,".RData",sep="")
saveRDS(pgscomp,file=outfile)
#stopCluster()
quit(save="no")


#First run only one covariatie
#pgscomp=readRDS(file="MayOct.plink.log.pgscomp_allcovs.0.001.RData")

##First testing whether including PCs 1 and 2 actually make a difference
withPCS=fread(file="MayOct2.plink.all.PrePost.glm.logistic")
noPCS=fread(file="MayOct2.plink.all.noPCs.PrePost.glm.logistic")
withPCS=withPCS %>% subset(TEST=="ADD" & !is.na(P))
noPCS=noPCS %>% subset(TEST=="ADD" & !is.na(P))
colnames(withPCS)=paste0("withPCS",colnames(withPCS))
df=data.frame(cbind(withPCS, noPCS))
head(df)
p_cor <- cor(as.numeric(df$withPCSP), as.numeric(df$P), method = "spearman")
p_corp <- cor.test(as.numeric(df$withPCSP), as.numeric(df$P), method = "spearman")
#^pvalue of 0, so will report at p<0.001

beta_cor <- cor(df$withPCSOR, df$OR)
beta_cor


set.seed(123)  # For reproducibility
sample_df <- df[sample(1:nrow(df), 50000), ]  # Sample 10,000 points

p1 <- ggplot(sample_df,aes(y=withPCSP, x=P)) + 
  geom_point(pch=21, cex=2) + theme_classic()+
  geom_smooth(method='lm',color = "red",show.legend=F )+
  ylab("GWAS p-values (including PC1 & PC2 covariates")+
  xlab("GWAS p-values (excluding PC1 & PC2 covariates")+
  annotate("text", x = 0.05, y = .95, label = paste0("r = ",signif(p_cor,3)))+
  annotate("text", x = 0.05, y = .91, label = paste0("p < 0.001 "))
p1
p2 <- ggplot(sample_df,aes(y=withPCSOR, x=OR)) + 
  geom_point(pch=21, cex=2) + theme_classic()+
  geom_smooth(method='lm',color = "red",show.legend=F )+
  ylab("GWAS beta-values (including PC1 & PC2 covariates)")+
  xlab("GWAS beta-values (excluding PC1 & PC2 covariates)")+
  annotate("text", x = 1, y = 20, label = paste0("r = ",signif(beta_cor,3)))+
  annotate("text", x = 1, y = 18, label = paste0("p < 0.001 "))

p2

p_cor <- cor(as.numeric(withPCS$P), as.numeric(noPCS$P), method = "spearman")
beta_cor <- cor(withPCS$OR, noPCS$OR)

directory_path <- "/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/" 
# Get the list of file names with the specified pattern
file.names<-basename(list.files(path = directory_path,pattern = glob2rx("MayOct2.plink.log.pgscomp_allcovs.*.RData"))) 
file.names
allpgscomps=list()
for (i in 1:length(file.names)) {
  #i=5
  pgscomp <- readRDS(file.path(directory_path, file.names[i]))
  df_name <- gsub("MayOct2.plink.log.pgscomp_allcovs|\\.RData", "", file.names[i]) 
  rsquared_values <- numeric(length(pgscomp))
  lmrsquared_values <- numeric(length(pgscomp))
  for (j in 1:length(pgscomp)) {
   # j=2
    df <- pgscomp[[j]]
    fitbetalm <- lm(as.numeric(df$trait) ~ as.numeric(df$sumbeta))
    lmrsquared_values[j] <- summary(fitbetalm)$r.squared
  }
  allpgscomps[[df_name]] <- lmrsquared_values
  #rsquared_values_null[i] <- summary(fit_null)$r.squared
}
means <- numeric(length(allpgscomps))
sds <- numeric(length(allpgscomps))
df_names <- names(allpgscomps)

for (i in 1:length(allpgscomps)) {
  values <- unlist(allpgscomps[[df_names[i]]])
  means[i] <- mean(values)
  sds[i] <- sd(values)
}

# Create a dataframe to store the results
summary_df <- data.frame(DataFrame = df_names, Mean = means, SD = sds)

ggplot(summary_df, aes(x = DataFrame, y = Mean)) +
  geom_point(size = 4, color = "skyblue") +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
  theme_minimal()+ylab("R2")+xlab("pval cutoff")
#updated 0.01 has highest r2


library(pROC)

# Get the list of file names with the specified pattern
file.names<-basename(list.files(path = directory_path,pattern = glob2rx("MayOct2.plink.log.pgscomp_allcovs.*.RData"))) 
file.names
allpgscomps <- list()
roc_curves <- list()  # List to store ROC curves

for (i in 1:length(file.names)) {
  #i=2
  pgscomp <- readRDS(file.path(directory_path, file.names[i]))
  df_name <- gsub("MayOct.plink.log.pgscomp_allcovs|\\.RData", "", file.names[i]) 
  rsquared_values <- numeric(length(pgscomp))
  roc_data_list <- list()  # List to store ROC data
  
  for (j in 1:length(pgscomp)) {
    #j=8
    df <- pgscomp[[j]]
    fitbeta <- glm(as.numeric(df$trait) ~ as.numeric(df$sumbeta), data = df, family = binomial)
    null_model <- glm(as.numeric(df$trait) ~ 1, data = df, family = binomial)
    
    null_loglik <- logLik(null_model)
    fitted_loglik <- logLik(fitbeta)
    fitted_loglik
    null_loglik
    mcfadden_r2 <- 1 - (fitted_loglik / null_loglik)
    rsquared_values[j] <- mcfadden_r2
    
    # Calculate ROC curve and store it
    if (any(is.na(df$trait))) {
      print("Warning: NAs detected")
    }
    roc_data <- roc(na.omit(as.numeric(df$trait)), predict(fitbeta, type = "response"))
    print(j)
    roc_data_list[[j]] <- roc_data
  }
  
  allpgscomps[[df_name]] <- rsquared_values
  roc_curves[[df_name]] <- roc_data_list
}

means <- numeric(length(allpgscomps))
sds <- numeric(length(allpgscomps))
auc_values <- numeric(length(roc_curves))
df_names <- names(allpgscomps)

for (i in 1:length(allpgscomps)) {
  values <- unlist(allpgscomps[[df_names[i]]])
  means[i] <- mean(values)
  sds[i] <- sd(values)
  
  # Calculate mean AUC for each subset
  auc_values[i] <- mean(sapply(roc_curves[[df_names[i]]], function(roc_data) auc(roc_data)))
}

# Create a dataframe to store the results
summary_df <- data.frame(DataFrame = df_names, Mean_R2 = means, SD_R2 = sds, Mean_AUC = auc_values)

ggplot(summary_df, aes(x = DataFrame)) +
  geom_col(aes(y = Mean_R2), fill = "skyblue", width = 0.7, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = Mean_R2 - SD_R2, ymax = Mean_R2 + SD_R2), width = 0.4, position = position_dodge(width = 0.8)) +
  geom_line(aes(y = Mean_AUC ), color = "red", group = 1, size = 1) +
  geom_point(aes(y = Mean_AUC ), color = "red", size = 4) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "Mean AUC")) +
  theme_minimal() + ylab("R2") + xlab("pval cutoff")
#This code will calculate ROC curves and AUC values for each subset, and then create a combined plot showing both the mean R2 and mean AUC for different p-value cutoffs. The AUC values are scaled to percentage for better visualization on a secondary y-axis.



#Taking p<0.01
pgscomp <- readRDS(file="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/MayOct2.plink.log.pgscomp_allcovs.0.2.RData")
#merge with trait data to get other covariates to compare
#header_line <- readLines("/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/MayOct2plink_covars")
#header <- strsplit(header_line[grep("^#", header_line)], "\\s+")[[1]]  # Split by spaces

# Read the data (skip the header row)
covars <- read.table(file="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/MayOct2plink_covars", 
                     header=T, 
                     comment.char="#")  # Skip the row with #
# Check the data
head(covars)
colnames(covars)[1]="Adapter"

merge_with_covars <- function(df) {
  merged_df <- merge(df, covars, by = "Adapter")
  return(merged_df)
}
# Use lapply to apply the merge_with_new function to each data frame in the list
merged_dataframes_list <- lapply(pgscomp, merge_with_covars)

#I think its better to create a list and then each i is a dataframe with different columns including the different predictors
rsquared_values=list()
mcrsquared_values=list()
roc_data_list=list()
for (i in 1:length(merged_dataframes_list)) {
  #i=1
  df <- merged_dataframes_list[[i]]
  #df$trait=df$BS_normalized
  fitPCs <- lm(as.numeric(df$trait) ~ as.numeric(df$PC1)+as.numeric(df$PC2))
  fitPCsEnv <- lm(as.numeric(df$trait) ~ as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$depth)+as.numeric(df$area))
  fitPCsEnvSym<- lm(as.numeric(df$trait) ~  as.numeric(df$PropA)+ as.numeric(df$PropC)+as.numeric(df$depth)+as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$area))
  fitPCsEnvBeta<- lm(as.numeric(df$trait) ~  as.numeric(df$sumbeta)+as.numeric(df$depth)+as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$area))
  fitfull <- lm(as.numeric(df$trait) ~ as.numeric(df$sumbeta) + as.numeric(df$PropA)+ as.numeric(df$PropC)+as.numeric(df$depth)+as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$area))
  newdf=data.frame(cbind(summary(fitfull)$r.squared,summary(fitPCsEnvSym)$r.squared,summary(fitPCsEnvBeta)$r.squared,summary(fitPCsEnv)$r.squared,summary(fitPCs)$r.squared))
  colnames(newdf)=c("fitfull","fitPCsEnvSym","fitPCsEnvBeta","fitPCsEnv","fitPCs")
  rsquared_values[[i]] <- newdf
  
  fitlogPCs <- glm(as.numeric(df$trait) ~ as.numeric(df$PC1)+as.numeric(df$PC2), data = df, family = binomial)
  fitlogPCsEnv <- glm(as.numeric(df$trait) ~ as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$depth)+as.numeric(df$area), data = df, family = binomial)
  fitlogPCsEnvSym<- glm(as.numeric(df$trait) ~  as.numeric(df$PropA)+ as.numeric(df$PropC)+as.numeric(df$depth)+as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$area), data = df, family = binomial)
  fitlogPCsEnvBeta<- glm(as.numeric(df$trait) ~  as.numeric(df$sumbeta)+as.numeric(df$depth)+as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$area), data = df, family = binomial)
  fitlogfull <- glm(as.numeric(df$trait) ~ as.numeric(df$sumbeta) + as.numeric(df$PropA)+ as.numeric(df$PropC)+as.numeric(df$depth)+as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$area), data = df, family = binomial)
  null_model <- glm(as.numeric(df$trait) ~ 1, data = df, family = binomial)
  null_loglik <- logLik(null_model)
  fitted_loglikPCs <- logLik(fitlogPCs)
  fitted_loglikPCsEnv<- logLik(fitlogPCsEnv)
  fitted_loglikPCsEnvSym<- logLik(fitlogPCsEnvSym)
  fitted_loglikPCsEnvBeta<- logLik(fitlogPCsEnvBeta)
  fitted_loglikfull<- logLik(fitlogfull)

  mcfadden_r2PCs <- 1 - (fitted_loglikPCs / null_loglik)
  mcfadden_r2PCsEnv <- 1 - (fitted_loglikPCsEnv / null_loglik)
  mcfadden_r2PCsEnvSym <- 1 - (fitted_loglikPCsEnvSym / null_loglik)
  mcfadden_r2PCsEnvBeta <- 1 - (fitted_loglikPCsEnvBeta / null_loglik)
  mcfadden_r2full <- 1 - (fitted_loglikfull / null_loglik)
  
  mcrsquared_values[[i]] <- data.frame(cbind(mcfadden_r2PCs,mcfadden_r2PCsEnv,mcfadden_r2PCsEnvSym,mcfadden_r2PCsEnvBeta,mcfadden_r2full))
  
   #roc_dataPCs <- roc(na.omit(as.numeric(df$trait)), predict(fitlogPCs, type = "response"))
   #roc_dataPCsEnv <- roc(na.omit(as.numeric(df$trait)), predict(fitlogPCsEnv, type = "response"))
   #roc_dataPCsEnvBeta <- roc(na.omit(as.numeric(df$trait)), predict(fitlogPCsEnvBeta, type = "response"))
   #roc_dataPCsEnvSym <- roc(na.omit(as.numeric(df$trait)), predict(fitlogPCsEnvSym, type = "response"))
   #roc_datafull <- roc(na.omit(as.numeric(df$trait)), predict(fitlogfull, type = "response"))
  # 
   #roc_data_list[[i]] <- roc_data
  #rsquared_values_null[i] <- summary(fit_null)$r.squared
}
r2df=do.call("rbind",rsquared_values) 
r2df_melt <- melt(do.call("rbind",rsquared_values))

r2df=do.call("rbind",mcrsquared_values) 
r2df_melt <- melt(do.call("rbind",mcrsquared_values))

#auc_values[i] <- mean(sapply(roc_curves[[df_names[i]]], function(roc_data) auc(roc_data)))

ggplot(r2df_melt %>%
         filter(value > 0), aes(x=variable, y=value))+geom_boxplot()+labs(title = " beta cutoff p<0.1",
                                                                y = "R2",x="Model")
#Pretty 
r2df$id= seq_len(nrow(r2df))
r2df_melt <- melt(r2df, id.vars = "id")

p1=ggplot(subset(r2df_melt, variable !="fitPCsEnvBeta"), aes(x=variable, y=value, fill=variable))+
  labs(y = expression("R"^2),x="")+ theme_classic()+theme(text=element_text(size = 15),axis.text.x = element_text(angle = 40, vjust =1, hjust=1))+
  geom_point(alpha = 0.2, color = "white") +  # Adjust alpha for transparency
  geom_line(aes(group = id), alpha = 0.3, color = "grey")+  # Adjust alpha for transparency
  geom_boxplot()+scale_x_discrete(limits = c("fitPCs", "fitPCsEnv", "fitPCsEnvSym", "fitfull")
                                  , labels = c("Genetic PC1 & PC2", "+ Surface Area & Depth", "+ Prop. Symbiont", "+ Survival PGS"))+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), guide="none")+  # Set your custom colors
   scale_y_continuous(breaks = seq(0, 1, by = .25))+ylim(0,1)+
   theme(aspect.ratio=2)
p1
#
p2=ggplot(subset(r2df_melt, variable =="fitPCsEnvBeta" | variable =="fitPCsEnv"), aes(x=variable, y=value, fill=variable))+
  labs(y = expression("R"^2),x="")+ theme_classic()+theme(text=element_text(size = 15))+
  geom_point(alpha = 0.2, color = "white") +  # Adjust alpha for transparency
  geom_line(aes(group = id), alpha = 0.3, color = "grey")+  # Adjust alpha for transparency
  geom_boxplot()+scale_x_discrete(limits = c("fitPCsEnv", "fitPCsEnvBeta")
                                  , labels = c("", "+ Survival PGS"))+
  scale_fill_manual(values = c( "#D55E00","#009E73"), guide="none")+ # Set your custom colors
   scale_y_continuous(breaks = seq(0, 1, by = .25))+ylim(0,1)+
  theme(aspect.ratio=2)
# Plot
figure <- ggarrange(p1,NULL, p2,
                    labels = c("", ""), align="h", nrow=1,widths = c(1, -.5,1))
figure

# Plot

result <- wilcox.test(r2df$fitfull, r2df$fitPCsEnvSym, alternative = "greater")
result[["p.value"]]
#0.05248861
result <- wilcox.test(r2df$fitPCsEnvBeta,r2df$fitPCsEnv,  alternative = "greater")
result[["p.value"]]
#[1] 0.001353182

result <- wilcox.test(r2df$mcfadden_r2full, r2df$mcfadden_r2PCsEnvSym, alternative = "greater")
result[["p.value"]]
#0.05248861
result <- wilcox.test(r2df$mcfadden_r2PCsEnvBeta,r2df$mcfadden_r2PCsEnv,  alternative = "greater")
result[["p.value"]]
#[1] 0.001353182

rsquared_values=list()
for (i in 1:length(merged_dataframes_list)) {
  df <- merged_dataframes_list[[i]]
  #df$trait=df$BS_normalized
  fitPCs <- lm(as.numeric(df$trait) ~ as.numeric(df$PC1)+as.numeric(df$PC2))
  fitPCsEnv <- lm(as.numeric(df$trait) ~ as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$depth)+as.numeric(df$area))
  fitPCsEnvSym<- lm(as.numeric(df$trait) ~  as.numeric(df$PropA)+ as.numeric(df$PropC)+as.numeric(df$depth)+as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$area))
  fitPCsEnvBeta<- lm(as.numeric(df$trait) ~  as.numeric(df$sumbeta)+as.numeric(df$depth)+as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$area))
  fitfull <- lm(as.numeric(df$trait) ~ as.numeric(df$sumbeta) + as.numeric(df$PropA)+ as.numeric(df$PropC)+as.numeric(df$depth)+as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$area))
  newdf=data.frame(cbind(summary(fitfull)$r.squared,summary(fitPCsEnvSym)$r.squared,summary(fitPCsEnvBeta)$r.squared,summary(fitPCsEnv)$r.squared,summary(fitPCs)$r.squared))
  colnames(newdf)=c("fitfull","fitPCsEnvSym","fitPCsEnvBeta","fitPCsEnv","fitPCs")
  rsquared_values[[i]] <- newdf
  #rsquared_values_null[i] <- summary(fit_null)$r.squared
}
r2df=do.call("rbind",rsquared_values) 
r2df_melt <- melt(do.call("rbind",rsquared_values))


ggplot(r2df_melt, aes(x=variable, y=value))+geom_boxplot()+labs(title = " beta cutoff p<0.1",
                                                                y = "R2",x="Model")
######In response to reviewers concern about site effects

#Taking p<0.1
pgscomp <- readRDS(file="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/MayOct.plink.log.pgscomp_allcovs.0.01.RData")
#merge with trait data to get other covariates to compare
#alternatively do GPS
GPS=read.csv(file="../SRC/GPS.csv")
#going to use projected coordinates but lat/long is probably fine 
new$site.habitat=paste(new$Site,new$Habitat, sep=".")
new2=merge(new, GPS[,c(1,4,5)], by.x = "site.habitat", by.y="X")

#adding depth (habitat)
new2$depth=ifelse(new$Habitat=='backreef', "1", ifelse(new$Habitat=='forereef', "10", "5"))


merge_with_covars <- function(df) {
  merged_df <- merge(df, new2, by.x = "Adapter", by.y="NewAdapter")
  return(merged_df)
}
merged_dataframes_list <- lapply(pgscomp, merge_with_covars)

#I think its better to create a list and then each i is a dataframe with different columns including the different predictors
rsquared_values=list()
mcrsquared_values=list()
roc_data_list=list()
for (i in 1:length(merged_dataframes_list)) {
  #i=1
  df <- merged_dataframes_list[[i]]
  count_LTER2 <- sum(df$site.habitat == "LTER2.forereef")
  
  #df$trait=df$BS_normalized
  #fitPCs <- lm(as.numeric(df$trait) ~ as.numeric(df$PC1)+as.numeric(df$PC2))
  fitPCsEnv <- lm(as.numeric(df$trait) ~ as.numeric(df$x)+as.numeric(df$y))
 # fitPCsEnvSym<- lm(as.numeric(df$trait) ~  as.numeric(df$PropA)+ as.numeric(df$PropC)+as.numeric(df$depth)+as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$area))
  fitPCsEnvBeta<- lm(as.numeric(df$trait) ~  as.numeric(df$sumbeta)+as.numeric(df$depth)+as.numeric(df$x)+as.numeric(df$y)+as.numeric(df$area))
  fitBeta<- lm(as.numeric(df$trait) ~  as.numeric(df$sumbeta))
  #fitfull <- lm(as.numeric(df$trait) ~ as.numeric(df$sumbeta) + as.numeric(df$PropA)+ as.numeric(df$PropC)+as.numeric(df$depth)+as.numeric(df$PC1)+as.numeric(df$PC2)+as.numeric(df$area))
  newdf=data.frame(cbind(summary(fitBeta)$r.squared,summary(fitPCsEnvBeta)$r.squared,summary(fitPCsEnv)$r.squared), count_LTER2)
  colnames(newdf)=c("fitBeta","fitPCsEnvBeta","fitPCsEnv","countLTER2")
  rsquared_values[[i]] <- newdf
  
  
  fitlogPCsEnv <- glm(as.numeric(df$trait) ~ as.numeric(df$x)+as.numeric(df$y)+as.numeric(df$depth)+as.numeric(df$area))
  fitlogPCsEnvBeta<- glm(as.numeric(df$trait) ~  as.numeric(df$sumbeta)+as.numeric(df$depth)+as.numeric(df$x)+as.numeric(df$y)+as.numeric(df$area))
  fitlogBeta<- glm(as.numeric(df$trait) ~  as.numeric(df$sumbeta))
  
  null_model <- glm(as.numeric(df$trait) ~ 1, data = df, family = binomial)
  null_loglik <- logLik(null_model)
  fitted_loglikPCsEnv<- logLik(fitlogPCsEnv)
  fitted_loglikPCsEnvBeta<- logLik(fitlogPCsEnvBeta)
  fitted_loglikBeta<- logLik(fitlogBeta)

  mcfadden_r2PCsEnv <- 1 - (fitted_loglikPCsEnv / null_loglik)
  mcfadden_r2PCsEnvBeta <- 1 - (fitted_loglikPCsEnvBeta / null_loglik)
  mcfadden_r2Beta<- 1 - (fitted_loglikBeta / null_loglik)
  mcrsquared_values[[i]] <- data.frame(cbind(mcfadden_r2PCsEnv,mcfadden_r2PCsEnvBeta,mcfadden_r2Beta))
  
  #roc_dataPCs <- roc(na.omit(as.numeric(df$trait)), predict(fitlogPCs, type = "response"))
  #roc_dataPCsEnv <- roc(na.omit(as.numeric(df$trait)), predict(fitlogPCsEnv, type = "response"))
  #roc_dataPCsEnvBeta <- roc(na.omit(as.numeric(df$trait)), predict(fitlogPCsEnvBeta, type = "response"))
  #roc_dataPCsEnvSym <- roc(na.omit(as.numeric(df$trait)), predict(fitlogPCsEnvSym, type = "response"))
  #roc_datafull <- roc(na.omit(as.numeric(df$trait)), predict(fitlogfull, type = "response"))
  # 
  #roc_data_list[[i]] <- roc_data
  #rsquared_values_null[i] <- summary(fit_null)$r.squared
}
r2df=do.call("rbind",rsquared_values) 
r2df_melt <- melt(do.call("rbind",rsquared_values))
COR=cor.test(r2df$fitBeta,r2df$countLTER2)
COR

r2df=do.call("rbind",mcrsquared_values) 
r2df_melt <- melt(do.call("rbind",mcrsquared_values))

#auc_values[i] <- mean(sapply(roc_curves[[df_names[i]]], function(roc_data) auc(roc_data)))

ggplot(r2df_melt %>%
         filter(value > 0), aes(x=variable, y=value))+geom_boxplot()+labs(title = " beta cutoff p<0.1",
                                                                          y = "R2",x="Model")
#Pretty 
r2df$id= seq_len(nrow(r2df))
r2df_melt <- melt(r2df, id.vars = "id")

p1=ggplot(r2df_melt, aes(x=variable, y=value, fill=variable))+
  labs(y = expression("R"^2),x="")+ theme_classic()+theme(text=element_text(size = 15),axis.text.x = element_text(angle = 40, vjust =1, hjust=1))+
  geom_point(alpha = 0.2, color = "white") +  # Adjust alpha for transparency
  geom_line(aes(group = id), alpha = 0.3, color = "grey")+  # Adjust alpha for transparency
  geom_boxplot()+
  scale_x_discrete(limits = c("fitPCsEnv", "fitBeta", "fitPCsEnvBeta")
                                  , labels = c("Lat/Long+SA+Depth", "PGS Only", "Lat/Long+SA+Depth+PGS"))+
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"), guide="none")+  # Set your custom colors
  scale_y_continuous(breaks = seq(0, 1, by = .25))+ylim(0,1)+
  theme(aspect.ratio=2)
p1

result <- wilcox.test(r2df$fitPCsEnv, r2df$fitBeta)
result[["p.value"]]

result <- wilcox.test(r2df$fitBeta,r2df$fitPCsEnvBeta)
result[["p.value"]]









########################################Troubleshooting
##Option 1
file.names<-basename(list.files(path = directory_path,pattern = glob2rx("MayOct.plink.log.pgscomp_allcovs.*.RData"))) 
file.names

  pgscomp <- readRDS(file.path(directory_path, file.names[5]))
  df_name <- gsub("MayOct.plink.log.pgscomp_allcovs|\\.RData", "", file.names[i]) 
  lmrsquared_values <- numeric(length(pgscomp))
  for (j in 1:length(pgscomp)) {
    j=1
    df <- pgscomp[[j]]
    fitbetalm <- lm(as.numeric(df$trait) ~ as.numeric(df$sumbeta))
    lmrsquared_values[j] <- summary(fitbetalm)$r.squared
  }
  allpgscomps[[df_name]] <- lmrsquared_values
  #rsquared_values_null[i] <- summary(fit_null)$r.squared
}
means <- numeric(length(allpgscomps))
sds <- numeric(length(allpgscomps))
df_names <- names(allpgscomps)

for (i in 1:length(allpgscomps)) {
  values <- unlist(allpgscomps[[df_names[i]]])
  means[i] <- mean(values)
  sds[i] <- sd(values)
}

# Create a dataframe to store the results
summary_df <- data.frame(DataFrame = df_names, Mean = means, SD = sds)

#####vs 
##Option 2
#Taking p<0.01
pgscomp <- readRDS(file="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/MayOct.plink.log.pgscomp_allcovs.0.01.RData")
#merge with trait data to get other covariates to compare
covars=read.table(file="/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/MayOctplink_covars", header=T)
colnames(covars)[1]="Adapter"
merge_with_covars <- function(df) {
  merged_df <- merge(df, covars, by = "Adapter")
  return(merged_df)
}
# Use lapply to apply the merge_with_new function to each data frame in the list
merged_dataframes_list <- lapply(pgscomp, merge_with_covars)

rsquared_values=list()
for (i in 1:length(merged_dataframes_list)) {
  df <- merged_dataframes_list[[i]]
  fitBeta <- lm(as.numeric(df$trait) ~ as.numeric(df$sumbeta))
  summary(fitBeta)$r.squared
  newdf=data.frame(cbind(summary(fitBeta)$r.squared))
  print(summary(fitBeta)$r.squared)
  colnames(newdf)=c("fitBeta")
  rsquared_values[[i]] <- newdf
}
r2df=do.call("rbind",rsquared_values) 
r2df_melt <- melt(do.call("rbind",rsquared_values))

mean(r2df$fitBeta)


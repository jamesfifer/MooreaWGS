library(stringr)
library(ggplot2)
library(dplyr)
setwd("/projectnb/wgsgrp/CatFastq/All3Lanes/Analyses/Intermediate/")

file.names<-basename(list.files(pattern = glob2rx("*afreq"))) 
file.names
remove_after_third_period <- function(x) {
  str_replace(x, "(^([^.]*)\\.([^.]*)\\.([^.]*)).*", "\\1")
}
locs <- unique(sapply(file.names, remove_after_third_period))


freqchanges=list()
for(i in 1:length(locs)){
  #filename=paste0(locs,"")
  pre=read.table(file=paste0(locs[i],".May.afreq"))
  post=read.table(file=paste0(locs[i],".October.afreq"))
  FreqChange=post$V5-pre$V5
  freqchanges[[locs[i]]] =data.frame(FreqChange)
  dfname=paste0(locs[i],"FreqChange")
  colnames(freqchanges[[locs[i]]])=dfname
}

combined_df <- na.omit(do.call(cbind, freqchanges))
#checkign the two BR "major effect" loci 
combined_df$locus=pre$V2
huh=combined_df %>% filter(locus=="chr12_12179676")
#chr2_1350173
#chr12_12179676


#checking FR sites
FRsites=read.table(file="MayOct2FR.plink.all.PrePost.glm.logistic.0.1.sites")
FRsites$snp=paste(FRsites$V1,FRsites$V2,sep="_")
freqchanges=list()
for(i in 1:length(locs)){
  #filename=paste0(locs,"")
  pre=read.table(file=paste0(locs[i],".May.afreq"))
  post=read.table(file=paste0(locs[i],".October.afreq"))
  presub=subset(pre,V2 %in% FRsites$snp)
  postsub=subset(post,V2 %in% FRsites$snp)
  FreqChange=postsub$V5-presub$V5
  freqchanges[[locs[i]]] =data.frame(FreqChange)
  dfname=paste0(locs[i],"FreqChange")
  colnames(freqchanges[[locs[i]]])=dfname
}

combined_df <- na.omit(do.call(cbind, freqchanges))

#Mean frequency changes 
col_means <- colMeans(combined_df, na.rm = TRUE)
print(col_means)


#cor_LTER1v2=cov(combined_df$plink2.forereef.LTER1FreqChange,combined_df$plink2.forereef.LTER2FreqChange)/sqrt(var(combined_df$plink2.forereef.LTER1FreqChange)*var(combined_df$plink2.forereef.LTER2FreqChange))
#Pinsky paper does ^, but thats the same thing as cor...
#cor_LTER1v2=cor(combined_df$plink2.forereef.LTER1FreqChange,combined_df$plink2.forereef.LTER2FreqChange)

# Get the list of column names in combined_df
column_names <- colnames(combined_df)

# Generate all possible pairwise combinations
combinations <- combn(column_names, 2, simplify = FALSE)
# Compute correlations for all combinations
cor_results <- lapply(combinations, function(cols) {
  x <- combined_df[[cols[1]]]
  y <- combined_df[[cols[2]]]
  # Ensure no NA values for correlation calculation
  valid_indices <- complete.cases(x, y)
  x_valid <- x[valid_indices]
  y_valid <- y[valid_indices]
  
  if (length(x_valid) > 1) {
    cor(x_valid, y_valid)
  } else {
    NA
  }
})

# Combine results into a data frame
cor_df <- data.frame(
  Combination = sapply(combinations, paste, collapse = " vs "),
  Correlation = unlist(cor_results)
)


# Number of bootstrap iterations
n_boot <- 1000

# Create a list to store bootstrap results
boot_cor_results <- list()

# Perform bootstrapping for all pairwise combinations
for (comb in combinations) {
  x <- combined_df[[comb[1]]]
  y <- combined_df[[comb[2]]]
  
  # Ensure no NA values for correlation calculation
  valid_indices <- complete.cases(x, y)
  x_valid <- x[valid_indices]
  y_valid <- y[valid_indices]
  
  if (length(x_valid) > 1) {
    # Create a vector to store correlations for each bootstrap sample
    boot_cor <- numeric(n_boot)
    
    # Bootstrap loop
    for (i in 1:n_boot) {
      # Resample data with replacement
      boot_indices <- sample(seq_along(x_valid), replace = TRUE)
      x_boot <- x_valid[boot_indices]
      y_boot <- y_valid[boot_indices]
      
      # Compute correlation for bootstrap sample
      boot_cor[i] <- cor(x_boot, y_boot)
    }
    
    # Store results: mean and 95% confidence interval
    boot_cor_results[[paste(comb[1], "vs", comb[2])]] <- c(
      mean_cor = mean(boot_cor),
      ci_lower = quantile(boot_cor, 0.025),
      ci_upper = quantile(boot_cor, 0.975)
    )
  } else {
    # If there are not enough valid values to compute correlation
    boot_cor_results[[paste(comb[1], "vs", comb[2])]] <- c(NA, NA, NA)
  }
}

# Convert list to data frame
boot_cor_df <- do.call(rbind, boot_cor_results)
boot_cor_df <- as.data.frame(boot_cor_df)

# Add row names as a column for pairwise combination names
boot_cor_df$Combination <- rownames(boot_cor_df)
# Create a function to simplify the comparison strings
simplify_names <- function(comparison) {
  # Use regular expression to extract the 'environment.site' part from both terms in the comparison
  gsub("plink2\\.(.*?)FreqChange", "\\1", comparison)
}

# Apply the simplification to the Comparison column
boot_cor_df$Comparison <- sapply(boot_cor_df$Combination, simplify_names)

#saveRDS(boot_cor_df, file="FRboot_cor_df.RData")
#saveRDS(boot_cor_df, file="boot_cor_df.RData")
boot_cor_df=readRDS(file="boot_cor_df.RData")
boot_cor_df=readRDS(file="FRboot_cor_df.RData")

# Plot with error bars and median points
boot_cor_df <- boot_cor_df %>%
  mutate(SecondLabel = ifelse(grepl("forereef.*forereef", Comparison), "DFRvsDFR", 
                              ifelse(grepl("crest.*forereef", Comparison), "SFRvsDFR",
                                     ifelse(grepl("backreef.*forereef", Comparison), "BRvsDFR",
                                           "BRvsSFR"))))

boot_cor_df$Comparison_Label <- sapply(boot_cor_df$Comparison, function(x) {
  x <- gsub("backreef.", "", x)
  x <- gsub("forereef.", "", x)
  x <- gsub("crest.", "", x)
  return(x)
})


p=ggplot(boot_cor_df, aes(x = factor(Comparison, levels = unique(Comparison)))) +
  geom_linerange(aes(ymin = boot_cor_df$`ci_lower.2.5%`, ymax = boot_cor_df$`ci_upper.97.5%`), color = "black", linewidth = 0.75) +
  geom_point(aes(y = mean_cor), color = "black",fill="white", size = 2, pch=21) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.75) + # Add dotted line at y = 0
  scale_x_discrete(labels = boot_cor_df$Comparison_Label) +
  labs(title = "",
       x = "",
       y = "Convergent Correlation") +
  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p



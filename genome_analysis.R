## KETAKI HARDIKAR THESIS- s358950 ################################

## SCRIPT : Comparative genome analysis

#### 1. SETTING ENVIRONMENT ########################################################################

# Clear workspace
rm(list=ls())

# Close any open graphics devices
graphics.off()

# Working folder
folder <- "C:/Users/Dell/OneDrive/Desktop/Anopheles/malariagen"
setwd(folder)


#### 2. INSTALLING AND LOADING PACKAGES ############################################################

# Package names

packages <- c("ggplot2", "ggfortify" , "readr", "plotly", 
              "amap", "pbapply", "data.table", "devtools", "vhcub" )

# Install packages not yet installed

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading

invisible(lapply(packages, library, character.only = TRUE))

#### 3. SNP DATA ####################################################################################

## Read raw SNP data as data frame 

#(Raw data contains alleles per position for all samples,
# where 0 = allele absent, 1 = allele present)

raw_data<-as.data.frame(fread("matrix.tsv.gz"))
View(raw_data)
raw_data[1:5, 3081:3086]  

## Extract SNPs per position for all the samples from the raw data

sample_column_start<-2
sample_column_end<-3082

# Loop through each position with sapply

genotypes<-pbsapply(unique(raw_data$Pos),function(pos){
  
  # Perform another apply that returns a vector of 0, 1 or NA (one for each sample) per position
  
  apply(X = raw_data[raw_data$Pos==pos,sample_column_start:sample_column_end],MARGIN = 2,FUN = function(x){
    # x is equal to the 5 cells per sample per position representing 
    # REF, ALT, ALT, ALT, No call in that order
    # We find out which cell contains 1 (indicating allele present) and convert it into:
    # REF: 0
    # ALT: 1
    # No call: NA
    
    # first check there is only one cell==1 (1 in both the REFERENCE and ALT positions are set to NA)
    if (sum(x==1)==1) {
      cell_with_one<-which(x==1)
      if (cell_with_one==1){
        0
      } else if (cell_with_one==5){
        NA
      } else {
        1
      }
    } else {
      # return NA due to cells with more than one cell==1
      NA
    }
  })
},cl = 6)

nsamps<-nrow(genotypes)

#Filter data to only contain non-monomorphic SNPs

aaf<-apply(genotypes,2,sum,na.rm=T)
genotypes.filt<-genotypes[,which(aaf>0)]
dim(genotypes) #Check data dimensions
dim(genotypes.filt)


#### 4. PCA by classical MDS ####################################################################################

dists<-Dist(genotypes.filt,method="manhattan",nbproc=60)
pca<-cmdscale(dists)


#### METADATA ########

mt_meta <- read.csv(file = "meta.csv", sep=",", header=TRUE)
View(mt_meta)

#Match row names by sample id to colour the pca plot per taxon/region
meta<-mt_meta[match(rownames(pca),mt_meta$sample_id),] 
head(meta)

#Replace "taxon" with "region" to visualise respective distribution on PCA

species<-unique(meta$taxon)  
species.cols<-rainbow(length(species))
point.cols<-species.cols[match(meta$taxon,species)]

#Plot PCA

x11()

par(mar = c(5, 5, 5, 9), xpd = TRUE)

plot(pca[,1],pca[,2], col=point.cols, cex=0.5, 
     main="PCA by MDS (per taxon)", 
     xlab = "PC1", ylab="PC2", pch = 20)

legend("topright", cex=0.5, inset=c(-0.44,0),pch = 20, 
       col=unique(point.cols), legend = unique(meta$taxon),       
       border = "black")


## 5. Codon bias analysis #########################################################################################

#Codon bias analysis for dominant malaria vector mitogenomes by vhcub R package 

# Read DNA sequences from fasta file
fasta =  fasta.read("Dominant_vectors.fasta", "Dominant_vectors.fasta")
fasta.v <- fasta[[1]]


# RSCU analysis ###################################################
rscu <- RSCU.values(fasta.m)
View(rscu)
rscu_t <- t(rscu)

#Write RSCU values to csv

#Csv file with malaria vectors' names and ids extracted from fasta file sequences
malaria = read.csv(file = "Dominant_vectors.csv", sep=",", header=TRUE) 
mala_org = malaria[,2]
rscu_values =dplyr::as_tibble(rscu_t, rownames = "codons")
colnames(rscu_t) <-(mala_org)

write.csv(rscu_values, "RSCU_values.csv", row.names=FALSE)


#ENc-GC3 plot ################################################################

enc.df <- ENc.values(fasta.v)
gc.df <- GC.content(fasta.v)

ENc.GC3plot(enc.df, gc.df)

#GC content
gc <- GC.content(fasta.m)
View(gc)

#Write ENC and GC content values to csv

enc <- enc.df[2]
enc.gc3 = gc.df
enc.gc3 <- cbind(gc.df, new_col = enc) 

write.csv(enc.gc3, "ENC_GC3.csv", row.names=FALSE)



#############################################################################################################################################

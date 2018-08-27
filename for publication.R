# FINAL VERSION

source("R/read MaxQuant.R")
source("R/Digest_modified.R")
source("R/plot-scripts.R")

cv <- function(x){sd(x, na.rm = T)/mean(x, na.rm = T)}
cv.geometric <- function(x){sqrt((exp(1)^sd(x, na.rm = T)^2)-1)}

########################################################
### consider non-secreted proteins only based on hpa ###
########################################################

gene.info   <- read.table("data/gene-info_hpa18.xls", header = T, fill = T, sep = "\t")
genes.secr  <- gene.info[grep("Predicted secreted", gene.info$Protein.class),] # only considering non-secreted proteins

#############################################
### Read Protein-data from Maxquant input ###
#############################################

prot <- read.MaxQuant("data/gianluca/9 cell lines/proteinGroups.txt") # read MaxQuant output
prot <- prot[grep(";", prot$Gene.names, invert = T),] # continue with protein data mapping to one gene only
prot <- prot[grep("^REV|^CON", prot$Protein.IDs, invert = T),] # exclude reverse and contaminant sequences from analysis and normalization
temp <- read.table("data/translate-gene-names.tsv", sep = "\t", header = T) # translate from uniprot it ensembl
prot <- merge(prot, temp, by.x = "Majority.protein.IDs", by.y = "Majority.protein.IDs", all = F, all.x = T, all.y =F) # introduce HPA classifications
prot$Gene.names <- prot$Gene
prot <- prot[,1:189] # delete additional columns added by bind
prot <- prot[!is.na(prot$Gene.names) & !duplicated(prot$Gene.names),] # exclude genes without any gene-name from HPA, cannot make comparisons with RNA levels
row.names(prot) <- prot$Gene.names # add gene-names as row names

#############################
### median center dataset ###
#############################

repl.median         <- apply(prot[,79:108],2, function(x){median(x, na.rm = T)}) # median intensity of tmt-dataset per sample
prot.norm           <- prot
prot.norm[,79:108]  <- t(apply(prot[,79:108],1, function(x){x/repl.median})) # normalize to median intensity
prot.norm$Gene.names <- paste(prot.norm$Gene.names) 
prot.norm           <- prot.norm[!prot.norm$Gene.names %in% genes.secr$Gene,] # exclude genes predicted to be secreted according to HPA classification v18
norm.12             <- median(prot[,79], na.rm = T) / median(prot[,89], na.rm = T) # normalize channels, transform ch1 to ch2 (median intensity)
norm.32             <- median(prot[,99], na.rm = T) / median(prot[,89], na.rm = T) # normalize channels, transform ch3 to ch2

#############################
### RNA data from RNA-seq ###
#############################

exp.tpm <- read.csv("data/cheng/P7652_TPM_Kallisto.csv") # median numbers from Kallisto, including ENSG + genes
hpa.tpm <- read.table("data/tpm_cell-lines(56)_genes(19599).xls", sep = "\t", header = T) # Read RNA data
tpm.all <- read.table("data/cheng/TPM_for_all_cell_line_samples.txt", sep = "\t", header = T, fill = T) # rna-dataset
layout  <- read.csv("data/cheng/layout.csv") # layout for combining files with MaxQuant output format
tpm.all <- tpm.all[,c("X", as.matrix(layout$sample[layout$order]))] # re-order data-frame based on MaxQuant layout
tpm.all <- merge(tpm.all, exp.tpm[,c("EnsemblIDs", "Gene.name")], by.x = "X", by.y = "EnsemblIDs") # introduce gene names
tpm.all <- merge(tpm.all, hpa.tpm[,c("Gene", "Hep.G2..tpm.")],  by.x = "Gene.name", by.y ="Gene", all.x = F, all.y = F)
tpm.all[, c("Hep2G.P4", "Hep2G.P5", "Hep2G.P7")] <- tpm.all$Hep.G2..tpm. # add old hpa data instead of HepG2 that failed the RNAseq experiment
tpm.all <- tpm.all[, c(2,3:29,1)] # delete extra column introduced from above
tpm.all <- tpm.all[!duplicated(tpm.all$Gene.name),]
row.names(tpm.all) <- tpm.all$Gene.name
tpm.all <- fill.na(tpm.all, 2:28, cutoff = 1, fill = NA ) # set all values without RNA-values to NA
tpm.all <- tpm.all[,-1]
tpm.all <- tpm.all[,1:27]
colnames(tpm.all) <- c(paste0(cells, ".tpm1"), paste0(cells, ".tpm2"), paste0(cells, ".tpm3"))

######################################
### calculate geometric cv for RNA ###
######################################

df.tpm <- tpm.all[row.names(tpm.all) %in% row.names(prot.norm),]
cv.summary.tpm <- data.frame(matrix(nrow = nrow(df.tpm), ncol = 12, NA))
colnames(cv.summary.tpm) <- c(paste0(cells, ".tpm.cv"), "rep1.tpm.cv", "rep2.tpm.cv", "rep3.tpm.cv")
row.names(cv.summary.tpm) <- row.names(df.tpm)
for(i in 1:9){
  cv.summary.tpm[,i] <- apply(log2(df.tpm[,c(0+i,9+i, 18+i)]), 1, cv.geometric)
}
cv.summary.tpm$HEPG2.tpm.cv <- NA # three unique replicates
cv.summary.tpm[,10] <- apply(log2(df.tpm[,c(1:9)]), 1, cv.geometric)
cv.summary.tpm[,11] <- apply(log2(df.tpm[,c(10:18)]), 1, cv.geometric)
cv.summary.tpm[,12] <- apply(log2(df.tpm[,c(19:27)]), 1, cv.geometric)
cv.summary.tpm$bio <- apply(cv.summary.tpm[,10:12], 1, median, na.rm = T) 
cv.summary.tpm$tech <- apply(cv.summary.tpm[,1:9], 1, median, na.rm = T) 

#calculate geometric cv for proteins
df.prot <- data.frame(prot.norm[,80:88]/norm.12 ,prot.norm[,80:88+10],prot.norm[,80:88+20]/norm.32) # normalize each channel to represent ion intensites for rep2 (median)
df.prot <- fill.na(df.prot, cutoff = 0, fill = NA)
colnames(df.prot) <- c(paste(cells, "tmt1", sep = "."), paste(cells, "tmt2", sep = "."), paste(cells, "tmt3", sep = "."))
cv.summary.tmt <- data.frame(matrix(nrow = nrow(df.prot), ncol = 12, NA))
row.names(cv.summary.tmt) <- row.names(df.prot)
colnames(cv.summary.tmt) <- c(paste0(cells, ".tmt.cv"), "rep1.tmt.cv", "rep2.tmt.cv", "rep3.tmt.cv")

for(i in 1:9){
  cv.summary.tmt[,i] <- abs(apply(log(df.prot[,c(0+i,9+i, 18+i)]), 1, cv.geometric))
}
cv.summary.tmt[,10] <- abs(apply(log(df.prot[,c(1:9)]), 1, cv.geometric))
cv.summary.tmt[,11] <- abs(apply(log(df.prot[,c(10:18)]), 1, cv.geometric))
cv.summary.tmt[,12] <- abs(apply(log(df.prot[,c(19:27)]), 1, cv.geometric))

cv.summary.tmt$tmt.bio.cv <- apply(cv.summary.tmt[,10:12], 1, median, na.rm = T) 
cv.summary.tmt$tmt.tech.cv <- apply(cv.summary.tmt[,1:9], 1, median, na.rm = T)

############################################
### calculate gene-specific correlations ###
############################################

comb <- merge(df.prot, df.tpm, by = 0, all = F)
rownames(comb) <- comb$Row.names
comb <- comb[,-1]

comb.median <- data.frame(matrix(NA, nrow = nrow(comb), ncol = 18))
row.names(comb.median) <- row.names(comb)
for(i in 1:9){
  comb.median[,i] <- apply(comb[,c(i+0, i+9, i+18)], 1, median, na.rm = T)
}
for(i in 1:9){
  comb.median[,i+9] <- apply(comb[,c(i+27, i+36, i+45)], 1, median, na.rm = T)
}
colnames(comb.median) <- c(paste(cells, "tmt.median", sep = "."), paste(cells, "tpm.median", sep = "."))

#######################################################
### correlate median rna with median protein levels ###
#######################################################

correlate <- as.vector(which(!rowSums(is.na(comb.median[,1:9] / comb.median[,10:18])) > 4)) # must have measurements in 5 or more cell-lines

comb.median$cor.r <- NA
comb.median$cor.rho <- NA
comb.median$p.r <- NA
comb.median$p.rho <- NA

comb.median$log.cor.r <- NA
comb.median$log.cor.rho <- NA
comb.median$log.p.r <- NA
comb.median$log.p.rho <- NA
for (i in correlate){
  pearson <- cor.test(as.numeric(comb.median[i,1:9]),as.numeric(comb.median[i,10:18]), use = "pairwise.complete.obs", na.action = "na.exclude", method = "pearson")
  comb.median$cor.r[i] <- pearson$estimate
  comb.median$p.r[i] <- pearson$p.val
  spearman <- cor.test(as.numeric(comb.median[i,1:9]),as.numeric(comb.median[i,10:18]), use = "pairwise.complete.obs", na.action = "na.exclude", method = "spearman")
  comb.median$cor.rho[i] <- spearman$estimate
  comb.median$p.rho[i] <- spearman$p.val
  
  pearson.log <- cor.test(as.numeric(log2(comb.median[i,1:9])),as.numeric(log2(comb.median[i,10:18])), use = "pairwise.complete.obs", na.action = "na.exclude", method = "pearson")
  comb.median$log.cor.r[i] <- pearson.log$estimate
  comb.median$log.p.r[i] <- pearson.log$p.val
  spearman.log <- cor.test(as.numeric(log2(comb.median[i,1:9])),as.numeric(log2(comb.median[i,10:18])), use = "pairwise.complete.obs", na.action = "na.exclude", method = "spearman")
  comb.median$log.cor.rho[i] <- spearman.log$estimate
  comb.median$log.p.rho[i] <- spearman.log$p.val
}

comb.median$p.r.adj <- p.adjust(comb.median$p.r, method = "fdr")
comb.median$p.rho.adj <- p.adjust(comb.median$p.rho, method = "fdr")

comb.median$log.p.r.adj <- p.adjust(comb.median$log.p.r, method = "fdr")
comb.median$log.p.rho.adj <- p.adjust(comb.median$log.p.rho, method = "fdr")





comb.median$fc.prot <- apply(comb.median[,1:9], 1, fold)
comb.median$fc.rna <- apply(comb.median[,10:18], 1, fold)


comb.all <- merge(comb.median, df.prot, by.x = "row.names", by.y = "row.names")
row.names(comb.all) <- comb.all[,1]
comb.all <- comb.all[,-1]
comb.all <- merge(comb.all, df.tpm, by.x = "row.names", by.y = "row.names")
row.names(comb.all) <- comb.all[,1]
comb.all <- comb.all[,-1]
comb.all <- merge(comb.all, cv.summary.tmt, by.x = "row.names", by.y = "row.names")
row.names(comb.all) <- comb.all[,1]
comb.all <- comb.all[,-1]
comb.all <- merge(comb.all, cv.summary.tpm, by.x = "row.names", by.y = "row.names")



write.csv(comb.all, "rtp_data-output.csv", row.names = F)



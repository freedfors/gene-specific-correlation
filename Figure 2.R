
df <- read.csv("combined-data.csv") # all protein values in one data-frame
input$cutoff <- 1

df$include.tpm <- 1
df$include.tpm[which(df$bio.tpm.cv.median/df$tech.tpm.cv.median < input$cutoff)] <- 0
exclude.tpm <- df[which(df$include.tpm == 0), "Gene"] # list of genes to exclude based on tpm cutoff

df$include.tmt <- 1
df$include.tmt[which(df$bio.tmt.cv.median/df$tech.tmt.cv.median < input$cutoff)] <- 0
exclude.tmt <- df[which(df$include.tmt == 0), "Gene"] 

df$include <- 1
df$include[which(!df$include.tpm == 1)] <- 2
df$include[which(!df$include.tmt == 1)] <- 0

dsub <- df[!df$Gene %in% exclude.tmt,]
dsub <- dsub[!dsub$Gene %in% exclude.tpm,]


x.all  <- dsub$log.cor.r
x.all <- x.all[!is.na(x.all)]
x.rna  <- dsub$log.cor.r[dsub$fc.rna >= 5]

hist(x.all, breaks = 30, col  = "gray", border = 'white',
     xlab = "Pearson's r", main = paste0("Gene-specific correlation betweeen\nRNA and protein levels for ", length(x.all), " genes"),  ylim = c(0,600))
hist(x.rna, add = T, breaks = 30, col  = "#265CA8", border = 'white')
abline(h = 0, v = median(x.all, na.rm = T), col = "gray")
abline(h = 0, v = median(x.rna, na.rm = T), col = "#265CA8")


#DEDFE0
x <- x[which(x > input$cor[1] & x<= input$cor[2])]
x <- x[!is.na(x)]

hist(x, breaks = 40, col  = "#265CA8", border = 'white',
     xlab = "Pearson's r", main = paste0("Gene-specific correlation betweeen RNA and protein levels for ", length(x), " genes"))

# plot CDK2

dplot <- data.frame(t(df[1010, 2:19]), type = c(rep("prot", 9), rep("rna", 9)))
dplot$cell <- gsub("\\..*", "", row.names(dplot))
dplot$order <- 1:9
ggplot(dplot, aes(x=order, y = X1010, col = type)) + geom_bar(stat = "identity", position = "dodge")


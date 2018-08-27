library(Metrics)
library(ggplot2)
library(grid)
library(gridExtra)
library(ClassDiscovery)
library(gplots)
library(reshape2)
library(PerformanceAnalytics)
library(Rtsne)
library(scatterplot3d)
library(ggrepel)
library(openxlsx)
library(functional)
library(RColorBrewer)
library(mvtboost)


# useful functions ####
cv <- function(x){sd(x, na.rm = T)/mean(x, na.rm = T)}

fill.na <- function(dataframe, columns, cutoff = 0, fill = NA){
  subset <- dataframe[,columns]
  subset[subset <= cutoff] <- fill
  dataframe[,columns] <- subset
  return(dataframe)
}

col.9 <- c("#337FC1", "#70C6E7", "#9FC554",
           "#6D579E", "#DA368D", "#E47639",
           "#327567", "#4EB2C4", "#A3CEA3")

cell.col <- c(col.9, addalpha(col.9, 0.8), addalpha(col.9, 0.6))

# add alpha
addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

# colorRampPaletteAlpha()
colorRampPaletteAlpha <- function(colors, n=32, interpolate='linear') {
  # Create the color ramp normally
  cr <- colorRampPalette(colors, interpolate=interpolate)(n)
  # Find the alpha channel
  a <- col2rgb(colors, alpha=T)[4,]
  # Interpolate
  if (interpolate=='linear') {
    l <- approx(a, n=n)
  } else {
    l <- spline(a, n=n)
  }
  l$y[l$y > 255] <- 255 # Clamp if spline is > 255
  cr <- addalpha(cr, l$y/255.0)
  return(cr)
}

# calculate fold change
fold <- function(x){max(x, na.rm = T)/min(x, na.rm = T)}


# order samples appear in result files (TMT and RNAseq), ie sample order from MaxQuant output (reporter 1-9)
cells   <- c("SKMEL24",
             "SKMEL28", 
             "U2OS", 
             "A549", 
             "U251", 
             "BJ", 
             "RT4", 
             "HEPG2", 
             "CACO2") 

theme <- theme_minimal() + theme(panel.grid = element_blank(),
                                 panel.background = element_rect(fill = "white", colour = "grey50"),
                                 legend.title= element_blank())


# create folder for today ###
today <- format(Sys.time(), "%Y%m%d")
path <- paste0("output/", today, "/")
dir.create(path) # create path for output files
dir.create(paste0(path, "/tables")) # create path for output files

# define custom palettes
my_palette_rna <- colorRampPalette(c("#8D4A99", "#F2F0F7"))(n = 32) # RNA heatmap
my_palette_prot <- colorRampPalette(c("#137EC2", "#F2F0F7"))(n = 32) # RNA heatmap



# NEW DATA FROM HERE

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


dsub <- dsub[!is.na(dsub$log.cor.r),] # Should be same number of genes gues

dsub$cv <- apply(dsub[,11:19], 1, cv)

most.var <- dsub[order(dsub$bio.tpm.cv.median, decreasing = T),][1:ceiling(nrow(dsub)*0.1),]


# Select input
protein <- most.var[2:10]
rna <- most.var[11:19]
row.names(protein) <- most.var$Gene
stop()

# all proteins
protein <- dsub[,2:10]
rna <- dsub[,11:19]
row.names(protein) <- dsub$Gene
stop()


#
rep = "motvar"
circular = F

df <- data.frame(protein, rna)
row.names(df) <- row.names(protein)
df$gene <- rownames(df)
df <- data.frame(df, df[,1:9]/df[,10:18]) # add rtp
colnames(df) <- c(colnames(df)[1:18], "gene", paste0(cells, ".rtp"))
  
empty <-matrix(NA, nrow = 9, ncol = 9)
  
  
rna.cor <- empty
ptr.cor <- empty
  prmed.cor <- empty
  rna.pval <- empty
  met <- "pearson"
  error.free <- NULL
  error.ptr <- NULL
  rmse.all.ptr <- NULL
  rmse.all.free <- NULL
  tab.ptr.error <- data.frame(matrix(NA, ncol = 9, nrow = nrow(protein)))
  tab.free.error <- data.frame(matrix(NA, ncol = 9, nrow = nrow(protein)))
  row.names(tab.ptr.error) <- row.names(df)
  row.names(tab.free.error) <- row.names(df)
  p1 <- NULL
  p2 <- NULL
  p3 <- NULL
  for(i in 1:9){
    if(circular == F){
      cols <- 19+(1:9)[!1:9 %in% i]
      df$rtp.median <- apply(df[,cols], 1, median, na.rm = T) # här kan vi lägga till en uteslutning av respektive cell-linje och anropa det när vi räknar ut median rtp  
      df$prot.median <- apply(df[,cols-19], 1, median, na.rm = T)
    }else{
      df$rtp.median <- apply(df[,20:(20+8)], 1, median, na.rm = T) # här kan vi lägga till en uteslutning av respektive cell-linje och anropa det när vi räknar ut median rtp
      df$prot.median <- apply(df[,1:9], 1, median, na.rm = T)
    }
  
    for(j in 1:9){
      rna.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      ptr.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]*df$rtp.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      prmed.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df$prot.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw for protein median prediction
      dplot.1 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1))]
      dplot.2 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 29)]
      dplot.2[,2] <- dplot.2[,2]*dplot.2[,3] # multiply with PTR
      dplot.3 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 30)] # only median protein
      xlab = colnames(df)[(1 + (i-1) * 1)]
      ylab = colnames(df)[(10 + (j-1) * 1)]
      colnames(dplot.1) <- c("x", "y")
      colnames(dplot.2) <- c("x", "y")
      colnames(dplot.3) <- c("x", "y", "z")
      
      dplot.2$error.log2 <- log2(dplot.2$y/dplot.2$x) # error ptr
      dplot.3$error.log2 <- log2(dplot.3$z/dplot.3$x) # error mrna-free
      
      
      # plot here
      
      if(i == j){
        alpha.v <- 0.35
        p1[[i]] <- ggplot(dplot.1, aes(x = log2(x), y = log2(y))) + geom_point(aes(colour=  log10(y)), alpha = alpha.v) + xlab(expression(log[2] *"(protein)")) + ylab(expression(log[2] *"(TPM)")) + theme + scale_color_distiller(palette = "RdPu") + theme(legend.position = "none") + ylim(2,12) + xlim(-10,10) +
          annotate("text", x = 7, y = 3, label = paste0(round(rna.cor[i,j], 2))) + ggtitle(cells[i])
        p2[[i]] <- ggplot(dplot.2, aes(x = log2(x), y = log2(y))) + geom_point(colour = col.9[i] , alpha = alpha.v) + xlab(expression(log[2] *"(protein)")) + ylab(expression(log[2] *"(PTR-based pred.)")) + theme + theme(legend.position = "none") + ylim(-10, 10) + xlim(-10,10) +
          annotate("text", x = 7, y = -8.4, label = paste0(round(ptr.cor[i,j], 2))) + geom_abline(intercept = 0, slope = 1, lty = 2)
        p3[[i]] <- ggplot(dplot.3, aes(x = log2(x), y = log2(z))) + geom_point(colour = col.9[i] , alpha = alpha.v) + xlab(expression(log[2] *"(protein)")) + ylab(expression(log[2] *"(RNA-free pred.)")) + theme +  theme(legend.position = "none") + ylim(-10,10) + xlim(-10,10) +
          annotate("text", x = 7, y = -8.4, label = paste0(round(prmed.cor[i,j], 2))) + geom_abline(intercept = 0, slope = 1, lty = 2)
        
        pdf(paste0(path, "error-regression_", cells[i], "_", rep, ".pdf"), height = 3, width = 3)
        hist(dplot.2$error.log2, breaks = seq(-100000, 100000, by = 0.5), xlim = c(-5,5), col = rgb(0.25,0.71,0.76,0.2), main = paste0(cells[i], " - log2 residual error"), freq = F)
        abline(v = median(dplot.2$error.log2, na.rm = T),col="cyan")
        hist(dplot.3$error.log2,breaks = seq(-100000, 100000, by = 0.5), xlim = c(-5,5), add = T, col = rgb(0.172,0.1,0.760,0.2), freq = F)
        abline(v = median(dplot.3$error.log2, na.rm = T),col="blue")
        dev.off()
        error.ptr <- c(error.ptr, dplot.2$error.log2)
        error.free <- c(error.free, dplot.3$error.log2)
        
        rmse.error.ptr <- dplot.2[!rowSums(is.na(dplot.2[,c("x", "y")])) > 0, c("x", "y")]
        rmse.error.free <- dplot.3[!rowSums(is.na(dplot.3[,c("x", "z")])) > 0, c("x", "z")]
        
        print(paste0(cells[i], ", RMSE, ptr=", 
                     round(rmse(log2(rmse.error.ptr[,1]), log2(rmse.error.ptr[,2])), 2), 
                     ", free=",
                     round(rmse(log2(rmse.error.free[,1]), log2(rmse.error.free[,2])), 2),
                     "; MAE, ptr = ",
                     round(mae(log2(rmse.error.ptr[,1]), log2(rmse.error.ptr[,2])), 2), 
                     ", free=",
                     round(mae(log2(rmse.error.free[,1]), log2(rmse.error.free[,2])), 2)))
        
        rmse.all.ptr <- rbind(rmse.all.ptr, rmse.error.ptr)
        rmse.all.free <- rbind(rmse.all.free, rmse.error.free)
      }
    }
  }
  error.ptr <- error.ptr[!is.na(error.ptr)]
  error.free <- error.free[!is.na(error.free)]
  table(abs(error.ptr) < 1) / length(error.ptr) * 100
  table(abs(error.free) < 1)/ length(error.free) * 100
  
  
  
  
  write.csv(prmed.cor, "temp.csv")
  pdf(paste0("prediction-mostvar.pdf"), height = 4, width = 4*9)
  do.call("grid.arrange", c(p1[1:9], nrow= 1))
  do.call("grid.arrange", c(p2[1:9], nrow= 1))
  do.call("grid.arrange", c(p3[1:9], nrow= 1))
  dev.off()
  
  # for caluclation only
  
  print(paste0("Total, RMSE, ptr=", 
               round(rmse(log2(rmse.all.ptr[,1]), log2(rmse.all.ptr[,2])), 2), 
               ", free=",
               round(rmse(log2(rmse.all.free[,1]), log2(rmse.all.free[,2])), 2),
               "; MAE, ptr = ",
               round(mae(log2(rmse.all.ptr[,1]), log2(rmse.all.ptr[,2])), 2), 
               ", free=",
               round(mae(log2(rmse.all.free[,1]), log2(rmse.all.free[,2])), 2)))
  



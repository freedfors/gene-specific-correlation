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

#add alpha
# addalpha()
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



# FUNKTION SOM PLOTTAR
rtp.calc.BA <- function(protein, rna, rep = "rep3", circular = F){ # for 9 cell-lines
  df <- data.frame(protein, rna)
  df$gene <- rownames(df)
  df <- data.frame(df, df[,1:9]/df[,10:18]) # add rtp
  colnames(df) <- c(colnames(df)[1:18], "gene", paste0(cells, ".rtp"))
  if(circular == F){
    cols <- 19+(1:9)[!1:9 %in% i]
    df$rtp.median <- apply(df[,cols], 1, median, na.rm = T) # här kan vi lägga till en uteslutning av respektive cell-linje och anropa det när vi räknar ut median rtp  
    df$prot.median <- apply(df[,cols-19], 1, median, na.rm = T)
  }else{
    df$rtp.median <- apply(df[,20:(20+8)], 1, median, na.rm = T) # här kan vi lägga till en uteslutning av respektive cell-linje och anropa det när vi räknar ut median rtp
    df$prot.median <- apply(df[,1:9], 1, median, na.rm = T)
  }
  write.csv(df, paste0(path, "tables/", rep, "labelfree.csv"))
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
  for(i in 1:9){
    for(j in 1:9){
      rna.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      ptr.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]*df$rtp.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      prmed.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df$prot.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw for protein median prediction
      #rna.pval[i,j] <- cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]))$p.val # rna raw
      #ptr.cor[i,j] <-  cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]*df$rtp.median))$estimate # rna raw
      #ptr.pval[i,j] <- cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]*df$rtp.median))$p.val # rna raw
      dplot.1 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1))]
      dplot.2 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 29)]
      dplot.2[,2] <- dplot.2[,2]*dplot.2[,3] # multiply with PTR
      dplot.3 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 30)] # only median protein
      xlab = colnames(df)[(1 + (i-1) * 1)]
      ylab = colnames(df)[(10 + (j-1) * 1)]
      colnames(dplot.1) <- c("x", "y")
      colnames(dplot.2) <- c("x", "y")
      colnames(dplot.3) <- c("x", "y", "z")
      
      dplot.2$error.log2 <- log2(dplot.2$y/dplot.2$x)
      dplot.3$error.log2 <- log2(dplot.3$z/dplot.3$x)
      if(i == j){
        alpha.v <- 0.3
        p1[[i]] <- ggplot(dplot.1, aes(x = log2(x), y = log2(y))) + geom_point(aes(colour=  log10(y)), alpha = alpha.v) + xlab(expression(log[2] *"(protein)")) + ylab(expression(log[2] *"(TPM)")) + theme + scale_color_distiller(palette = "RdPu") + theme(legend.position = "none") + ylim(2,12) + xlim(-10,10) +
          annotate("text", x = 7, y = 3, label = paste0(round(rna.cor[i,j], 2))) + ggtitle(cells[i])
        p2[[i]] <- ggplot(dplot.2, aes(x = log2(x), y = log2(y))) + geom_point(aes(colour=  log10(y)), alpha = alpha.v) + xlab(expression(log[2] *"(protein)")) + ylab(expression(log[2] *"(PTR-based pred.)")) + theme + scale_color_distiller(palette = "YlGnBu") + theme(legend.position = "none") + ylim(-10, 10) + xlim(-10,10) +
          annotate("text", x = 7, y = -8.4, label = paste0(round(ptr.cor[i,j], 2))) + geom_abline(intercept = 0, slope = 1, lty = 2)
        p3[[i]] <- ggplot(dplot.3, aes(x = log2(x), y = log2(z))) + geom_point(aes(colour=  log10(y)), alpha = alpha.v) + xlab(expression(log[2] *"(protein)")) + ylab(expression(log[2] *"(RNA-free pred.)")) + theme + scale_color_distiller(palette = "PuBu") + theme(legend.position = "none") + ylim(-10,10) + xlim(-10,10) +
          annotate("text", x = 7, y = -8.4, label = paste0(round(prmed.cor[i,j], 2))) + geom_abline(intercept = 0, slope = 1, lty = 2)
        pdf(paste0(path, "error-regression_", cells[i], "_", rep, ".pdf"), height = 3, width = 3)
        hist(dplot.2$error, breaks = seq(-100000, 100000, by = 0.5), xlim = c(-5,5), col = rgb(0.25,0.71,0.76,0.2), main = paste0(cells[i], " - residual error (%)"), freq = F)
        abline(v = median(dplot.2$error, na.rm = T),col="cyan")
        hist(dplot.3$error,breaks = seq(-100000, 100000, by = 0.5), xlim = c(-5,5), add = T, col = rgb(0.172,0.1,0.760,0.2), freq = F)
        abline(v = median(dplot.3$error, na.rm = T),col="blue")
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
  write.csv(prmed.cor, "temp.csv")
  pdf(paste0(path, rep, "_prediction.pdf"), height = 2, width = 2*9)
  do.call("grid.arrange", c(p1[1:9], nrow= 1))
  do.call("grid.arrange", c(p2[1:9], nrow= 1))
  do.call("grid.arrange", c(p3[1:9], nrow= 1))
  dev.off()
  
  # for caluclation only
  # length(error.ptr[which(error.ptr < 1 & error.ptr > -1)]) / length(which(!is.na(error.ptr)))
  # length(error.free[which(error.free < 1 & error.free > -1)]) / length(which(!is.na(error.free)))
  
  print(paste0("Total, RMSE, ptr=", 
               round(rmse(log2(rmse.all.ptr[,1]), log2(rmse.all.ptr[,2])), 2), 
               ", free=",
               round(rmse(log2(rmse.all.free[,1]), log2(rmse.all.free[,2])), 2),
               "; MAE, ptr = ",
               round(mae(log2(rmse.all.ptr[,1]), log2(rmse.all.ptr[,2])), 2), 
               ", free=",
               round(mae(log2(rmse.all.free[,1]), log2(rmse.all.free[,2])), 2)))
  
  pdf(paste0(path, "error-regression-all_", rep, ".pdf"), height = 4, width = 4)
  hist(error.ptr, breaks = seq(-100000, 100000, by = .5), xlim = c(-5,5), add = F, col = rgb(0.411,0.741,0.27, 0.2), freq = F)
  hist(error.free, breaks = seq(-100000, 100000, by = .5), xlim = c(-5,5), add = T, col = rgb(0.341,0.549,0.682,0.2), freq = F)
  lines(density(error.ptr[!is.na(error.ptr)], adjust = 2), col = rgb(0.411,0.741,0.27,1), lwd = 3)
  lines(density(error.free[!is.na(error.free)], adjust = 2), col = rgb(0.341,0.549,0.682,1), lwd = 3)
  abline(v = median(error.free, na.rm = T),col="blue")
  abline(v = -1 ,col="gray", lwd = 2, lty = 2)
  abline(v = 1 ,col="gray", lwd = 2, lty = 2)
  dev.off()
  
  pdf(paste(path, rep, "_rna-prediction_lf.pdf"), height = 5, width = 5)
  heatmap.2(rna.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  heatmap.2(ptr.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  heatmap.2(prmed.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  dev.off()
  
}


rtp.calc.abs <- function(protein, rna, rep = "rep3", circular = F){ # for 9 cell-lines
  df <- data.frame(protein, rna)
  df$gene <- rownames(df)
  df <- data.frame(df, df[,1:9]/df[,10:18]) # add rtp
  colnames(df) <- c(colnames(df)[1:18], "gene", paste0(cells, ".rtp"))
  if(circular == F){
    cols <- 19+(1:9)[!1:9 %in% i]
    df$rtp.median <- apply(df[,cols], 1, median, na.rm = T) # här kan vi lägga till en uteslutning av respektive cell-linje och anropa det när vi räknar ut median rtp  
    df$prot.median <- apply(df[,cols-19], 1, median, na.rm = T)
  }else{
    df$rtp.median <- apply(df[,20:(20+8)], 1, median, na.rm = T) # här kan vi lägga till en uteslutning av respektive cell-linje och anropa det när vi räknar ut median rtp
    df$prot.median <- apply(df[,1:9], 1, median, na.rm = T)
  }
  write.csv(df, paste0(path, "tables/", rep, "abs.csv"))
  empty <-matrix(NA, nrow = 9, ncol = 9)
  
  
  rna.cor <- empty
  ptr.cor <- empty
  prmed.cor <- empty
  met <- "pearson"
  error.free <- NULL
  error.ptr <- NULL
  
  
  for(i in 1:9){
    for(j in 1:9){
      rna.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      ptr.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]*df$rtp.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      prmed.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df$prot.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw for protein median prediction
      #rna.pval[i,j] <- cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]))$p.val # rna raw
      #ptr.cor[i,j] <-  cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]*df$rtp.median))$estimate # rna raw
      #ptr.pval[i,j] <- cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]*df$rtp.median))$p.val # rna raw
      dplot.1 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1))]
      dplot.2 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 29)]
      dplot.2[,2] <- dplot.2[,2]*dplot.2[,3] # multiply with PTR
      dplot.3 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 30)] # only median protein
      xlab = colnames(df)[(1 + (i-1) * 1)]
      ylab = colnames(df)[(10 + (j-1) * 1)]
      colnames(dplot.1) <- c("x", "y")
      colnames(dplot.2) <- c("x", "y")
      colnames(dplot.3) <- c("x", "y", "z")
      
      dplot.2$error.log2 <- log2(dplot.2$y/dplot.2$x)
      dplot.3$error.log2 <- log2(dplot.3$z/dplot.3$x)
      if(i == j){
        alpha.v <- 0.3
        p1[[i]] <- ggplot(dplot.1, aes(x = log10(x), y = log10(y))) + geom_point(aes(colour=  log10(y)), alpha = 0.2) + xlab(expression(log[] *"(protein)")) + ylab(expression(log[10] *"(TPM)")) + theme + scale_color_distiller(palette = "RdPu") + theme(legend.position = "none") + ylim(0,5) + xlim(2,8) +
          annotate("text", x = 20, y = 0, label = paste0("rho=", round(rna.cor[i,j], 2))) + geom_abline(intercept = 0, slope = 1, lty = 2)
        p2[[i]] <- ggplot(dplot.2, aes(x = log10(x), y = log10(y))) + geom_point(aes(colour=  log10(y)), alpha = 0.2) + xlab(expression(log[10] *"(protein copies)")) + ylab(expression(log[10] *"(RNA*PTR)")) + theme + scale_color_distiller(palette = "YlGnBu") + theme(legend.position = "none") + ylim(2.5,7.5) + xlim(2.5,7.5) +
          annotate("text", x = 20, y = 6, label = paste0("rho=", round(ptr.cor[i,j], 2))) + geom_abline(intercept = 0, slope = 1, lty = 2)
        p3[[i]] <- ggplot(dplot.3, aes(x = log10(x), y = log10(z))) + geom_point(aes(colour=  log10(y)), alpha = 0.2) + xlab(expression(log[10] *"(protein copies)")) + ylab(expression(log[10] *"(median, protein)")) + theme + scale_color_distiller(palette = "PuBu") + theme(legend.position = "none") + ylim(2.5,7.5) + xlim(2.5,7.5) +
          annotate("text", x = 20, y = 6, label = paste0("rho=", round(prmed.cor[i,j], 2))) + geom_abline(intercept = 0, slope = 1, lty = 2)
        pdf(paste0(path, "error-regression_", cells[i], "_", rep, ".pdf"), height = 3, width = 3)
        hist(dplot.2$error, breaks = seq(-100000, 100000, by = 0.5), xlim = c(-5,5), col = rgb(0.25,0.71,0.76,0.2), main = paste0(cells[i], " - residual error (%)"), freq = F)
        abline(v = median(dplot.2$error, na.rm = T),col="cyan")
        hist(dplot.3$error,breaks = seq(-100000, 100000, by = 0.5), xlim = c(-5,5), add = T, col = rgb(0.172,0.1,0.760,0.2), freq = F)
        abline(v = median(dplot.3$error, na.rm = T),col="blue")
        dev.off()
        error.ptr <- c(error.ptr, dplot.2$error.log2)
        error.free <- c(error.free, dplot.3$error.log2)
      }
    }
  }
  
  
  pdf(paste0(path, rep, "_prediction.pdf"), height = 2, width = 2*9)
  do.call("grid.arrange", c(p1[1:9], nrow= 1))
  do.call("grid.arrange", c(p2[1:9], nrow= 1))
  do.call("grid.arrange", c(p3[1:9], nrow= 1))
  dev.off()
  
  pdf(paste0(path, "error-regression-all_", rep, ".pdf"), height = 4, width = 4)
  hist(error.ptr, breaks = seq(-100000, 100000, by = .5), xlim = c(-5,5), col = rgb(0.25,0.71,0.76,0.2), main = paste0("All residual error (%)"), freq = F)
  abline(v = median(error.ptr, na.rm = T),col="cyan")
  hist(error.free, breaks = seq(-100000, 100000, by = .5), xlim = c(-5,5), add = T, col = rgb(0.172,0.1,0.760,0.2), freq = F)
  abline(v = median(error.free, na.rm = T),col="blue")
  abline(v = -1 ,col="red")
  abline(v = 1 ,col="red")
  dev.off()
  
  pdf(paste(path, rep, "_rna-prediction_lf.pdf"), height = 5, width = 5)
  heatmap.2(rna.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  heatmap.2(ptr.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  heatmap.2(prmed.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  dev.off()
  
}






BA.density.plot <- function(BA.x, BA.y, group.vector,BA.y.name = "y", BA.x.name = "x", x.axis.log = F, linreg=F, set.x.lim = NA, return.results = F){
  
  
  #Test for normal distribution in differences:
  #Kolmogorov-Smirnov test
  ks <- ks.test(BA.y, "pnorm", mean(BA.y), sd(BA.y))
  warning("p-value has been inactivated by Fredrik")
  # if (ks$p.value<0.05){
  print("Kolmogorov-Smirnov normality test failed (p<0.05). Mean and 2.5th and 97.5th percentiles will be used as limits of agreement.")
  limit.mean <- median(BA.y)
  limit.upper <- quantile(BA.y, .975)
  limit.upper.label <- "97.5th percentile"
  limit.lower <- quantile(BA.y, 0.025)
  limit.lower.label <- "2.5th percentile"
  # } else{
  #   print(paste0("Kolmogorov-Smirnov normality p value: ",ks$p.value,">0.05. Median and +/- 1.96 sd will be used as limits of agreement."))
  #   limit.mean <- mean(BA.y)
  #   limit.upper <- limit.mean + 1.96*sd(BA.y)
  #   limit.upper.label <- "+1.96 SD"
  #   limit.lower <- limit.mean - 1.96*sd(BA.y)
  #   limit.lower.label <- "-1.96 SD"
  # }
  # 
  DF <- data.frame(BA.x=BA.x, BA.y=BA.y, group.vector=group.vector, 
                   outside.limits=BA.y<limit.lower | BA.y>limit.upper)
  
  
  limit.df <- data.frame(limit = c(limit.upper.label, "Mean", limit.lower.label),
                         x = c(0.05,0.05,0.05),
                         y = c(limit.upper, limit.mean, limit.lower),
                         percent.points.outside.within.label = c(100*length(BA.y[BA.y>limit.upper])/length(BA.y),
                                                                 100*length(BA.y[BA.y<limit.upper & BA.y>limit.lower])/length(BA.y),
                                                                 100*length(BA.y[BA.y<limit.lower])/length(BA.y)))
  
  
  theme.BA <- function(...) theme(legend.position="none",
                                  axis.text.y=element_blank(),
                                  axis.ticks.y=element_blank(),
                                  axis.title.y=element_blank(),
                                  axis.text.x=element_blank(),
                                  axis.ticks.x=element_blank(),
                                  axis.title.x=element_blank(),
                                  panel.background = element_blank(),
                                  panel.spacing = unit(0,"null"),
                                  plot.margin = unit(c(0,0,0,0),"lines"),...)
  
  plot.y.axis <- ggplot()+
    ylim(c(2*limit.lower, 2*limit.upper))+
    ylab(BA.y.name)+
    theme(axis.line.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.line.y = element_line(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = element_blank(),
          panel.border = element_blank())
  
  plot.x.axis <- ggplot()+
    xlim(c(min(BA.x), max(BA.x)))+
    xlab(BA.x.name)+
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          axis.line.x = element_line(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = element_blank(),
          panel.border = element_blank())
  
  if(x.axis.log) plot.x.axis <- plot.x.axis + scale_x_log10(limits = c(min(BA.x), max(BA.x)))
  
  plot.top <- ggplot(DF, aes(BA.x,fill=group.vector))+
    geom_density(alpha=0.2)+
    xlim(c(min(BA.x), max(BA.x)))+
    theme.BA(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border=element_blank(),
      axis.line.y=element_blank(),
      axis.line.x=element_blank())
  
  if(x.axis.log) plot.top <- plot.top + scale_x_log10(limits = c(min(BA.x), max(BA.x)))
  
  plot.empty <- ggplot()+
    geom_point(aes(1,1), colour="white")+
    theme(legend.position = "none",
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0,"null"),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.length = unit(0,"null"))
  
  plot.right <- ggplot(DF, aes(BA.y,fill=group.vector))+
    geom_density(alpha=0.2)+
    geom_vline(xintercept = c(limit.upper, limit.mean, limit.lower),
               linetype = "dashed",
               color = c("gray", "black", "gray"))+
    xlim(c(limit.lower-abs(limit.lower), limit.upper+abs(limit.upper)))+
    coord_flip()+
    theme.BA(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border=element_blank(),
      axis.line.y=element_blank(),
      axis.line.x=element_blank())
  
  
  
  plot.BA <- ggplot(DF, aes(BA.x,BA.y,color=group.vector))+
    geom_point(alpha=0.05, size=0, stroke=2)+
    geom_point(data=DF[DF$outside.limits,],
               fill = DF[DF$outside.limits,]$group.vector, 
               alpha=0.2, size=0, stroke=2) +
    geom_hline(yintercept = c(0, limit.upper, limit.mean, limit.lower),
               linetype = c("dotted", "dashed","dashed","dashed"),
               color = c("lightgray","gray", "black", "gray"))
  
  if(linreg){
    plot.BA <- plot.BA +
      geom_smooth(method=lm)
    
    #Linear regression of values:
    linreg.parameters <- lm(formula = BA.y~BA.x, data = DF)
  }
  
  plot.BA <- plot.BA +
    annotate("text",
             label = paste0(limit.df$limit, ": ", round(limit.df$y, digits = 3)),
             x = limit.df$x , 
             y = limit.df$y,
             vjust=-0.5,
             color = c("darkgray", "black", "darkgray")) + 
    annotate("text",
             label = paste0(round(limit.df$percent.points.outside.within.label, digits = 2), " %"),
             x = limit.df$x , 
             y = limit.df$y,
             vjust=1.2,
             color = c("darkgray", "black", "darkgray")) + 
    ylim(c(limit.lower-abs(limit.lower), limit.upper+abs(limit.upper)))+
    xlim(c(min(BA.x), max(BA.x)))+
    ylab(BA.y.name)+
    theme.BA()
  
  if(x.axis.log) plot.BA <- plot.BA + scale_x_log10(limits = c(min(BA.x), max(BA.x)))
  
  if(!is.na(set.x.lim)){
    plot.BA <- plot.BA + xlim(set.x.lim)
    plot.top <- plot.top + xlim(set.x.lim)
    plot.x.axis<- plot.x.axis + xlim(set.x.lim)
  }
  
  
  grid.arrange(plot.empty, plot.top, plot.empty, 
               plot.y.axis, plot.BA, plot.right, 
               plot.empty,plot.x.axis,plot.empty,
               ncol=3, nrow=3, widths=c(0.5, 4, 1), heights=c(1, 4, 0.6))
  if(return.results == T){
    return(data.frame(limit.mean = median(BA.y),
                      limit.upper = quantile(BA.y, .975),
                      limit.lower <- quantile(BA.y, 0.025)))
    
  }
  
}

# compare all targets that don't pass BA-analysis compared to 



# function
rtp.calc <- function(protein, rna, rep = "rep3"){ # for 9 cell-lines
  df <- data.frame(protein, rna)
  df$gene <- rownames(df)
  df <- data.frame(df, df[,1:9]/df[,10:18])
  colnames(df) <- c(colnames(df)[1:18], "gene", paste0(cells, ".rtp"))
  df$rtp.median <- apply(df[,20:(20+8)], 1, median, na.rm = T)
  df$prot.median <- apply(df[,1:9], 1, median, na.rm = T)
  write.csv(df, paste0(path, "tables/", rep, "labelfree.csv"))
  empty <-matrix(NA, nrow = 9, ncol = 9)
  
  
  rna.cor <- empty
  ptr.cor <- empty
  prmed.cor <- empty
  met <- "spearman"
  for(i in 1:9){
    for(j in 1:9){
      rna.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      ptr.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]*df$rtp.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      prmed.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df$prot.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw for protein median prediction
      #rna.pval[i,j] <- cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]))$p.val # rna raw
      #ptr.cor[i,j] <-  cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]*df$rtp.median))$estimate # rna raw
      #ptr.pval[i,j] <- cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]*df$rtp.median))$p.val # rna raw
      dplot.1 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1))]
      dplot.2 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 29)]
      dplot.2[,2] <- dplot.2[,2]*dplot.2[,3]
      dplot.3 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 30)]
      xlab = colnames(df)[(1 + (i-1) * 1)]
      ylab = colnames(df)[(10 + (j-1) * 1)]
      colnames(dplot.1) <- c("x", "y")
      colnames(dplot.2) <- c("x", "y")
      colnames(dplot.3) <- c("x", "y", "z")
      if(i == j){
        p1[[i]] <- ggplot(dplot.1, aes(x = log2(x), y = log2(y))) + geom_point(aes(colour=  log10(y)), alpha = 0.2) + xlab(paste0("log2(protein) , ", xlab)) + ylab(paste0("log2(RNA) ,", ylab)) + theme + scale_color_distiller(palette = "RdPu") + theme(legend.position = "none") + ylim(-1,12) + xlim(-10,10) +
          annotate("text", x = 5, y = 0, label = paste0("rho=", round(rna.cor[i,j], 2)))
        p2[[i]] <- ggplot(dplot.2, aes(x = log2(x), y = log2(y))) + geom_point(aes(colour=  log10(y)), alpha = 0.2) + xlab(paste0("log2(protein) , ", xlab)) + ylab(paste0("log2(RNA*PTR) ,", ylab)) + theme + geom_smooth(method = "lm", fill = "#41b6c4", col = "#2c7fb8") + scale_color_distiller(palette = "YlGnBu") + theme(legend.position = "none") + ylim(-10, 10) + xlim(-10,10) +
          annotate("text", x = 5, y = -5, label = paste0("rho=", round(ptr.cor[i,j], 2)))
        p3[[i]] <- ggplot(dplot.3, aes(x = log2(x), y = log2(z))) + geom_point(aes(colour=  log10(y)), alpha = 0.2) + xlab(paste0("log2(protein) , ", xlab)) + ylab(paste0("log2(protein, median)")) + theme + geom_smooth(method = "lm", fill = "#a6bddb", col = "#2b8cbe") + scale_color_distiller(palette = "PuBu") + theme(legend.position = "none") + ylim(-10,10) + xlim(-10,10) +
          annotate("text", x = 5, y = -5, label = paste0("rho=", round(prmed.cor[i,j], 2)))
      }
    }
  }
  pdf(paste0(path, rep, "_prediction.pdf"), height = 3, width = 3*9)
  do.call("grid.arrange", c(p1[1:9], nrow= 1))
  do.call("grid.arrange", c(p2[1:9], nrow= 1))
  do.call("grid.arrange", c(p3[1:9], nrow= 1))
  dev.off()
  
  pdf(paste(path, rep, "_rna-prediction_lf.pdf"), height = 5, width = 5)
  heatmap.2(rna.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  heatmap.2(ptr.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  heatmap.2(prmed.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  dev.off()
  
  
  
  df[,30:38] <- df[,10:18]*df$rtp.median
  tsne <- log2(df[,c(1:9, 30:38)])
  tsne <- fill.na(tsne, columns = 1:18, cutoff = 0)
  tsne <- tsne[!rowSums(is.na(tsne)) > 0,] # exclude rows with na
  
  # TSNE
  tsne_out = Rtsne(log2(t(tsne)), perplexity = 5)
  tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], type = c(rep("Protein (absolute)", 9), rep("RNA (tpm*ptr)", 9)),  col = gsub("\\.t.*", "", colnames(df2[,c(1:18)])))
  # print PCA
  pdf(paste0("output/", today, "/", rep, "_tsne_lf_predicted.pdf"), width = 4.5, height = 3, useDingbats=FALSE)
  print(ggplot(tsne_plot, aes(x = x, y = y, col = col)) + geom_point(aes(shape = factor(type))) + theme + ggtitle("tSNE") + scale_color_brewer(palette = "Spectral"))
  dev.off()
}

rtp.calc.publication <- function(protein, rna, rep = "rep3"){ # for 9 cell-lines
  df <- data.frame(protein, rna)
  df$gene <- rownames(df)
  df <- data.frame(df, df[,1:9]/df[,10:18])
  colnames(df) <- c(colnames(df)[1:18], "gene", paste0(cells, ".rtp"))
  df$rtp.median <- apply(df[,20:(20+8)], 1, median, na.rm = T)
  df$prot.median <- apply(df[,1:9], 1, median, na.rm = T)
  write.csv(df, paste0(path, "tables/", rep, "labelfree.csv"))
  empty <-matrix(NA, nrow = 9, ncol = 9)
  
  
  rna.cor <- empty
  ptr.cor <- empty
  prmed.cor <- empty
  met <- "spearman"
  for(i in 1:9){
    for(j in 1:9){
      rna.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      ptr.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]*df$rtp.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      prmed.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df$prot.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw for protein median prediction
      #rna.pval[i,j] <- cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]))$p.val # rna raw
      #ptr.cor[i,j] <-  cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]*df$rtp.median))$estimate # rna raw
      #ptr.pval[i,j] <- cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]*df$rtp.median))$p.val # rna raw
      dplot.1 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1))]
      dplot.2 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 29)]
      dplot.2[,2] <- dplot.2[,2]*dplot.2[,3]
      dplot.3 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 30)]
      xlab = colnames(df)[(1 + (i-1) * 1)]
      ylab = colnames(df)[(10 + (j-1) * 1)]
      colnames(dplot.1) <- c("x", "y")
      colnames(dplot.2) <- c("x", "y")
      colnames(dplot.3) <- c("x", "y", "z")
      if(i == j){
        p1[[i]] <- ggplot(dplot.1, aes(x = log2(x), y = log2(y))) + geom_point(aes(colour=  log10(y)), alpha = 0.1) + xlab(paste0("log2(protein)")) + ylab(paste0("log2(RNA)")) + theme + scale_color_distiller(palette = "RdPu") + theme(legend.position = "none") + ylim(2,12) + xlim(-10,10) +
          annotate("text", x = 7, y = 3, label = paste0(round(rna.cor[i,j], 2))) + ggtitle(cells[i])
        p2[[i]] <- ggplot(dplot.2, aes(x = log2(x), y = log2(y))) + geom_point(aes(colour=  log10(y)), alpha = 0.1) + xlab(paste0("log2(protein)")) + ylab(paste0("log2(RNA*PTR)")) + theme + geom_smooth(method = "lm", fill = "#41b6c4", col = "#2c7fb8") + scale_color_distiller(palette = "YlGnBu") + theme(legend.position = "none") + ylim(-10, 10) + xlim(-10,10) +
          annotate("text", x = 7, y = -8.4, label = paste0(round(ptr.cor[i,j], 2)))
        p3[[i]] <- ggplot(dplot.3, aes(x = log2(x), y = log2(z))) + geom_point(aes(colour=  log10(y)), alpha = 0.1) + xlab(paste0("log2(protein)")) + ylab(paste0("log2(protein, median)")) + theme + geom_smooth(method = "lm", fill = "#a6bddb", col = "#2b8cbe") + scale_color_distiller(palette = "PuBu") + theme(legend.position = "none") + ylim(-10,10) + xlim(-10,10) +
          annotate("text", x = 7, y = -8.4, label = paste0(round(prmed.cor[i,j], 2)))
      }
    }
  }
  pdf(paste0(path, rep, "_prediction.pdf"), height = 2, width = 2*9)
  do.call("grid.arrange", c(p1[1:9], nrow= 1))
  do.call("grid.arrange", c(p2[1:9], nrow= 1))
  do.call("grid.arrange", c(p3[1:9], nrow= 1))
  dev.off()
  
  pdf(paste(path, rep, "_rna-prediction_lf.pdf"), height = 5, width = 5)
  heatmap.2(rna.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  heatmap.2(ptr.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  heatmap.2(prmed.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  dev.off()
  
  
  
  df[,30:38] <- df[,10:18]*df$rtp.median
  tsne <- log2(df[,c(1:9, 30:38)])
  tsne <- fill.na(tsne, columns = 1:18, cutoff = 0)
  tsne <- tsne[!rowSums(is.na(tsne)) > 0,] # exclude rows with na
  
  # TSNE
  tsne_out = Rtsne(log2(t(tsne)), perplexity = 5)
  tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], type = c(rep("Protein (absolute)", 9), rep("RNA (tpm*ptr)", 9)),  col = gsub("\\.t.*", "", colnames(df2[,c(1:18)])))
  # print PCA
  pdf(paste0("output/", today, "/", rep, "_tsne_lf_predicted.pdf"), width = 4.5, height = 3, useDingbats=FALSE)
  print(ggplot(tsne_plot, aes(x = x, y = y, col = col)) + geom_point(aes(shape = factor(type))) + theme + ggtitle("tSNE") + scale_color_brewer(palette = "Spectral"))
  dev.off()
}


# Function below
rtp.calc.abs <- function(protein, rna, rep = "rep3"){ # for 9 cell-lines
  df <- data.frame(protein, rna)
  df$gene <- rownames(df)
  df <- data.frame(df, df[,1:9]/df[,10:18])
  colnames(df) <- c(colnames(df)[1:18], "gene", paste0(cells, ".rtp"))
  df$rtp.median <- apply(df[,20:(20+8)], 1, median, na.rm = T)
  df$prot.median <- apply(df[,1:9], 1, median, na.rm = T)
  write.csv(df, paste0(path, "tables/", rep, "_ptr.csv"))
  empty <-matrix(NA, nrow = 9, ncol = 9)
  
  
  rna.cor <- empty
  ptr.cor <- empty
  prmed.cor <- empty
  random.rna <- empty
  met <- "spearman"
  for(i in 1:9){
    for(j in 1:9){
      rna.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      ptr.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]*df$rtp.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      prmed.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df$prot.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw for protein median prediction
      #rna.pval[i,j] <- cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]))$p.val # rna raw
      #ptr.cor[i,j] <-  cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]*df$rtp.median))$estimate # rna raw
      #ptr.pval[i,j] <- cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]*df$rtp.median))$p.val # rna raw
      dplot.1 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1))]
      dplot.2 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 29)]
      dplot.2[,2] <- dplot.2[,2]*dplot.2[,3]
      dplot.3 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 30)]
      xlab = colnames(df)[(1 + (i-1) * 1)]
      ylab = colnames(df)[(10 + (j-1) * 1)]
      colnames(dplot.1) <- c("x", "y")
      colnames(dplot.2) <- c("x", "y")
      colnames(dplot.3) <- c("x", "y", "z")
      if(i == j){
        p1[[i]] <- ggplot(dplot.1, aes(x = log2(x), y = log2(y))) + geom_point(aes(colour=  log10(y)), alpha = 0.2) + xlab(paste0("log2(protein)")) + ylab(paste0("log2(RNA)")) + theme + scale_color_distiller(palette = "RdPu") + theme(legend.position = "none") + ylim(2,12) + xlim(8,24) +
          annotate("text", x = 20, y = 3, label = paste0(round(rna.cor[i,j], 2)))
        p2[[i]] <- ggplot(dplot.2, aes(x = log2(x), y = log2(y))) + geom_point(aes(colour=  log10(y)), alpha = 0.2) + xlab(paste0("log2(protein)")) + ylab(paste0("log2(RNA*PTR)")) + theme + geom_smooth(method = "lm", fill = "#41b6c4", col = "#2c7fb8") + scale_color_distiller(palette = "YlGnBu") + theme(legend.position = "none") + ylim(8,24) + xlim(8,24) +
          annotate("text", x = 20, y = 9, label = paste0(round(ptr.cor[i,j], 2)))
        p3[[i]] <- ggplot(dplot.3, aes(x = log2(x), y = log2(z))) + geom_point(aes(colour=  log10(y)), alpha = 0.2) + xlab(paste0("log2(protein)")) + ylab(paste0("log2(protein median)")) + theme + geom_smooth(method = "lm", fill = "#a6bddb", col = "#2b8cbe") + scale_color_distiller(palette = "PuBu") + theme(legend.position = "none") + ylim(8,24) + xlim(8,24) +
          annotate("text", x = 20, y = 9, label = paste0(round(prmed.cor[i,j], 2)))
      }
    }
  }
  pdf(paste(path, rep, "_rna-prediction.pdf"), height = 5, width = 5)
  heatmap.2(rna.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  heatmap.2(ptr.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  heatmap.2(prmed.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  dev.off()
  
  pdf(paste0(path, rep, "_prediction.pdf"), height = 2, width = 2*9)
  do.call("grid.arrange", c(p1[1:9], nrow= 1))
  do.call("grid.arrange", c(p2[1:9], nrow= 1))
  do.call("grid.arrange", c(p3[1:9], nrow= 1))
  dev.off()
  
  df[,30:38] <- df[,10:18]*df$rtp.median
  tsne <- log2(df[,c(1:9, 30:38)])
  tsne <- fill.na(tsne, columns = 1:18, cutoff = 0)
  tsne <- tsne[!rowSums(is.na(tsne)) > 0,] # exclude rows with na
  
  # TSNE
  tsne_out = Rtsne(log2(t(tsne)), perplexity = 5)
  tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], type = c(rep("Protein (absolute)", 9), rep("RNA (tpm*ptr)", 9)),  col = gsub("\\.t.*", "", colnames(df2[,c(1:18)])))
  # print PCA
  pdf(paste0("output/", today, "/", rep, "_tsne_absolute_tpm_predicted.pdf"), width = 4.5, height = 3, useDingbats=FALSE)
  print(ggplot(tsne_plot, aes(x = x, y = y, col = col)) + geom_point(aes(shape = factor(type))) + theme + ggtitle("tSNE") + scale_color_brewer(palette = "Spectral"))
  dev.off()
}




# FINAL




  
  stop()
  
  
  for(i in 1:9){
    for(j in 1:9){
      rna.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      ptr.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df[,(10 + (j-1) * 1)]*df$rtp.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw
      prmed.cor[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df$prot.median), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw for protein median prediction
      random.rna[i,j] <- cor.test(log2(df[,(1 + (i-1) * 1)]),log2(df$rand), use = "pairwise.complete.obs", na.action = "na.exclude", method = met)$estimate # rna raw for protein median prediction
      #rna.pval[i,j] <- cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]))$p.val # rna raw
      #ptr.cor[i,j] <-  cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]*df$rtp.median))$estimate # rna raw
      #ptr.pval[i,j] <- cor.test(log2(df[,1*(i-1)+1]),log2(df[,1*(i-1)+(9 + (1*(j-1)))]*df$rtp.median))$p.val # rna raw
      dplot.1 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1))]
      dplot.2 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 29)]
      dplot.2[,2] <- dplot.2[,2]*dplot.2[,3] # multiply with PTR
      dplot.3 <- df[,c((1 + (i-1) * 1), (10 + (j-1) * 1), 30)] # only median protein
      xlab = colnames(df)[(1 + (i-1) * 1)]
      ylab = colnames(df)[(10 + (j-1) * 1)]
      colnames(dplot.1) <- c("x", "y")
      colnames(dplot.2) <- c("x", "y")
      colnames(dplot.3) <- c("x", "y", "z")
      
      dplot.2$error.log2 <- log2(dplot.2$y/dplot.2$x)
      dplot.3$error.log2 <- log2(dplot.3$z/dplot.3$x)
      if(i == j){
        alpha.v <- 0.3
        p1[[i]] <- ggplot(dplot.1, aes(x = log2(x), y = log2(y))) + geom_point(aes(colour=  log10(y)), alpha = alpha.v) + xlab(expression(log[2] *"(protein)")) + ylab(expression(log[2] *"(TPM)")) + theme + scale_color_distiller(palette = "RdPu") + theme(legend.position = "none") + ylim(2,12) + xlim(-10,10) +
          annotate("text", x = 7, y = 3, label = paste0(round(rna.cor[i,j], 2))) + ggtitle(cells[i])
        p2[[i]] <- ggplot(dplot.2, aes(x = log2(x), y = log2(y))) + geom_point(aes(colour=  log10(y)), alpha = alpha.v) + xlab(expression(log[2] *"(protein)")) + ylab(expression(log[2] *"(PTR-based pred.)")) + theme + scale_color_distiller(palette = "YlGnBu") + theme(legend.position = "none") + ylim(-10, 10) + xlim(-10,10) +
          annotate("text", x = 7, y = -8.4, label = paste0(round(ptr.cor[i,j], 2))) + geom_abline(intercept = 0, slope = 1, lty = 2)
        p3[[i]] <- ggplot(dplot.3, aes(x = log2(x), y = log2(z))) + geom_point(aes(colour=  log10(y)), alpha = alpha.v) + xlab(expression(log[2] *"(protein)")) + ylab(expression(log[2] *"(RNA-free pred.)")) + theme + scale_color_distiller(palette = "PuBu") + theme(legend.position = "none") + ylim(-10,10) + xlim(-10,10) +
          annotate("text", x = 7, y = -8.4, label = paste0(round(prmed.cor[i,j], 2))) + geom_abline(intercept = 0, slope = 1, lty = 2)
        pdf(paste0(path, "error-regression_", cells[i], "_", rep, ".pdf"), height = 3, width = 3)
        hist(dplot.2$error, breaks = seq(-100000, 100000, by = 0.5), xlim = c(-5,5), col = rgb(0.25,0.71,0.76,0.2), main = paste0(cells[i], " - residual error (%)"), freq = F)
        abline(v = median(dplot.2$error, na.rm = T),col="cyan")
        hist(dplot.3$error,breaks = seq(-100000, 100000, by = 0.5), xlim = c(-5,5), add = T, col = rgb(0.172,0.1,0.760,0.2), freq = F)
        abline(v = median(dplot.3$error, na.rm = T),col="blue")
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
  write.csv(prmed.cor, "temp.csv")
  pdf(paste0(path, rep, "_prediction.pdf"), height = 2, width = 2*9)
  do.call("grid.arrange", c(p1[1:9], nrow= 1))
  do.call("grid.arrange", c(p2[1:9], nrow= 1))
  do.call("grid.arrange", c(p3[1:9], nrow= 1))
  dev.off()
  
  # for caluclation only
  # length(error.ptr[which(error.ptr < 1 & error.ptr > -1)]) / length(which(!is.na(error.ptr)))
  # length(error.free[which(error.free < 1 & error.free > -1)]) / length(which(!is.na(error.free)))
  
  print(paste0("Total, RMSE, ptr=", 
               round(rmse(log2(rmse.all.ptr[,1]), log2(rmse.all.ptr[,2])), 2), 
               ", free=",
               round(rmse(log2(rmse.all.free[,1]), log2(rmse.all.free[,2])), 2),
               "; MAE, ptr = ",
               round(mae(log2(rmse.all.ptr[,1]), log2(rmse.all.ptr[,2])), 2), 
               ", free=",
               round(mae(log2(rmse.all.free[,1]), log2(rmse.all.free[,2])), 2)))
  
  pdf(paste0(path, "error-regression-all_", rep, ".pdf"), height = 4, width = 4)
  hist(error.ptr, breaks = seq(-100000, 100000, by = .5), xlim = c(-5,5), add = F, col = rgb(0.411,0.741,0.27, 0.2), freq = F)
  hist(error.free, breaks = seq(-100000, 100000, by = .5), xlim = c(-5,5), add = T, col = rgb(0.341,0.549,0.682,0.2), freq = F)
  lines(density(error.ptr[!is.na(error.ptr)], adjust = 2), col = rgb(0.411,0.741,0.27,1), lwd = 3)
  lines(density(error.free[!is.na(error.free)], adjust = 2), col = rgb(0.341,0.549,0.682,1), lwd = 3)
  abline(v = median(error.free, na.rm = T),col="blue")
  abline(v = -1 ,col="gray", lwd = 2, lty = 2)
  abline(v = 1 ,col="gray", lwd = 2, lty = 2)
  dev.off()
  
  pdf(paste(path, rep, "_rna-prediction_lf.pdf"), height = 5, width = 5)
  heatmap.2(rna.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  heatmap.2(ptr.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  heatmap.2(prmed.cor, Rowv = F, Colv = F, xlab = "Cell-lines (protein copies)", ylab = "Cell-lines (TPM)", main = "Direct correlation between datasets (protein and rna)",
            key.title = "Spearman rho", key.xlab = "Spearman rho", margins=c(7,7), cexRow = 1, cexCol = 1, bk <- seq(0, 1, by=0.025), symkey = F, symm = F, symbreaks = F, scale = "none")
  dev.off()
  



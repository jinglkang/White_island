# calculate FDR P values from EVE LRT test,  plot results and save results to .xlsx and .pdf files
# takes as input a *TestLRTs.res result file from EVE
# load xlsx package to write results to excel file
require(xlsx)

# set directory to save output files to
setwd("~/Documents/2019/香港大学/co2_seeps/EVE_release/CO2_seep/total_exclu_mid")
# usage: results <- EVEresults(genes=<vector_of_gene_names>, LRT=<betaTestLRTs.res>, rate=<FDR>, output =<file_name>)


# define EVEresults function
EVEresults <- function(genes, LRT, rate, output){
  
  LRT <- t(as.matrix(LRT))
  df <- 1
  title <- expression(paste(beta[i]!=beta[shared]))
  P_diverge <- 1 - pchisq(LRT, df = df)
  P_diverse <- pchisq(LRT, df = df)
  FDR_diverge <- p.adjust(P_diverge, "fdr")
  FDR_diverse <- p.adjust(P_diverse, "fdr")
  EVE_diverge <- as.data.frame(cbind(genes, LRT, P_diverge, FDR_diverge))
  EVE_diverse <- as.data.frame(cbind(genes, LRT, P_diverse, FDR_diverse))
  colnames(EVE_diverge) <- c("Gene", "LRT", "P", "FDR")
  colnames(EVE_diverse) <- c("Gene", "LRT", "P", "FDR")
  EVE_diverge.sub <- subset(EVE_diverge, EVE_diverge$FDR < rate)
  EVE_diverge.sub <- as.data.frame(EVE_diverge.sub)
  colnames(EVE_diverge.sub)<- c("Gene", "LRT", "P", "FDR")
  EVE_diverse.sub <- subset(EVE_diverse, EVE_diverse$FDR < rate)
  EVE_diverse.sub <- as.data.frame(EVE_diverse.sub)
  colnames(EVE_diverse.sub)<- c("Gene", "LRT", "P", "FDR")
  FDR_diverge.len <- length(EVE_diverge.sub[,1])
  EVE_diverge.len <- length(EVE_diverge[,1])
  FDR_diverse.len <- length(EVE_diverse.sub[,1])
  EVE_diverse.len <- length(EVE_diverse[,1])
  cat(paste0(FDR_diverge.len, " genes have higher variance among than within lineages at ", rate, " FDR."),"\n")
  cat(paste0(FDR_diverse.len, " genes have higher variance within than among lineages at ", rate, " FDR."),"\n")
  diverge_result <- list(EVE_diverge, EVE_diverge.sub)
  diverse_result <- list(EVE_diverse, EVE_diverse.sub)
  wd <- getwd()
  cat(paste0("Writing results to ", wd),"\n")
  
  #write results to 2x .csv files if too many genes for .xslx, otherwise write to 2 sheets of .xlsx
  #Note. I haven't tested exactly how many rows write.xlsx can cope with, it is likely more than the 5000 I have limited this to.
  if (FDR_diverge.len >5000) {
    write.csv(diverge_result[1], file=paste0(output,"_EVE_diverge_Pvalues.csv"))
    write.csv(diverge_result[2], file=paste0(output,"_EVE_diverge_FDRgenes.csv"))
    
  } else {
    if (FDR_diverge.len >0) {
      write.xlsx(diverge_result[1], file=paste0(output,"_EVE_diverge_result.xlsx"), sheetName="Gene Pvalues")
      write.xlsx(diverge_result[2], file=paste0(output,"_EVE_diverge_result.xlsx"), sheetName="FDRgenes", append=TRUE)
    } else {
      write.xlsx(diverge_result[1], file=paste0(output,"_EVE_diverge_result.xlsx"), sheetName="Gene Pvalues")
    }
  }
  
  if (FDR_diverse.len >5000) {
    write.csv(diverse_result[1], file=paste0(output,"_EVE_diverse_result_Pvalues.csv"))
    write.csv(diverse_result[2], file=paste0(output,"_EVE_diverse_result_FDRgenes.csv"))
    
  } else {
    if (FDR_diverse.len >0) {
      write.xlsx(diverse_result[1], file=paste0(output,"_EVE_diverse_result.xlsx"), sheetName="Gene Pvalues")
      write.xlsx(diverse_result[2], file=paste0(output,"_EVE_diverse_result.xlsx"), sheetName="FDRgenes", append=TRUE)
    } else {
      write.xlsx(diverse_result[1], file=paste0(output,"_EVE_diverse_result.xlsx"), sheetName="Gene Pvalues")
    }
  }
  
  index.lim <- length(EVE_diverge[,1])
  index.lim <- ceiling(index.lim/500)*500
  #fig <- plot.default(EVE_diverge$LRT, col=ifelse(EVE_diverge$FDR < rate, "red", "black"), main = title, xlim = c(0, index.lim), xaxt = "n")
  #axis(side=1,at=pretty(seq(0, index.lim, by=1000)),labels=pretty(seq(0, index.lim, by=1000)))
  pdf(file = paste0(output,"_EVE_diverge_result.pdf"), height = 6, width = 6)
  plot.default(EVE_diverge$LRT, col=ifelse(EVE_diverge$FDR < rate, "red", "black"), main = title, xlim = c(0, index.lim), xaxt = "n")
  axis(side=1,at=pretty(seq(0, index.lim, by=1000)),labels=pretty(seq(0, index.lim, by=1000)))
  dev.off()
  
  pdf(file = paste0(output,"_EVE_diverse_result.pdf"), height = 8.27, width = 11.69)
  plot.default(EVE_diverse$LRT, col=ifelse(EVE_diverse$FDR < rate, "red", "black"), main = title, xlim = c(0, index.lim), xaxt = "n")
  axis(side=1,at=pretty(seq(0, index.lim, by=1000)),labels=pretty(seq(0, index.lim, by=1000)))
  dev.off()
  
  res <- list()
  res$diverge_result <- diverge_result
  res$diverse_result <- diverse_result
  
  cat("\n", "Your analysis has finished, have a nice day.")
  
  return(res)
}


### test function using the example data
setwd("/path/to/test_data")

# load gene names and EVE result file
gene_names <- read.table("allgenes_names.txt")
LRT.res <- read.table("betaTestLRTs_total_exclu_mid.res")

# run analysis
results.LRT.res <- EVEresults(genes=gene_names, LRT=LRT.res, rate=0.05, output = "allgenes_data_results")






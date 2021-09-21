# White_island.TPM.TMM.sqrt.rename.matrix
# coldata_rename.txt
# coldata_rename_trait.txt
# coldata_rename_trait_record.txt
# coldata_record_info.txt

library(WGCNA)
setwd("~/Desktop/WGCNA/Blue_eyed/")
# Don't treat string as Factors, but my data has no string
options(stringsAsFactors = FALSE)

# read expression data
####################################################
white_island = read.table("Blue_eyed.TPM.TMM.sqrt.matrix", header = T, row.names = 1)
# transposed the expression data for the hclust
datExpr0=t(white_island)

# remove the bad genes and bad smaples
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average")
# plot the dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 500, col = "red")
# remove the outliers from the data, notice the cutHeight which is used to check the outlies
# minSize: minimum number of object on a branch to be considered a cluster
clust = cutreeStatic(sampleTree, cutHeight = 500, minSize = 10)
# show the individual number in each cluster after the cut
table(clust)
# select the cluter ID wanted to keep samples
keepSamples = (clust==1)
# extract these samples from the original data
datExpr = datExpr0[keepSamples, ]
####################################################

# read the trait data, and adjust it according to the column order of dataframe that removed outliers
####################################################
# the missing values are indicated as "NA"
traitData = read.table("coldata_Blue_eyed_record.txt", header = T, row.names = 1)
# match the order between the left dataframe that removed the outlier and traitData
White_island_Samples = rownames(datExpr)
traitRows = match(White_island_Samples, row.names(traitData))
# select the trait data according the index, and not include the first column [, -1]
datTraits = traitData[traitRows, -1]

# plot Sample dendrogram and trait heatmap
####################################################
traitColors = numbers2colors(datTraits, signed = FALSE)
sampleTree2 = hclust(dist(datExpr), method = "average")
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
####################################################




# network construction and module detection
####################################################
# Choose a set of soft-thresholding powers
par(mfrow=c(1, 2))
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
cex1 = 0.9
# Plot the results
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# add a line to check the lowest power for which the scale-free topology fit index reaches 0.90
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# calculate the adjacencies
####################################################
# select the lowest power value
softPower = 6
# Calculates (correlation or distance) network adjacency from the expression data
adjacency = adjacency(datExpr, power = softPower)
####################################################

# transform the adjacency into Topological Overlap Matrix, to minimize effects of noise and spurious associations
####################################################
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
# Clustering using TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree
sizeGrWindow(12, 9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
minModuleSize = 30
# deepSplit: integer (0-4), higher will produce the more and small clusters
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
# show the genes number in each cluster, modules labeled 1-[the end one] were from largest to smallest
# label 0 is reserved for unassigned genes
table(dynamicMods)
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")




# Merging of modules whose expression profiles are very similar
####################################################
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes before merge",
     xlab = "", sub = "")
# choose a height cut of 0.2, corresponding to correlation of 0.8
MEDissThres = 0.2
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs
# Calculate dissimilarity of merged module eigengenes
MEDiss_merged = 1-cor(mergedMEs)
# Cluster module eigengenes
METree_merged = hclust(as.dist(MEDiss_merged), method = "average")
plot(METree_merged, main = "Clustering of module eigengenes after merge",
     xlab = "", sub = "")


# plot the gene dendrogram again
sizeGrWindow(12, 9)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
####################################################




# Relating modules to external clinical traits
# identify modules that are significantly associated with the measured clinical traits
####################################################
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

# cor(): calculate pearson correlation between "MEs" and "datTraits"
# use="p": pairwise.complete.obs
moduleTraitCor = cor(MEs, datTraits, use = "p")
# Calculates Student asymptotic p-value for given correlations
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# plot the correlation and p-value of each module
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix_1 = paste(signif(moduleTraitCor, 2), " (",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
dim(textMatrix_1) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
textMatrix_1=as.data.frame(textMatrix_1)
names(textMatrix_1)=names(datTraits)
row.names(textMatrix_1)= names(MEs)
write.csv(textMatrix_1, file="textMatrix_1.csv")
write.csv(table(mergedColors),file = "mergedColors.csv")
####################################################

# Gene relationship to trait and important modules
# Gene Significance and Module Membership
####################################################
# Define variable "Species" containing the "Species" column of datTrait
Species = as.data.frame(datTraits$Species)
names(Species) = "Species"
# substring(x, start, stop): Extract or replace substrings in a character vector
# extract the string starting from the third position
modNames = substring(names(MEs), 3)

# extract the correlation and p-value between trait and the specific module into the dataframe
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, Species, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Species), sep="")
names(GSPvalue) = paste("p.GS.", names(Species), sep="")

# Identify genes with high GS (Gene significance) and MM (module membership)
# here select "greenyellow" module which shows the highest positive correlation with "Species"
module = "greenyellow"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), 
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Species",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# summary output of module analysis results
geneInfo0 = data.frame(moduleColor=moduleColors, geneTraitSignificance, GSPvalue)
nrow(geneInfo0)
# Order modules by their significance for "Species"
modOrder = order(-abs(cor(MEs, Species, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the gene
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Species));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file="geneInfo_species.csv")
###################################################################


# Visualizing the gene network: I will not do this
###################################################################
# NEVER try to cluster all genes
# heatmap
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 5);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

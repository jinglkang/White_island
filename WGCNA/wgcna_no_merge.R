library(WGCNA)
setwd("/media/HDD/white_island/wgcna/Blenny")
# Don't treat string as Factors, but my data has no string
options(stringsAsFactors = FALSE)

# read expression data
####################################################
white_island = read.table("Blenny_matrix_wgcna.xls", header = T, row.names = 1)
# transposed the expression data for the hclust
datExpr0=t(white_island)

sampleTree = hclust(dist(datExpr0), method = "average")
# plot the dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

datExpr = datExpr0
####################################################

# read the trait data, and adjust it according to the column order of dataframe that removed outliers
####################################################
# the missing values are indicated as "NA"
traitData = read.table("coldata_Blenny_wgcna.txt", header = T, row.names = 1)
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
                            deepSplit = 4, pamRespectsDendro = FALSE,
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



# Relating modules to external clinical traits
# identify modules that are significantly associated with the measured clinical traits
####################################################
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, dynamicColors)$eigengenes
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
write.csv(table(dynamicColors),file = "dynamicColors.csv")
####################################################

# Gene relationship to trait and important modules
# Gene Significance and Module Membership
####################################################
# Define variable "Species" containing the "Species" column of datTrait
pH = as.data.frame(datTraits$pH)
names(pH) = "pH"
# substring(x, start, stop): Extract or replace substrings in a character vector
# extract the string starting from the third position
modNames = substring(names(MEs), 3)

# extract the correlation and p-value between trait and the specific module into the dataframe
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, pH, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(pH), sep="")
names(GSPvalue) = paste("p.GS.", names(pH), sep="")

# Identify genes with high GS (Gene significance) and MM (module membership)
# here select "greenyellow" module which shows the highest positive correlation with "Species"
module = "orangered"
column = match(module, modNames)
moduleGenes = dynamicColors ==module
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

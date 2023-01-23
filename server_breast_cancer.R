# ##################################
# # WGCNA: module construction
# # identify which modules OBP2A/B are assigned
# # identify gene set(s) highly co-expressed with OBP2A/B
# # using WGCNA tutorial as reference to write the code and scripts
# #################
# # WGCNA preprocess data
# 
# preprocess_explore_data <- function(workdir = "D:/major_project",
#                                     input_file = "00_TCGA breast Cancer_norm_all_ExpectedCnt.csv")
# {
#   # dir.create("/home/lyang73/data/breast_cancer")
#   # data input and pre-process
#   setwd(workdir)
# 
#   # setwd("/home/lyang73/data")
#   library(GO.db)
#   # Load the WGCNA package
#   library(WGCNA)
# 
#   # check for genes and samples with too many missing values
#   # If the last statement returns TRUE, all genes have passed the cuts.
#   #
# 
#   exprFinal = read.csv(input_file,
#                        check.names = FALSE)
#   exprdata = as.data.frame(t(exprFinal[,-1:-3]))
#   colnames(exprdata) = exprFinal[,"gene"]
# 
#   expr_cancer = exprdata[!startsWith(rownames(exprdata), "GTEX"),]
#   expr_normal = exprdata[startsWith(rownames(exprdata), "GTEX"),]
# 
#   exprdata = expr_cancer
#   # Remove genes that are lowly-expressed and/or genes with low variation between samples.
#   # For WGCNA on TCGA gene expression data, we removed genes that has normalized
#   # expression value< O in all of our samples
# 
# 
#   #Add a new variable "Count" that counts the number of samples that have normalized expression value <0 for gene x
#   count = colSums(exprdata < 0);
#   table(count);
# 
# 
# 
# 
# 
#   # This function checks data for missing entries, entries with weights below a
#   # threshold, and zero-variance genes, and returns a list of samples and genes
#   # that pass criteria on maximum number of missing or low weight values. If
#   # necessary, the filtering is iterated.
#   gsg = goodSamplesGenes(exprdata, verbose = 3);
#   # BiocManager::install("genefilter")
#   # library(genefilter)
#   # test = t(exprdata)
#   # selProbes <- genefilter(test, filterfun(pOverA(0.20, log2(100)),
#   #                                         function(x) (IQR(x) > 0.25)))
#   # eset <- ALL[selProbes, ]
# 
# 
#   gsg$allOK
#   # If not, we remove the offending genes and samples
#   # from the data:
#   if (!gsg$allOK)
#   { # Optionally, print the gene and sample names that were removed:
#     if (sum(!gsg$goodGenes)>0)
#       printFlush(paste("Removing genes:", paste(names(exprdata)[!gsg$goodGenes], collapse = ", ")));
#     if (sum(!gsg$goodSamples)>0)
#       printFlush(paste("Removing samples:", paste(rownames(exprdata)[!gsg$goodSamples], collapse = ", ")));
#     # Remove the offending genes and samples from the data:
#     exprdata = exprdata[gsg$goodSamples, gsg$goodGenes]
#   }
# 
# 
#   # cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious
#   # outliers.
#   sampleTree = hclust(dist(exprdata), method = "average");
#   # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
#   # The user should change the dimensions if the window is too large or too small.
#   # sizeGrWindow(12,9)
#   pdf(file = "WGCNA_sample_Clustering.pdf", width = 20, height = 12);
#   par(cex = 0.6);
#   par(mar = c(0,4,2,0))
#   plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#        cex.axis = 1.5, cex.main = 2)
#   abline(h = 75, col = "red");
#   # Determine cluster under the line
#   clust = cutreeStatic(sampleTree, cutHeight = 75, minSize = 10)
# 
#   dev.off()
# 
#   table(clust)
#   # It appears there is one outlier (sample F2_221, see Fig. 1). One can remove it by hand, or use an automatic approach.
#   # Choose a height cut that will remove the offending sample, say 15 (the red line in the plot), and use a branch cut at
#   # that height.
#   # Plot a line to show the cut
# 
# 
# 
#   # Choose a set of soft-thresholding powers
#   powers = c(c(1:10), seq(from = 12, to=20, by=2))
#   # Call the network topology analysis function
#   sft = pickSoftThreshold(exprdata, powerVector = powers, verbose = 5)
#   # Plot the results:
#   pdf(file = "WGCNA_Analysis_of_network_topology.pdf", width = 18, height = 12);
# 
#   # sizeGrWindow(9, 5)
#   par(mfrow = c(1,2));
#   cex1 = 0.9;
#   # Scale-free topology fit index as a function of the soft-thresholding power
#   plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#        xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#        main = paste("Scale independence"));
#   text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#        labels=powers,cex=cex1,col="red");
#   # this line corresponds to using an R^2 cut-off of h
#   abline(h=0.8,col="red") # Mean connectivity as a function of the soft-thresholding power
#   plot(sft$fitIndices[,1], sft$fitIndices[,5],
#        xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#        main = paste("Mean connectivity"))
#   text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#   dev.off()
#   save(exprdata, file = "WGCNA_01-dataInput_mydata.RData")
#   return(exprdata)
# }
# # exprdata = preprocess_explore_data("D:/major_project", "00_TCGA breast Cancer_norm_all_ExpectedCnt.csv")
# 
# exprdata = preprocess_explore_data("/home/lyang73/data/breast_cancer",
#                                    "00_TCGA breast Cancer_norm_all_ExpectedCnt.csv")
# 
#########
# WGCNA Automatic network construction and module detection
# Constructing the gene network and identifying modules using simple function call:

# # disableWGCNAThreads()
# net = blockwiseModules(exprdata, power = 10,
#                        TOMType = "unsigned", minModuleSize = 30,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE, maxBlockSize = 20000,
#                        saveTOMFileBase = "femaleMouseTOM",
#                        verbose = 3)
# # 18232 genes, 1393 samples
#
# table(net$colors)
#
# # open a graphics window
# pdf(file = "WGCNA_Clustering_dendrogram_of_genes.pdf", width = 18, height = 12);
# # Convert labels to colors for plotting
# mergedColors = labels2colors(net$colors) # Plot the dendrogram and the module colors underneath
# plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# dev.off()
#
# moduleLabels = net$colors
# moduleColors = labels2colors(net$colors)
# MEs = net$MEs;
# geneTree = net$dendrograms[[1]];
# save(MEs, moduleLabels, moduleColors, geneTree,
#      file = "WGCNA_networkConstruction-auto_mydata.RData")
# head(moduleLabels)


for (val in seq(4,10))
{
  print(val)
  enableWGCNAThreads()
  net = blockwiseModules(exprdata, power = val,
                         TOMType = "unsigned", minModuleSize = 10,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE, maxBlockSize = 20000,
                         saveTOMFileBase = "breast_cancer",
                         verbose = 3)
  # 18232 genes, 1393 samples

  table(net$colors)

  # open a graphics window

  pdf(file = paste0("WGCNA_Clustering_dendrogram_of_genes",val,".pdf"),
      width = 18, height = 12);
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors) # Plot the dendrogram and the module colors underneath
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()

  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]]
  if(val == 10)
  {
    save(MEs, moduleLabels, moduleColors, geneTree,
         file = paste0("WGCNA_networkConstruction-auto_mydata",val,".RData"))
  }else
  {
    save(MEs, moduleLabels, moduleColors, geneTree,
         file = paste0("WGCNA_networkConstruction-auto_mydata0",val,".RData"))
  }
  print(paste0("WGCNA_networkConstruction-auto_mydata0",val,".RData"))
}


##########
# summary which modules interest genes (OBPs) are locate/assigned
# the gene numbers contained in each module under different runs (soft threshold 5-10)
# setwd("D:/major_project/breast_cancer")
setwd("/home/lyang73/data/breast_cancer")

summarize_module_information <- function(
  input_file = "WGCNA_networkConstruction-auto_mydata",
  output_symbol = "auto")
{
  load("WGCNA_01-dataInput_mydata.RData")
  gene_assign <- data.frame("power" = seq(4,10), "OBP2B" = rep(0,7),
                            "OBP2A" = rep(0,7),"LCN2" = rep(0,7))

  module_num <- as.data.frame(matrix(rep(0,14*32),nrow=14,ncol=32,byrow=TRUE))
  colnames(module_num) = c("power",rep("module",31))
  module_num$power = c(4,4,5,5,6,6,7,7,8,8,9,9,10,10)


  for (val in seq(4,10))
  {
    print(val)
    if(val == 10){
      load(file = paste0(input_file,val,".RData"))
    }else{
      load(file = paste0(input_file,"0",val,".RData"))
    }
    names(moduleLabels) = colnames(exprdata)

    gene_assign[gene_assign$power==val,"OBP2B"] = moduleLabels[names(moduleLabels) == "OBP2B"]
    gene_assign[gene_assign$power==val,"OBP2A"] = moduleLabels[names(moduleLabels) == "OBP2A"]
    gene_assign[gene_assign$power==val,"LCN2"] = moduleLabels[names(moduleLabels) == "LCN2"]
    if(length(unique(moduleLabels)) <= 31)
    {
      module_num[module_num$power==val,1:length(unique(moduleLabels))+1] = t(as.data.frame(table(moduleLabels)))

      # module_num[module_num$power==val,1:length(unique(moduleLabels))+1] = as.data.frame(table(moduleLabels))$Freq

    }else{
      module_num[module_num$power==val,2:32] = t(as.data.frame(table(moduleLabels)[1:31]))

      # module_num[module_num$power==val,1:32] = as.data.frame(table(moduleLabels))$Freq

    }
    print(table(moduleLabels))
    print("finished")
    rm(moduleLabels,moduleColors)
  }
  write.csv(gene_assign,file = paste0("WGCNA_gene_assign_",output_symbol,".csv"),
            row.names = FALSE)

  write.csv(module_num,file = paste0("WGCNA_module_summary_",output_symbol,".csv"),
            row.names = FALSE)
  output_file_names = c(paste0("WGCNA_gene_assign_",output_symbol,".csv"),
                        paste0("WGCNA_module_summary_",output_symbol,".csv"))
  return(output_file_names)

}
# stepbystep
output_file_names = summarize_module_information(input_file = "WGCNA_networkConstruction-auto_mydata",
                                                 output_symbol = "auto")
print(output_file_names)

# load(file = "WGCNA_01-dataInput_mydata.RData")
# library(WGCNA)
# rm(moduleLabels)
# load(file = "02-networkConstruction-stepByStep09.RData")
# names(moduleLabels) = colnames(exprdata)
# moduleLabels[names(moduleLabels) == "OBP2B"]
# moduleLabels[names(moduleLabels) == "OBP2A"]
# table(moduleLabels)





########
# Exporting a list of  genes within a module and their connectivity scores
# scores used for gene network visualization and enrichment analysis
# Exporting a gene network within one module to the file formate compatible with network visualization software

# D:/major_project/breast_cancer  /home/lyang73/data/breast_cancer

retrive_gene_module_scores <- function(interest_module = 10,
                                       workdir = "D:/major_project/breast_cancer",
                                       interest_gene = "OBP2B",
                                       WGCNA_result_file = "WGCNA_networkConstruction-auto_mydata08.RData")
{
  setwd(workdir)
  library(WGCNA)


  load(file = WGCNA_result_file)
  load(file = "WGCNA_01-dataInput_mydata.RData")
  enableWGCNAThreads()
  # Recalculate topological overlap
  # TOM = TOMsimilarityFromExpr(exprdata, power = 8);
  # save(TOM,file = "TOM.RData")
  load("TOM.RData")

  test = (names(moduleLabels) == interest_gene )
  test1 = TOM[test, ];
  names(test1) = names(moduleLabels)
  interest_gene_score = test1
  save(interest_gene_score,file = paste0(interest_gene,"_score.RData"))


  # Select module
  module = interest_module;
  # Select module probes
  probes = names(exprdata)
  inModule = (moduleLabels==module);
  # inModule = (moduleColors==module);
  modProbes = probes[inModule];
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];

  dimnames(modTOM) = list(modProbes, modProbes)
  # Export the network into an edge list file VisANT can read
  vis = exportNetworkToVisANT(modTOM,
                              file = paste("VisANTInput-", module, ".txt", sep=""),
                              weighted = TRUE,
                              threshold = 0)
  #
  # Because the brown module is rather large, we can restrict the genes in
  # the output to say the 30 top hub genes in the module:


  nTop = 30;
  IMConn = softConnectivity(exprdata[, modProbes],power = 8);
  top = (rank(-IMConn) <= nTop)
  vis = exportNetworkToVisANT(modTOM[top, top],
                              file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                              weighted = TRUE,
                              threshold = 0 )

  output_file_names = c(paste0(interest_gene,"_score.RData"),
                        paste("VisANTInput-", module, ".txt", sep=""),
                        paste("VisANTInput-", module, "-top30.txt", sep=""))
  return(output_file_names)
}


output_file_names = retrive_gene_module_scores(interest_gene = "OBP2B",
                                               interest_module = 7,
                                               WGCNA_result_file = "WGCNA_networkConstruction-auto_mydata05.RData" )
print(output_file_names)
output_file_names = retrive_gene_module_scores(interest_module = 7,
                                               interest_gene = "LCN2",
                                               WGCNA_result_file = "WGCNA_networkConstruction-auto_mydata05.RData" )

output_file_names = retrive_gene_module_scores(interest_module = 11,
                                               interest_gene = "OBP2A",
                                               WGCNA_result_file = "WGCNA_networkConstruction-auto_mydata05.RData" )





#
# Export the network into an edge list file Cytoscape can read
# # Recalculate topological overlap if needed
# TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# # Read in the annotation file
# annot = read.csv(file = "GeneAnnotation.csv");
# # Select modules
# modules = c("brown", "red");
# # Select module probes
# probes = names(datExpr)
# inModule = is.finite(match(moduleColors, modules));
# modProbes = probes[inModule];
# modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# # # Select the corresponding Topological Overlap
# # modTOM = TOM[inModule, inModule];
# # dimnames(modTOM) = list(modProbes, modProbes)
# # Export the network into edge and node list files Cytoscape can read
# cyt = exportNetworkToCytoscape(modTOM,
#                                edgeFile = paste("CytoscapeInput-edges-",".txt", sep=""),
#                                nodeFile = paste("CytoscapeInput-nodes-",".txt", sep=""),
#                                weighted = TRUE,
#                                threshold = 0.02,
#                                nodeNames = modProbes,
#                                nodeAttr = moduleColors[inModule])



###########
# plot gene co-expression network, pick top 30 co-expressed genes 
# based on connectivity score between interest gene with all other genes
# plot connectivity scores distribution histogram 

# install.packages("networkD3")

# Read a data set. 
# Data format: dataframe with 3 variables; variables 1 & 2 correspond to
# interactions; variable 3 is weight of interaction


plot_network <- function (input_file = "VisANTInput-10.txt",
                          interest_gene = "OBP2B", interest_module = 10,
                          workdir = "/home/lyang73/data/breast_cancer")
{
  setwd(workdir)
  library(networkD3)
  library(ggplot2)
  edgeList <- read.table(input_file, header = FALSE, sep = " ")
  edgeList = edgeList[,c("V1" ,"V2", "V5")]
  colnames(edgeList) <- c("SourceName", "TargetName", "Weight")
  edgeList_raw =  edgeList
  
  # only contain edges connected with interest_gene
  edgeList = edgeList[edgeList$SourceName == interest_gene |
                        edgeList$TargetName == interest_gene,]
  edgeList_OBP =  edgeList
  # only contain top 30 edges strongly connected with interest_gene
  edgeList = edgeList[sort(edgeList$Weight, decreasing = TRUE, 
                           index.return=TRUE)$ix[1:30],]
  
  write.csv(edgeList,paste0(interest_gene,"_top30_coexpressed_genes.csv"), 
            row.names = FALSE)
  
  n = 1
  while (n < 3) 
  {
    if(n == 1)
    {
      df = edgeList_raw
      n = n+1
      mytitle = paste0("Histogram of topological overlap score of edges within module ",
                       as.character(interest_module) )
      output_file_names = mytitle
      
    }else if(n == 2)
    {
      df = edgeList_OBP
      n = n+1
      mytitle = paste0("Histogram of topological overlap score of edges connected with ",
                       interest_gene)
      output_file_names = c(output_file_names, mytitle)
    }
    qts <- quantile(df$Weight,probs=c(.05,.95))
    
    # Histogram of scores of all edges within interest module 
    ggplot(df, aes(x=Weight)) + 
      geom_histogram(aes(y=after_stat(count)), colour="black", fill="white")+
      geom_density(alpha=.2, fill="#FF6666")+ xlab("topological overlap score")+
      geom_vline(aes(xintercept=qts[2]),color="red", linetype="dashed", 
                 linewidth=1)   +
      ggtitle(mytitle)
    # + scale_x_continuous(breaks=seq(0, 0.3, 0.05))
    ggsave(paste0(mytitle,".pdf"))
  }
  
  
  
  
  
  # Create a graph. Use simplyfy to ensure that there are no duplicated edges or self loops
  gD <- igraph::simplify(igraph::graph.data.frame(edgeList, directed=FALSE))
  
  # Create a node list object (actually a data frame object) that will contain information about nodes
  nodeList <- data.frame(ID = c(0:(igraph::vcount(gD) - 1)), # because networkD3 library requires IDs to start at 0
                         nName = igraph::V(gD)$name)
  nodeList$group = 1
  nodeList[nodeList$nName==interest_gene,"group"] = 2
  
  
  # Map node names from the edge list to node IDs
  getNodeID <- function(x){
    which(x == igraph::V(gD)$name) - 1 # to ensure that IDs start at 0
  }
  # And add them to the edge list
  edgeList <- plyr::ddply(edgeList, .variables = c("SourceName", "TargetName", "Weight"),
                          function (x) data.frame(SourceID = getNodeID(x$SourceName),
                                                  TargetID = getNodeID(x$TargetName)))
  # Let's create a network
  
  D3_network_interest_gene <- forceNetwork(Links = edgeList, # data frame that contains info about edges
                                           Nodes = nodeList, # data frame that contains info about nodes
                                           Source = "SourceID", # ID of source node
                                           Target = "TargetID", # ID of target node
                                           Value = "Weight", # value from the edge list (data frame) that will be used to value/weight relationship amongst nodes
                                           NodeID = "nName", # value from the node list (data frame) that contains node description we want to use (e.g., node name)
                                           height = 1000, # Size of the plot (vertical)
                                           width = 1000,  # Size of the plot (horizontal)
                                           Group = "group", # Font size
                                           opacity = 1, # opacity
                                           zoom = TRUE, # ability to zoom when click on the node
                                           opacityNoHover = 0.95,fontSize = 8,
                                           radiusCalculation = networkD3::JS("Math.sqrt(d.nodesize/5000)+1"),
                                           linkDistance = networkD3::JS("function(d) { return 8000*d.value; }"),
                                           linkWidth = networkD3::JS("function(d) { return d.value*70; }")) # opacity of labels when static
  
  # Plot network
  D3_network_interest_gene
  
  # Save network as html file
  
  networkD3::saveNetwork(D3_network_interest_gene, 
                         paste0("D3_network_top30_",interest_gene,".html"), 
                         selfcontained = TRUE)
  output_file_names = c(qts[2],output_file_names, paste0("D3_network_top30_",interest_gene,".html"))
  return(output_file_names)
}


significant_threshold = rep(0,3)
input_file = c("VisANTInput-11.txt","VisANTInput-7.txt","VisANTInput-7.txt")
interest_genes = c("OBP2A","OBP2B","LCN2")
module_assignment = c(11,7,7)
n = 1
while (n < 4) 
{
  output_file_names = plot_network(input_file[n],interest_genes[n],
                                   module_assignment[n],
                                   "D:/major_project/breast_cancer")
  significant_threshold[n] = output_file_names[1]
  n = n+1
}
# setwd("D:/major_project/breast_cancer")
save(significant_threshold,file = "significant_threshold.RData")
# D:/major_project/breast_cancer  /home/lyang73/data/breast_cancer

# output_file_names = plot_network("VisANTInput-7.txt","OBP2B",7,
#                                  "/home/lyang73/data/breast_cancer")
# # 0.005300881 OBP2B
# # print(output_file_names)
# # significant_threshold[2] = output_file_names[1]
# 
# output_file_names = plot_network("VisANTInput-11.txt","OBP2A",11,
#                                  "/home/lyang73/data/breast_cancer")
# 
# # 0.014264764778568


#############
# topGO: gene ontology enrichment analysis
# use the paper script as reference, plot the top 20 enriched GO entities  
# Interfacing genes in a module with gene ontology enrichment analysis


# setwd("/home/lyang73/data")

# D:/major_project/breast_cancer  /home/lyang73/data/breast_cancer

GO_enrichment_scores <- function(input_file = "OBP2A_score.RData",
                                 interest_genes = "OBP2A",
                                 significant_threshold = 0.01,
                                 workdir = "D:/major_project/breast_cancer",
                                 moduleLabels_file = "WGCNA_networkConstruction-auto_mydata05.RData",
                                 annot_file = "00_TCGA breast Cancer_norm_all_ExpectedCnt.csv")
{
  setwd(workdir)
  library(WGCNA)
  library(data.table)
  library(grex);
  library(biomaRt);
  library(topGO);
  library(dplyr);
  library(ggplot2)
  if( file.exists("topGO_moduleGeneAnnotated.csv"))
  {
    annotComplete = read.csv("topGO_moduleGeneAnnotated.csv")
  }else 
  {
    # Read in the probe annotation
    annot = read.csv(file = annot_file ,header = TRUE)
    annot = annot[,2:3]
    # Match probes in the data set to the probe IDs in the annotation file
    
    load(file = moduleLabels_file)
    
    table(moduleLabels)
    
    genes_module = names(moduleLabels)
    module_labels = as.data.frame(moduleLabels)
    module_labels$gene = rownames(module_labels)
    genes_module = merge(annot,module_labels,by = "gene")
    # write.csv(genes_module,file = "genes_module.csv", row.names = FALSE)
    # genes_module = read.csv("genes_module.csv")
    
    # # Annotate the Ensembl id of WGCNA result, with Entrez IDs, via function grex
    annotComplete = grex(cleanid(genes_module$id))
    genes_module$id = cleanid(genes_module$id)
    annotComplete = merge(annotComplete,genes_module[,-1],by.x = "ensembl_id",by.y = "id")
    
    write.csv(annotComplete, "topGO_moduleGeneAnnotated.csv", row.names = FALSE)
    
  }
  
  # Connect to the ENSEMBL_MART_ENSEMBL BioMart database to query GO IDs for our Entrez IDs
  # filter the annotComplete, delete the interest gene itself
  if(interest_genes == "OBP2A" | interest_genes == "OBP2B")
  {
    genes_bg = annotComplete[annotComplete$hgnc_symbol != "OBP2A" &
                               annotComplete$hgnc_symbol != "OBP2B", ]$entrez_id
  }else if(interest_genes == "LCN2")
  {
    genes_bg = annotComplete[annotComplete$hgnc_symbol != "LCN2" , ]$entrez_id
    
  }
  
  tot_background = length(genes_bg);
  
  db = useMart("ENSEMBL_MART_ENSEMBL",  dataset = "hsapiens_gene_ensembl",
               host = "https://www.ensembl.org");
  
  go_ids = getBM(attributes = c("go_id","entrezgene_id", "namespace_1003"),
                 filters = "entrezgene_id", values = genes_bg, mart = db)
  
  
  # Build a topGO data object with the base function new, then run GO analysis.                  
  gene2GO = unstack(go_ids[,c(1,2)])
  if(interest_genes == "OBP2A" | interest_genes == "OBP2B")
  {
    load(input_file)
    interest_gene_score = interest_gene_score[names(interest_gene_score) != "OBP2A" &
                                                names(interest_gene_score) != "OBP2B"]
  }else if(interest_genes == "LCN2")
  {
    load(input_file)
    interest_gene_score = interest_gene_score[names(interest_gene_score) != "LCN2"]
  }
  
  
  
  topDiffGenes <- function(allScore) {
    return(allScore > significant_threshold )
  }
  names(interest_gene_score) = genes_bg
  GOdata_interst_gene = new("topGOdata",
                            ontology = c("BP"), allGenes = interest_gene_score,
                            geneSel = topDiffGenes,
                            annot = annFUN.gene2GO,
                            gene2GO = gene2GO,  nodeSize = 5)
  save(GOdata_interst_gene,file = paste0("GOdata_",interest_genes,".RData"))  
  output_file = paste0("GOdata_",interest_genes,".RData")
  
  
  GOdata = GOdata_interst_gene
  allGO = usedGO(GOdata)
  # Test significance of GO terms using the function runTest.
  # algorithm and statistic specify the method for dealing with the GO graph 
  # structure and the test statistic, respectively.
  
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultks <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  
  
  # Generate a summary table of results obtained from topGO enrichment analysis with the function
  # GenTable.
  all_res_ks = GenTable(GOdata,  KS = resultks, orderBy = "KS", 
                        topNodes = length(allGO))
  
  all_res_ks$OR = log2((all_res_ks$Significant/numSigGenes(GOdata))/
                         (all_res_ks$Annotated/tot_background))
  
  # SigGenes(GOdata) return the genes that the method think is significant
  # there it return the top 36 genes have highest weights  
  # Significant column means the number of significant genes that been annotated 
  # in this GO entity  
  # the larger the value of OR, the more significant genes annotated into certain 
  # GO entity, the fewer the total annotated genes (including non-significant genes)
  
  # test = sigGenes(GOdata) 
  # test1 = annotComplete[annotComplete$entrez_id %in% test,]
  # annotComplete$weights = OBP2A_score
  # write.csv(all_res,"topGO_OBP_GOAnnotated.csv")
  
  # visualization of top 20 GO entities 
  # Generate a summary figure of topGO results via the ggplot2 R package
  #This sets the threshold p<0.05 and #.annotated.background.genes >= 30. 
  
  GO_bar = all_res_ks
  sapply(GO_bar,class);
  GO_bar = transform(GO_bar,KS = as.numeric(KS));
  sapply(GO_bar,class)
  
  GO_bar[is.na(GO_bar$KS),"KS"] =  1e-30
  GO_bar$adjust_KS = p.adjust(GO_bar$KS)
  
  GO_bar[is.nan(GO_bar$OR),"OR"] = 0
  GO_bar[GO_bar$OR == -Inf,"OR"] =  0
  GO_bar = GO_bar %>%
    arrange(adjust_KS) 
  GO_bar = GO_bar[1:20,]
  
  ggplot(GO_bar,
         aes(x = reorder(Term, -adjust_KS), y = Annotated, fill = adjust_KS))+
    geom_bar(stat = "identity", width = 0.9, color = "black") +
    coord_flip() +
    scale_fill_gradient(low="#feff2b",high="#fe0100")+
    ggtitle(paste("Enriched GO Biological Processes, classic ks",interest_genes))+
    labs( X = NULL,
          y = "number of annotated genes", fill = "p.value")+theme_bw()+
    theme(plot.title.position = "plot")+
    theme(plot.title = element_text(size = 12))+
    theme(axis.title.x = element_text(size = 5,face = "bold"),
          axis.title.y = element_text(size = 5,face = "bold"))+
    theme(axis.text.x = element_text(size = 5,face = "bold"),
          axis.text.y = element_text(size = 5,face = "bold"))+
    theme(legend.position = "right")
  ggsave(paste0("top 20 enriched GO for ",interest_genes," classic ks",".pdf"))
  output_file = c(output_file,paste0("top 20 enriched GO for ",interest_genes," classic ks",".pdf"))
  
  
  ## fisher results
  # sigTest = runTest(GOdata,
  #                   algorithm = "elim", statistic = "fisher")
  
  all_res_fisher = GenTable(GOdata,  Fis = resultFisher, orderBy = "Fis", 
                            topNodes = length(allGO))
  
  all_res_fisher$OR = log2((all_res_fisher$Significant/numSigGenes(GOdata))/
                             (all_res_fisher$Annotated/tot_background))
  
  GO_bar = all_res_fisher
  GO_bar = transform(GO_bar,Fis = as.numeric(Fis));
  GO_bar$adjust_fisher = p.adjust(GO_bar$Fis)
  
  GO_bar = GO_bar %>%
    arrange(adjust_fisher) 
  GO_bar = GO_bar[1:20,]
  
  ggplot(GO_bar,
         aes(x = reorder(Term, -adjust_fisher), y = OR, fill = adjust_fisher))+
    geom_bar(stat = "identity", width = 0.9, color = "black") +
    coord_flip() +
    scale_fill_gradient(low="#feff2b",high="#fe0100")+
    ggtitle(paste("Enriched GO Biological Processes, classic fisher",interest_genes))+
    labs(X = NULL,
         y = "Odds Ratio", fill = "p.value")+theme_bw()+
    theme(plot.title.position = "plot")+
    theme(plot.title = element_text(size = 12))+
    theme(axis.title.x = element_text(size = 5,face = "bold"),
          axis.title.y = element_text(size = 5,face = "bold"))+
    theme(axis.text.x = element_text(size = 5,face = "bold"),
          axis.text.y = element_text(size = 5,face = "bold"))+
    theme(legend.position = "right")
  ggsave(paste0("top 20 enriched GO for ",interest_genes," classic fisher",".pdf"))
  output_file = c(output_file,
                  paste0("top 20 enriched GO for ",interest_genes," classic fisher",".pdf"))
  return(output_file)
}

load("significant_threshold.RData")
output_file_names = GO_enrichment_scores(input_file = "OBP2A_score.RData",
                                         interest_genes = "OBP2A",
                                         significant_threshold = significant_threshold[1])

output_file_names = GO_enrichment_scores(input_file = "OBP2B_score.RData",
                                         interest_genes = "OBP2B",
                                         significant_threshold = significant_threshold[2])

output_file_names = GO_enrichment_scores(input_file = "LCN2_score.RData",
                                         interest_genes = "LCN2",
                                         significant_threshold = significant_threshold[3])

# D:/major_project/breast_cancer  /home/lyang73/data/breast_cancer

GO_enrichment_module <- function(interest_module = 5,
                                 workdir = "D:/major_project/breast_cancer",
                                 moduleLabels_file = "WGCNA_networkConstruction-auto_mydata05.RData")
{
  # # Read in the probe annotation
  # annot = read.csv(file = "00_TCGA breast Cancer_norm_all_ExpectedCnt.csv",header = TRUE)
  # annot = annot[,2:3]
  # # Match probes in the data set to the probe IDs in the annotation file
  # 
  # load(file = "WGCNA_networkConstruction-auto_mydata08.RData")
  # 
  # table(moduleLabels)
  # 
  # genes_module = names(moduleLabels)
  # module_labels = as.data.frame(moduleLabels)
  # module_labels$gene = rownames(module_labels)
  # genes_module = merge(annot,module_labels,by = "gene")
  # write.csv(genes_module,file = "genes_module.csv", row.names = FALSE)
  # genes_module = read.csv("genes_module.csv")
  # 
  # # # Annotate the Ensembl id of WGCNA result, with Entrez IDs, via function grex
  # annotComplete = grex(cleanid(genes_module$id))
  # genes_module$id = cleanid(genes_module$id)
  # annotComplete = merge(annotComplete,genes_module[,-1],by.x = "ensembl_id",by.y = "id")
  # 
  # write.csv(annotComplete, "topGO_moduleGeneAnnotated.csv", row.names = FALSE)
  annotComplete = read.csv("topGO_moduleGeneAnnotated.csv")
  # Connect to the ENSEMBL_MART_ENSEMBL BioMart database to query GO IDs for our Entrez IDs
  genes_bg = annotComplete$entrez_id
  tot_background = length(genes_bg);
  
  db = useMart("ENSEMBL_MART_ENSEMBL", 
               dataset = "hsapiens_gene_ensembl",
               host = "https://www.ensembl.org");
  
  go_ids = getBM(attributes = c("go_id",
                                "entrezgene_id", "namespace_1003"),
                 filters = "entrezgene_id", values = genes_bg, mart = db)
  
  # Define the WGCNA module of interest then set up named factors for genes
  # located within and outside of the module of interest
  modInt = as.factor(annotComplete$module)
  annotSplit = split(annotComplete,modInt)
  candidate_list = annotSplit[[interest_module+1]]$entrez_id
  tot_candidate = length(candidate_list);
  
  keep = candidate_list %in%  go_ids[,2];
  keep = which(keep == TRUE);
  candidate_list = candidate_list[keep];
  
  geneList = factor(as.integer(genes_bg %in% candidate_list))
  names(geneList) = genes_bg
  
  # Build a topGO data object with the base function new, then run GO analysis.                  
  gene2GO = unstack(go_ids[,c(1,2)])
  
  GOdata_interest_module = new("topGOdata",
                               ontology = c("BP"), allGenes = geneList,
                               annot = annFUN.gene2GO,
                               gene2GO = gene2GO,  nodeSize = 5)
  
  save(GOdata_interest_module,file = 
         paste0("GOdata_module_",interest_module,".RData"))
  
  # load(paste0("GOdata_module_",interest_module,".RData"))
  GOdata = GOdata_interest_module
  allGO = usedGO(GOdata)
  
  # Test significance of GO terms using the function runTest.
  # algorithm and statistic specify the method for dealing with the GO graph 
  # structure and the test statistic, respectively.
  
  # sigTest = runTest(GOdata,
  #                   algorithm = "elim", statistic = "fisher")
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  # resultks <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  
  
  # Generate a summary table of results obtained from topGO enrichment analysis with the function
  # GenTable.
  
  all_res = GenTable(GOdata, weightFisher = resultFisher, orderBy = "weightFisher", 
                     topNodes = length(allGO))
  
  
  # goID <- all_res_ks[1, "GO.ID"]
  # print(showGroupDensity(GOdata, goID, ranks = TRUE))
  # 
  # goID <- all_res_ks[10, "GO.ID"]
  # sg <- sigGenes(GOdata)
  # str(sg)
  # numSigGenes(GOdata)
  
  
  # Calculate odds ratios 
  all_res$OR = log2((all_res$Significant/numSigGenes(GOdata))/
                      (all_res$Annotated/tot_background))
  
  
  # visualization of top 20 GO entities 
  # Generate a summary figure of topGO results via the ggplot2 R package
  #This sets the threshold p<0.05 and #.annotated.background.genes >= 30. 
  # GO_bar = all_res %>%
  #   filter(weightFisher < 0.05) %>% filter(Annotated >= 30);
  GO_bar = all_res 
  
  # GO_bar = GO_bar %>%
  #   select(Term, weightFisher, OR);
  
  GO_bar = transform(GO_bar,weightFisher = as.numeric(weightFisher));
  GO_bar$adjusted_Fisher = p.adjust(GO_bar$weightFisher)
  #This selects the top 20 most significant(i,e,,orcler by ascending p-values) GO terms for presentation. 
  GO_bar = GO_bar %>%
    arrange(adjusted_Fisher) %>% slice(1:20);
  
  
  ggplot(GO_bar,
         aes(x = reorder(Term, -adjusted_Fisher), y = OR, fill = adjusted_Fisher))+
    geom_bar(stat = "identity", width = 0.9, color = "black") +
    coord_flip() +
    scale_fill_gradient(low="#feff2b",high="#fe0100")+
    ylim(0,8)+
    ggtitle(paste("Enriched GO Biological Processes, classic fisher for module",interest_module))+
    labs(X = NULL,
         y = "Odds Ratio", fill = "p.value")+theme_bw()+
    theme(plot.title.position = "plot")+
    theme(plot.title = element_text(size = 10))+
    theme(axis.title.x = element_text(size = 7,face = "bold"),
          axis.title.y = element_text(size = 7,face = "bold"))+
    theme(axis.text.x = element_text(size = 7,face = "bold"),
          axis.text.y = element_text(size = 7,face = "bold"))+
    theme(legend.position = "right");
  ggsave(paste0("top 20 enriched GO for module",interest_module,"classic fisher",".pdf"))
  
  
  
  sigTest = runTest(GOdata,
                    algorithm = "elim", statistic = "fisher")
  
  # Generate a summary table of results obtained from topGO enrichment analysis with the function
  # GenTable.
  
  all_res = GenTable(GOdata, weightFisher = sigTest, orderBy = "weightFisher", 
                     topNodes = length(allGO))
  
  # Calculate odds ratios 
  all_res$OR = log2((all_res$Significant/numSigGenes(GOdata))/
                      (all_res$Annotated/tot_background))
  
  
  # visualization of top 10 GO entities 
  # Generate a summary figure of topGO results via the ggplot2 R package
  #This sets the threshold p<0.05 and #.annotated.background.genes >= 30. 
  
  GO_bar = all_res 
  
  GO_bar = transform(GO_bar,weightFisher = as.numeric(weightFisher));
  GO_bar$adjusted_Fisher = p.adjust(GO_bar$weightFisher)
  
  #This selects the top 10 most significant(i,e,,orcler by ascending p-values) GO terms for presentation. 
  GO_bar = GO_bar %>%
    arrange(adjusted_Fisher) %>% slice(1:20);
  
  
  ggplot(GO_bar,
         aes(x = reorder(Term, -adjusted_Fisher), y = OR, fill = adjusted_Fisher))+
    geom_bar(stat = "identity", width = 0.9, color = "black") +
    coord_flip() +
    scale_fill_gradient(low="#feff2b",high="#fe0100")+
    ylim(0,8)+
    ggtitle(paste("Enriched GO Biological Processes,elim fisher for module",interest_module))+
    labs(X = NULL,
         y = "Odds Ratio", fill = "p.value")+theme_bw()+
    theme(plot.title.position = "plot")+
    theme(plot.title = element_text(size = 9))+
    theme(axis.title.x = element_text(size = 5,face = "bold"),
          axis.title.y = element_text(size = 5,face = "bold"))+
    theme(axis.text.x = element_text(size = 5,face = "bold"),
          axis.text.y = element_text(size = 5,face = "bold"))+
    theme(legend.position = "right");
  ggsave(paste0("top 20 enriched GO for module",interest_module,"elim fisher.pdf"))
  
  output_file_names = c(paste0("GOdata_module_",interest_module,".RData"),
                        paste0("top 20 enriched GO for module",interest_module,"classic fisher",".pdf"),
                        paste0("top 20 enriched GO for module",interest_module,"elim fisher",".pdf"))
  return(output_file_names)
}
output_file_names = GO_enrichment_module(7)
output_file_names = GO_enrichment_module(11)



##########
# sanity check, pick the top 10 co-expressed genes with OBP2A, plot the pearson 
# correlation between these genes, check whether their correlation scores are high
sanity_check <- function(workdir = "/home/lyang73/data/breast_cancer",
                         input_file = "OBP2A_score.RData",interest_gene = "OBP2A",
                         counts_file = "00_TCGA breast Cancer_norm_all_ExpectedCnt.csv" )
{
  library(dplyr)
  library(limma) 
  library(edgeR)
  library(ggplot2)
  library(corrplot)
  setwd(workdir)
  load(input_file)
  
  # use topological overlap matrix as weights
  edgeList_OBP = data.frame(Weight = interest_gene_score, genes = names(interest_gene_score))
  
  
  candidates = edgeList_OBP[edgeList_OBP$genes %in% c("LCN2","OBP2A","OBP2B","OR51E1","OR2B6","OR51E2"),]
  write.csv(candidates,file =paste0(interest_gene,"_candidate_genes_scores.csv"),
            row.names = FALSE)
  
  edgeList_OBP = edgeList_OBP[sort(edgeList_OBP$Weight, decreasing = TRUE, 
                                   index.return=TRUE)$ix[1:30],]  
  write.csv(edgeList_OBP,file =paste0(interest_gene,"_top30_correlated_genes.csv"),
            row.names = FALSE)
  
  
  edgeList_OBP = edgeList_OBP[sort(edgeList_OBP$Weight, decreasing = TRUE, 
                                   index.return=TRUE)$ix[1:10],]
  
  # 
  # edgeList_OBP$TargetName[edgeList_OBP$SourceName != "OBP2A"] = edgeList_OBP[edgeList_OBP$SourceName != "OBP2A","SourceName"]
  # edgeList_OBP$SourceName = "OBP2A"
  interest_genes = edgeList_OBP$genes
  
  norm_exprFinal = read.csv(counts_file)
  
  dim(norm_exprFinal)
  norm_expr_correlation = norm_exprFinal[norm_exprFinal$gene %in% 
                                           interest_genes,]
  
  dim(norm_expr_correlation)
  
  expr_correlation = norm_expr_correlation
  expr_statist = expr_correlation[,-c(1:3)];
  rownames(expr_statist) = expr_correlation$gene
  
  expr = as.data.frame(t(expr_statist)) 
  expr_cancer = expr[!startsWith(rownames(expr), "GTEX"),]
  expr_normal = expr[startsWith(rownames(expr), "GTEX"),]
  exprSet = expr_cancer
  # exprSet = expr_cancer[,-ncol(expr_cancer)]
  # exprSet = sapply(exprSet,as.numeric)
  M <- stats::cor(x = exprSet) 
  
  pdf(file = paste0("sanity check of correlation of ",interest_gene," in breast cancer.pdf"))
  
  corrplot(M,order = "AOE",addCoef.col = "white", mar=c(0,0,1,0),type = "lower",
           title = paste("Pearson correlation between genes strongly co-expressed with",interest_gene))
  
  dev.off()
  
  if(file.exists("correlation plot of co-expression scores of genes in breast cancer.pdf"))
  {
    output_file_names = c( paste0(interest_gene,"_top30_correlated_genes.csv"),
                           paste0("sanity check of correlation of ",interest_gene," in breast cancer.pdf")
    )
    
  }else
  {
    library(corrplot)
    for (gene in c("OBP2A","OBP2B","LCN2")) 
    {
      
      load(paste0(gene,"_score.RData"))
      if(gene == "OBP2A")
      {
        score_all = as.data.frame(interest_gene_score)
        colnames(score_all) = gene
      }else
      {
        temp = as.data.frame(interest_gene_score)
        colnames(temp) = gene
        score_all = cbind(score_all,temp)
      }
      rm(interest_gene_score)
      
    }
    
    M <- stats::cor(x = score_all) 
    
    
    
    pdf(file = "correlation plot of co-expression scores of genes in breast cancer.pdf")
    
    corrplot(M,order = "AOE",addCoef.col = "white", mar=c(0,0,1,0),type = "lower",
             title = paste("correlation between co-expression scores of interest genes in breast cancer"))
    
    
    dev.off()
    output_file_names = c( paste0(interest_gene,"_top30_correlated_genes.csv"),
                           paste0("sanity check of correlation of ",interest_gene," in breast cancer.pdf"),
                           "correlation plot of co-expression scores of genes in breast cancer.pdf" )
    
  }
  
  return(output_file_names)
}

output_file_names = sanity_check(workdir = "D:/major_project/breast_cancer",
                                 input_file = "OBP2A_score.RData","OBP2A")
output_file_names = sanity_check(workdir = "D:/major_project/breast_cancer",
                                 input_file = "OBP2B_score.RData","OBP2B")
output_file_names = sanity_check(workdir = "D:/major_project/breast_cancer",
                                 input_file = "LCN2_score.RData","LCN2")

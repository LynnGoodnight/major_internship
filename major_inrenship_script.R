
##############################
# SRAT paper tutoral
# download data, normalize, do statistical test and pearson correlation plot
# install.packages("BiocManager")
# 
# install.packages("UCSCXenaTools")
# install.packages("R.utils")

initialize_parameter <- function()
{
  paraCohort = c("TCGA breast Cancer","TCGA Ovarian Cancer",
                 "TCGA Uterine Carcinosarcoma","TCGA Testicular Cancer",
                 "TCGA Stomach Cancer","TCGA Melanoma",
                 "TCGA Prostate Cancer","TCGA Lung Cancer","TCGA Glioblastoma",
                 "TCGA Lower Grade Glioma","TCGA Liver Cancer",
                 "TCGA Esophageal Cancer","TCGA Pancreatic Cancer",
                 "TCGA Colon Cancer","TCGA Adrenocortical Cancer");
  #Selecting the Colon Cancer cohort.
  paraDatasets = c("TCGA.BRCA.sampleMap/BRCA_clinicalMatrix",
                   "TCGA.OV.sampleMap/OV_clinicalMatrix",
                   "TCGA.UCS.sampleMap/UCS_clinicalMatrix",
                   "TCGA.TGCT.sampleMap/TGCT_clinicalMatrix",
                   "TCGA.STAD.sampleMap/STAD_clinicalMatrix",
                   "TCGA.SKCM.sampleMap/SKCM_clinicalMatrix",
                   "TCGA.PRAD.sampleMap/PRAD_clinicalMatrix",
                   "TCGA.LUNG.sampleMap/LUNG_clinicalMatrix",
                   "TCGA.GBM.sampleMap/GBM_clinicalMatrix",
                   "TCGA.LGG.sampleMap/LGG_clinicalMatrix",
                   "TCGA.LIHC.sampleMap/LIHC_clinicalMatrix",
                   "TCGA.ESCA.sampleMap/ESCA_clinicalMatrix",
                   "TCGA.PAAD.sampleMap/PAAD_clinicalMatrix",
                   "TCGA.COAD.sampleMap/COAD_clinicalMatrix",
                   "TCGA.ACC.sampleMap/ACC_clinicalMatrix");
  paraPrimarySiteTCGA = c("Breast","Ovary","Uterus","Testis","Stomach","Skin",
                          "Prostate","Lung","Brain","Brain","Liver","Esophagus",
                          "Pancreas","Colon","Adrenal gland");
  #Setting "Colon" as the primary site of interest.
  
  #Selecting the Colon Cancer clinical matrix.
  paraPrimarySiteGTEx = c("Breast","Ovary","Uterus","Testis","Stomach","Skin",
                          "Prostate","Lung","Brain","Brain","Liver","Esophagus",
                          "Pancreas","Colon","Adrenal Gland");#Setting "Colon" as the primary site of interest. 
  paraPrimaryTissueGTEx = c("^Breast","Ovary","Uterus","Testis","Stomach",
                            "^Skin","Prostate","Lung",
                            "Brain - Hypothalamus","Brain - Hypothalamus",
                            "Liver","^Esophagus","Pancreas","^Colon",
                            "Adrenal Gland");#Setting "Colon" as the primary tissue of interest.
  # summary_table = data.frame("Primary Site" = paraPrimarySiteGTEx,"cancer database" = paraCohort,
  #            "Primary Tissue GTEx" = paraPrimaryTissueGTEx)
  # write.csv(summary_table,"summary table.csv")
  newList <- list("paraCohort" = paraCohort, "paraDatasets" = paraDatasets,
                  "paraPrimarySiteTCGA" = paraPrimarySiteTCGA,
                  "paraPrimarySiteGTEx" = paraPrimarySiteGTEx,
                  "paraPrimaryTissueGTEx" = paraPrimaryTissueGTEx)
  return(newList)
  
}
parameter_ist = initialize_parameter()
paraCohort = parameter_ist$paraCohort
paraDatasets = parameter_ist$paraDatasets
paraPrimarySiteTCGA = parameter_ist$paraPrimarySiteTCGA
paraPrimarySiteGTEx = parameter_ist$paraPrimarySiteGTEx
paraPrimaryTissueGTEx = parameter_ist$paraPrimaryTissueGTEx



retrieve_data <- function(paraCohort, paraDatasets, paraPrimarySiteGTEx,
                          paraPrimaryTissueGTEx,paraPrimarySiteTCGA) 
{
  if(file.exists(paste0("00_",gsub("Cancer","",paraCohort),"_all_ExpectedCnt.csv")))
  {
    # expr_correlation = read.csv(paste0("00_",paraCohort,"selected_ExpectedCnt.csv"),
    #                             check.names = FALSE)
    # expr_correlation = expr_correlation[,-1]
    clinFinal= read.csv(paste0("00_",gsub("Cancer","",paraCohort),"_ClinTraits.csv")) 
    clinFinal = clinFinal[,-1]
    exprFinal = read.csv(paste0("00_",gsub("Cancer","",paraCohort),"_all_ExpectedCnt.csv"),
                         check.names = FALSE)
    
  }else{
    library(UCSCXenaTools)
    library(data.table)
    library(R.utils)
    library(dplyr)
    
    #4. Generate a record of datasets hosted on UCSC Xena Data Hubs via the UCSCXenaTools
    # R package
    data(XenaData);
    # write.csv(XenaData, "00_tblxenaHubInfo.csv")
    
    # 5. Retrieve gene expression, clinical, survival, and phenotype data from the UCSC Xena platform via
    # the UCSCXenaTools R package (Wang and Liu, 2019). For network analysis of TCGA colon cancer
    # gene expression, we selected and downloaded:
    #a. download The TcgaTargetGtex_gene_expected_count dataset for the TCGA TARGET GTEx cohort
    # from host toilHub
    GeneExpectedCnt_toil = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
      XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
      XenaFilter(filterDatasets = "TcgaTargetGtex_gene_expected_count");
    XenaQuery(GeneExpectedCnt_toil) %>%
      XenaDownload(destdir = "./")
    
    # b.The COAD_clinicalMatrix dataset for the TCGA Colon Cancer cohort from host
    # tcgaHub by setting the filters paraCohort ="TCGA Colon Cancer" and paraDatasets =
    #   "TCGA.COAD.sampleMap/COAD_clinicalMatrix".
    Clin_TCGA = XenaGenerate(subset = XenaHostNames == "tcgaHub") %>%
      XenaFilter(filterCohorts = paraCohort) %>%
      XenaFilter(filterDatasets = paraDatasets);
    XenaQuery(Clin_TCGA) %>%
      XenaDownload(destdir = "./")
    
    
    # c.The TCGA_survival_data dataset for the TCGA TARGET GTEx cohort from host toilHub.
    Surv_TCGA = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
      XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
      XenaFilter(filterDatasets = "TCGA_survival_data");
    XenaQuery(Surv_TCGA) %>%
      XenaDownload(destdir = "./")
    
    # d.The TcgaTargetGTEX_phenotype dataset for the TCGA TARGET GTEx cohort from host
    # toilHub.
    
    Pheno_GTEx = XenaGenerate(subset = XenaHostNames == "toilHub") %>%
      XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
      XenaFilter(filterDatasets = "TcgaTargetGTEX_phenotype");
    XenaQuery(Pheno_GTEx) %>%
      XenaDownload(destdir = "./")
    
    # 6. Subset the gene expression matrix to include only observations of desired tissue type(s).
    #a. The Genotype-Tissue Expression (GTEx) project provides gene expression data from healthy,
    # cancer-free individuals. For differential gene expression analysis, we selected GTEx normal
    # colon tissue samples by setting the filters paraStudy = "GTEX", paraPrimarySiteGTEx =
    #   "Colon" and paraPrimaryTissueGTEx = "^Colon" to the TcgaTargetGTEX_phenotype dataset. Table 3 can be used to set values for the arguments paraPrimarySiteGTEx
    # and paraPrimaryTissueGTEx.
    
    filterGTEx01 = fread("TcgaTargetGTEX_phenotype.txt.gz");
    names(filterGTEx01)=gsub("\\_","",  names(filterGTEx01));
    
    paraStudy = "GTEX";#Setting "GTEx" as the study of interest.
    
    filterGTEx02=subset(filterGTEx01,
                        study == paraStudy &
                          primarysite == paraPrimarySiteGTEx &
                          grepl(paraPrimaryTissueGTEx,filterGTEx01$`primary disease or tissue` ))
    
    # b.The Cancer Genome Atlas (TCGA) program provides gene expression data from
    # primary tumors. For differential gene expression analysis, we selected TCGA colon cancer
    # primary tumor samples by setting the filters paraSampleType = "Primary Tumor",
    # paraPrimarySiteTCGA ="Colon", and paraHistologicalType = "Colon Adenocarcinoma". Table 4 can be used to set values for the arguments paraPrimarySiteTCGA
    # and paraHistologicalType.
    
    filterTCGA01 = fread(paraDatasets);
    names(filterTCGA01)=gsub("\\_","", names(filterTCGA01));
    # paraSampleType = "Primary Tumor";#Setting "Primary Tumor" as the sample type of interest. 
    # paraHistologicalType = "Infiltrating Ductal Carcinoma";#Setting "Colon Adenocarcinoma" as the histological type of intere st.
    
    # filterTCGA02 = subset(filterTCGA01,
    #                       sampletype == paraSampleType &
    #                         primarysite == paraPrimarySiteTCGA &
    #                         grepl(paraHistologicalType,filterTCGA01$histologicaltype))
    
    filterTCGA02 = subset(filterTCGA01, primarysite == paraPrimarySiteTCGA)
    
    # c. The TcgaTargetGtex_gene_expected_count dataset from the toilHub data hub combines RNA-Seq data 
    # from TCGA and GTEx by uniformly realigning reads to the hg38 genome
    # and re-calling expressions using RSEM and Kallisto methods (Vivian et al., 2017). To compare
    # gene expression between GTEx normal and TCGA tumor for network analyses, subset the
    # TcgaTargetGtex_gene_expected_count dataset via lists generated during step 6.a and
    # step 6.b
    filterExpr = c(filterGTEx02$sample,filterTCGA02$sampleID,"sample");
    
    ExprSubsetBySamp = fread("TcgaTargetGtex_gene_expected_count.gz",
                             select = filterExpr)
    
    # filterExpr = c(filterGTEx02$sample,"sample");
    # 
    # ExprSubsetBySamp = fread("TcgaTargetGtex_gene_expected_count.gz",
    #                          select = filterExpr)
    
    # 7. Subset the gene expression matrix to include only protein-coding genes. This can be achieved by
    # utilizing the zz_gene.protein.coding.csv saved during step 1
    probemap = fread("zz_gencode.v23.annotation.csv",select = c(1,2));
    exprALL = merge(probemap,ExprSubsetBySamp,by.x = "id", by.y = "sample");
    genesPC = fread("zz_gene.protein.coding.csv");
    exprPC = subset(exprALL, gene %in%  genesPC$Gene_Symbol);
    #Remove duplicate gene symbols.
    exprFinal = exprPC[!(duplicated(exprPC$gene)|
                           duplicated(exprPC$gene,fromLast = TRUE)),  ]
    
    # expr_correlation = exprFinal[exprFinal$gene %in% c("LCN2","OBP2A","OBP2B","OR51E1",
    #                                                    "OR2B6","OR51E2"),]
    
    # 8. The gene expression matrix can now be saved for downstream analyses.
    # write.csv(expr_correlation, paste0("00_",paraCohort,"selected_ExpectedCnt.csv"))
    write.csv(exprFinal, paste0("00_",paraCohort,"_all_ExpectedCnt.csv"))
    
    # 9. The COAD_clinicalMatrix dataset contains 133 administrative and phenotypic annotations
    # (see the Methods S1 file ????R Markdown Code Script for UCSCXenaTool???? for a full list), ranging
    # from sample IDs to pathologic stage to KRAS mutation codon. Keep only the variable(s)
    # of interest. For example, to identify potential biomarkers for lymphatic invasion during network
    # analysis of TCGA colon cancer gene expression, we retained the following variables:
    names(filterTCGA02)
    #Keep variable "Lymphatic Invasion".
    if(paraPrimarySiteGTEx == "Skin"){
      varClinKeep = c("sampleID","sampletype")
    }else{
      varClinKeep = c("sampleID","sampletype","histologicaltype");
    }
    
    # varClinKeep = c("sampleID","sampletype","histologicaltype",
    #                 "AJCCStagenature2012");
    clinDF01 = as.data.frame(do.call(cbind,filterTCGA02));
    clinFinal = clinDF01[varClinKeep];
    
    # clinFinal$clinicalstage = clinFinal$AJCCStagenature2012
    # clinFinal$clinicalstage[grepl("^Stage IV", clinFinal$clinicalstage)] <- 4;
    # clinFinal$clinicalstage[grepl("^Stage III", clinFinal$clinicalstage)] <- 3;
    # clinFinal$clinicalstage[grepl("^Stage II", clinFinal$clinicalstage)] <- 2;
    # clinFinal$clinicalstage[grepl("^Stage I", clinFinal$clinicalstage)] <- 1
    
    #Identify observations/samples with no values assigned to the kept variables. 
    colSums(clinFinal =="");
    colSums(is.na(clinFinal));
    
    #Replace "no values" with "NA"
    
    NA -> clinFinal[clinFinal  == ""];
    
    colSums(is.na(clinFinal));
    
    
    clinFinal = clinFinal[!is.na(clinFinal$sampletype),]
    #Verify that the count of 1/0 mirrors previous count of YES/NO. 
    table(clinFinal$sampletype);
    
    
    # varClinKeep = c("sampleID",  "clinicalstage")
    # 
    
    
    # 10. The phenotype annotation matrix can now be saved for downstream analyses.
    
    write.csv(clinFinal,paste0("00_",paraCohort,"_ClinTraits.csv"))
    
  }
  newList <- list("clinFinal" = clinFinal, "exprFinal" = exprFinal)
  return(newList)
}
norm_library_size <- function(exprFinal,paraCohort)
{
  if(file.exists(paste0("00_",gsub("Cancer","",paraCohort),"_norm_all_ExpectedCnt.csv")))
  {
    norm_y = read.csv(paste0("00_",gsub("Cancer","",paraCohort),"_norm_all_ExpectedCnt.csv"))
  }else
  {
    expr_all = as.data.frame(exprFinal)
    expr_all = expr_all[,-1:-3]
    library(edgeR)
    library(reshape2)
    library(ggpubr)
    y <- DGEList(counts=expr_all)
    y <- calcNormFactors(y)
    norm_y <- cpm(y, log = TRUE, normalized.lib.sizes=TRUE)
    # test  <- cpm(y, log = FALSE, normalized.lib.sizes=FALSE)
    norm_y = as.data.frame(norm_y)
    norm_y = cbind(exprFinal[,c("id","gene")],norm_y)
    write.csv(norm_y, paste0("00_",
                             gsub("Cancer","",paraCohort),"_norm_all_ExpectedCnt.csv"),
              row.names = FALSE)
    
  }
    
  
  return(norm_y)
  # expr_all$sampleid = colnames(expr_all)[3:10]
  # install.packages("reshape2")
  
  # norm_y = melt(norm_y, measure.vars = 1:11)
  # colnames(norm_y)[1] = "sampleid"
  # ggboxplot(norm_y, x = "sampleid", y = "value") 
  
}
statistical_compare <- function(expr_correlation,clinFinal,paraPrimarySiteGTEx,
                                paraCohort,selected_gene ="OBP2A")
{
  
  library(ggpubr)
  # load(".RData")
  
  # exprFinal = read.csv("00_ExpectedCnt.csv");
  expr_statist = expr_correlation[,-c(1:2)];
  
  
  
  
  expr_statist = as.data.frame(t(expr_statist))
  colnames(expr_statist) = expr_correlation$gene;
  expr_statist$sampleID = rownames(expr_statist)
  
  # table(clinFinal$sampletype)
  
  expr_statist = merge(x=expr_statist,y=clinFinal,by="sampleID",all.x = TRUE)
  expr_statist[startsWith(expr_statist$sampleID, "GTEX"),"sampletype" ] = "normal"
  # types = c(unique(expr_statist$sampletype))
  summary_group = table(expr_statist$sampletype)
  types = names(summary_group[table(expr_statist$sampletype) > 5]) 
  if(length(types) <= 2){
    my_comparisons <- list(types)
  }else{
    my_comparisons <- combn(types, 2,simplify = FALSE)
  }
  
  
  # Visualize: Specify the comparisons you want
  # function for number of observations 
  give.n <- function(x){
    return(c(y = median(x)*1.05, label = length(x))) 
    # experiment with the multiplier to find the perfect position
  }  
  ggboxplot(expr_statist, x = "sampletype", y = selected_gene,
            color = "sampletype", palette = "jco")+ 
    stat_compare_means(comparisons = my_comparisons, label = "p.signif",
                       method = "t.test", p.adjust.method = "fdr") + # Add pairwise comparisons p-value
    ggtitle(paste(selected_gene,"expression in",paraPrimarySiteGTEx,"tissue")) +
    xlab("sample type") + ylab(paste(selected_gene,"log2(RSEM norm coun+1)") )+
    labs(color = "sample type")+
    stat_summary(fun.data = give.n, geom = "text", fun = median)
  # annotate("text",
  #          x = 1:length(table(expr_statist$sampletype)),
  #          y = aggregate(OBP2A ~ sampletype, expr_statist, median)[ , 2],
  #          label = table(expr_statist$sampletype),
  #          col = "red",
  #          vjust = 0)
  # ns=non significant 
  if(paraPrimarySiteGTEx == "Brain")
  {
    ggsave(paste(selected_gene,"expression in",paraCohort,"tissue",".png"), width = 20, height = 20, units = "cm")
    print(paste(selected_gene,"expression in",paraCohort,"tissue",".png"))
  }else{
    ggsave(paste(selected_gene,"expression in",paraPrimarySiteGTEx,"tissue",".png"), width = 20, height = 20, units = "cm")
    print(paste(selected_gene,"expression in",paraPrimarySiteGTEx,"tissue",".png"))
  }
  
}
correlation_between_genes <- function(expr_correlation,clinFinal,
                                      paraPrimarySiteGTEx,paraCohort) 
{
  library(dplyr)
  library(limma) 
  library(edgeR)
  library(ggplot2)
  library(corrplot)
  
  expr_statist = expr_correlation[,-c(1:2)];
  expr_statist = as.data.frame(t(expr_statist))
  colnames(expr_statist) = expr_correlation$gene;
  expr_statist$sampleID = rownames(expr_statist)
  
  
  # exprBT = round(((2^expr_statist[,1:6])-1),0)
  # write.csv(exprBT, "01_ExpectedCntBT.csv")    
  
  # ʹ??5???????ı???����????????ͼ
  expr = merge(x=expr_statist,y=clinFinal,by="sampleID",all.x = TRUE)
  
  
  expr_cancer = expr[!startsWith(expr$sampleID, "GTEX"),]
  expr_normal = expr[startsWith(expr$sampleID, "GTEX"),]
  
  exprSet = expr_cancer[,2:8]
  # ??????????
  M <- cor(exprSet)  #M??һ?u
  
  # resulting value will be NA whenever one of its contributing observations is NA.
  png(height=1800, width=1800, 
      filename = paste0("expression correlation in ",paraCohort,".png"), type = "cairo",
      pointsize = 40)   
  
  corrplot(M,order = "AOE",addCoef.col = "white", mar=c(0,0,1,0),type = "full",
           title = paste("expression correlation between genes in",paraCohort))
  
  dev.off()
  print(paste0("expression correlation in ",paraCohort,".png"))
  
  ggplot(data=exprSet, aes(x=OBP2A, y=OBP2B)) + geom_point() +
    ggtitle(paste("correlation dot plot in",paraCohort))
  
  ggsave(paste("correlation dot plot in",paraCohort,".png"), 
         width = 20, height = 20, units = "cm")
  print(paste("correlation dot plot in",paraCohort,".png"))
  ##############  
  
  M <- cor(expr_normal[,2:8]) 
  png(height=1800, width=1800,
      file=paste("expression correlation in normal",paraPrimarySiteGTEx,"tissue.png"),
      type = "cairo",pointsize = 40)
  corrplot(M,order = "AOE",addCoef.col = "white", mar=c(0,0,1,0),type = "full",
           title = paste("expression correlation in normal",paraPrimarySiteGTEx,"tissue"))
  
  dev.off()
  print(paste("expression correlation in normal",paraPrimarySiteGTEx,"tissue.png"))
  
  ggplot(data=expr_normal[,2:7], aes(x=OBP2A, y=OBP2B)) + geom_point()+
    ggtitle(paste("correlation dot plot in normal",paraPrimarySiteGTEx,"tissue"))
  ggsave(paste("correlation dot plot in normal",paraPrimarySiteGTEx,"tissue.png"), 
         width = 20, height = 20, units = "cm")
  print(paste("correlation dot plot in normal",paraPrimarySiteGTEx,"tissue.png"))
  
}

setwd("D:/major_project")
n = 1
while (n < 2) {
  # print(n)
  # 
  # data = retrieve_data(paraCohort[n],paraDatasets[n], paraPrimarySiteGTEx[n],
  #                      paraPrimaryTissueGTEx[n],paraPrimarySiteTCGA[n])
  # clinFinal = data$clinFinal
  # exprFinal = data$exprFinal
  # print(table(clinFinal$sampletype))
  # norm_exprFinal = norm_library_size(exprFinal,paraCohort[n])
  # norm_expr_correlation = norm_exprFinal[norm_exprFinal$gene %in%
  #                                          c("LCN2","OBP2A","OBP2B","OR51E1","OR2B6","OR51E2","LCN1"),]
  # write.csv(norm_expr_correlation, paste0("00_",gsub("Cancer","",paraCohort[n]),"_norm_selected_count.csv"))

  
  correlation_between_genes(norm_expr_correlation,clinFinal,paraPrimarySiteGTEx[n],
                            paraCohort[n])
  # my_comparisons <- list( c("normal", "Primary Tumor"), c("normal", "Solid Tissue Normal"), 
  #                         c("normal","Metastatic"),
  #                         c("Primary Tumor","Solid Tissue Normal"),
  #                         c("Primary Tumor","Metastatic"),
  #                         c("Solid Tissue Normal","Metastatic"))
  
  statistical_compare(norm_expr_correlation,clinFinal,paraPrimarySiteGTEx[n],
                      paraCohort[n],"OBP2A")
  statistical_compare(norm_expr_correlation,clinFinal,paraPrimarySiteGTEx[n],
                      paraCohort[n],"OBP2B")
  
  n = n+1
}


##################################
# WGCNA: module construction
# identify which modules OBP2A/B are assigned
# identify gene set(s) highly co-expressed with OBP2A/B
# using WGCNA tutorial as reference to write the code and scripts
#################
# WGCNA preprocess data

preprocess_explore_data <- function(workdir = "/home/lyang73/data/breast_normal",
                                    state = "normal",
                                    input_file = "00_TCGA breast _norm_all_ExpectedCnt.csv")
{
  # dir.create("/home/lyang73/data/breast_cancer")
  # data input and pre-process
  setwd(workdir)
  
  # setwd("/home/lyang73/data")
  library(GO.db)
  # Load the WGCNA package
  library(WGCNA)
  
  # check for genes and samples with too many missing values
  # If the last statement returns TRUE, all genes have passed the cuts. 
  # 
  
  exprFinal = read.csv(input_file,
                       check.names = FALSE)
  exprdata = as.data.frame(t(exprFinal[,-1:-3]))
  colnames(exprdata) = exprFinal[,"gene"]
  
  expr_cancer = exprdata[!startsWith(rownames(exprdata), "GTEX"),]
  expr_normal = exprdata[startsWith(rownames(exprdata), "GTEX"),]
  if(state == "cancer")
  {
    exprdata = expr_cancer
  }else if(state == "normal")
  {
    exprdata = expr_normal
  }
  
  # Remove genes that are lowly-expressed and/or genes with low variation between samples.
  # For WGCNA on TCGA gene expression data, we removed genes that has normalized 
  # expression value< O in all of our samples
  
  
  #Add a new variable "Count" that counts the number of samples that have normalized expression value <0 for gene x 
  count = colSums(exprdata < 0);
  table(count);
  
  
  
  
  
  # This function checks data for missing entries, entries with weights below a
  # threshold, and zero-variance genes, and returns a list of samples and genes
  # that pass criteria on maximum number of missing or low weight values. If 
  # necessary, the filtering is iterated.
  gsg = goodSamplesGenes(exprdata, verbose = 3);
  # BiocManager::install("genefilter")
  # library(genefilter)
  # test = t(exprdata)
  # selProbes <- genefilter(test, filterfun(pOverA(0.20, log2(100)),
  #                                         function(x) (IQR(x) > 0.25)))
  # eset <- ALL[selProbes, ]
  
  
  gsg$allOK
  # If not, we remove the offending genes and samples
  # from the data:
  if (!gsg$allOK)
  { # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(exprdata)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(exprdata)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    exprdata = exprdata[gsg$goodSamples, gsg$goodGenes]
  }
  
  
  # cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious
  # outliers.
  sampleTree = hclust(dist(exprdata), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  # sizeGrWindow(12,9)
  pdf(file = "WGCNA_sample_Clustering.pdf", width = 20, height = 12);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  abline(h = 75, col = "red");
  # Determine cluster under the line
  clust = cutreeStatic(sampleTree, cutHeight = 75, minSize = 10)
  
  dev.off()
  
  table(clust)
  # It appears there is one outlier (sample F2_221, see Fig. 1). One can remove it by hand, or use an automatic approach.
  # Choose a height cut that will remove the offending sample, say 15 (the red line in the plot), and use a branch cut at
  # that height.
  # Plot a line to show the cut
  
  
  
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(exprdata, powerVector = powers, verbose = 5)
  # Plot the results:
  pdf(file = "WGCNA_Analysis_of_network_topology.pdf", width = 18, height = 12);
  
  # sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  abline(h=0.8,col="red") # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  save(exprdata, file = "WGCNA_01-dataInput_mydata.RData")
  return(exprdata)
}
# exprdata = preprocess_explore_data("D:/major_project", "00_TCGA breast Cancer_norm_all_ExpectedCnt.csv")

# exprdata = preprocess_explore_data("/home/lyang73/data/breast_cancer", "breast",
#                                    "00_TCGA breast _norm_all_ExpectedCnt.csv")

exprdata = preprocess_explore_data("/home/lyang73/data/breast_normal", "normal",
                                   "00_TCGA breast _norm_all_ExpectedCnt.csv")

##########
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

library(WGCNA)
for (val in seq(4,10))
{
  print(val)
  enableWGCNAThreads()
  net = blockwiseModules(exprdata, power = val,
                         TOMType = "unsigned", minModuleSize = 20,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE, maxBlockSize = 20000,deepSplit = 2,
                         saveTOMFileBase = "mydata",
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
  geneTree = net$dendrograms[[1]];
  if(val >= 10)
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

###########
# WGCNA Step-by-step network construction and module detection
load(file = "WGCNA_01-dataInput_mydata.RData")

setwd("/home/lyang73/data"); 
# Load the WGCNA package
library(WGCNA)


# calculate the adjacencies, Co-expression similarity
softPower = 10;
adjacency = adjacency(exprdata, power = softPower);

# To minimize effects of noise and spurious associations, we transform the
# adjacency into Topological Overlap Matrix,
# and calculate the corresponding dissimilarity

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);

# calculate the corresponding dissimilarity
dissTOM = 1-TOM




# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
pdf(file = "WGCNA_Gene clustering on TOM−based dissimilarity.pdf", width = 18, height = 12);

plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()



# We like large modules, so we set the minimum module size relatively high:
# minModuleSize = 30;
minModuleSize = 1;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)




# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
# pdf(file = "WGCNA_Gene dendrogram and module colors.pdf", width = 18, height = 12);
# 
# plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Gene dendrogram and module colors")
# dev.off()



# Calculate eigengenes
MEList = moduleEigengenes(exprdata, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
# pdf(file = "WGCNA_Clustering of module eigengenes.pdf", width = 18, height = 12);
# 
# plot(METree, main = "Clustering of module eigengenes",
#      xlab = "", sub = "")
# Plot the cut line into the dendrogram

# abline(h=MEDissThres, col = "red")
# dev.off()


MEDissThres = 0

# Call an automatic merging function
merge = mergeCloseModules(exprdata, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


# pdf(file = "WGCNA_geneDendro compare merge difference.pdf", width = 18, height = 12);
# 
# 
# plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
#                     c("Dynamic Tree Cut", "Merged dynamic"),
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# dev.off()




# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "02-networkConstruction-stepByStep.RData")


############
# summary which modules interest genes (OBPs) are locate/assigned
# the gene numbers contained in each module under different runs (soft threshold 5-10)
# setwd("D:/major_project/breast_cancer")


summarize_module_information <- function(
  input_file = "/home/lyang73/data/breast_normal",
  output_symbol = "auto",
  workdir = "/home/lyang73/data/breast_cancer")
{
  setwd(workdir)
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
# output_file_names = summarize_module_information(input_file = "WGCNA_networkConstruction-auto_mydata",
#                                                  output_symbol = "auto",
#                                                  workdir = "/home/lyang73/data/breast_cancer")

# 
# output_file_names = summarize_module_information(input_file = "WGCNA_networkConstruction-auto_mydata",
#                                                  output_symbol = "auto",
#                                                  workdir = "/home/lyang73/data/breast_normal")

output_file_names = summarize_module_information(input_file = "WGCNA_networkConstruction-auto_mydata",
                                                 output_symbol = "auto",
                                                 workdir = "D:/major_project/Testis_cancer")


print(output_file_names)

# load(file = "WGCNA_01-dataInput_mydata.RData")
# library(WGCNA)
# rm(moduleLabels)
# load(file = "02-networkConstruction-stepByStep09.RData")
# names(moduleLabels) = colnames(exprdata)
# moduleLabels[names(moduleLabels) == "OBP2B"]
# moduleLabels[names(moduleLabels) == "OBP2A"]
# table(moduleLabels)



#########
# Exporting a list of  genes within a module and their connectivity scores
# scores used for gene network visualization and enrichment analysis
# Exporting a gene network within one module to the file formate compatible with network visualization software

retrive_gene_module_scores <- function(interest_module = 10, 
                                       workdir = "/home/lyang73/data/breast_cancer",
                                       interest_gene = "OBP2B",
                                       WGCNA_result_file = "WGCNA_networkConstruction-auto_mydata08.RData")
{
  setwd(workdir)
  library(WGCNA)
  
  
  load(file = WGCNA_result_file)
  load(file = "WGCNA_01-dataInput_mydata.RData")
  enableWGCNAThreads()
  # Recalculate topological overlap
  TOM = TOMsimilarityFromExpr(exprdata, power = 8);
  save(TOM,file = "TOM.RData")
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


output_file_names = retrive_gene_module_scores(interest_module = 10, 
                                               interest_gene = "OBP2B")
print(output_file_names)
output_file_names = retrive_gene_module_scores(interest_module = 10,
                                               interest_gene = "LCN2")

output_file_names = retrive_gene_module_scores(interest_module = 5,
                                               interest_gene = "OBP2A")





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
                          interest_gene = "OBP2B", interest_module = 10)
{
  setwd("D:/major_project")
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
output_file_names = plot_network("VisANTInput-10.txt","OBP2B",10)
# 0.005300881 OBP2B
# print(output_file_names)
# significant_threshold[2] = output_file_names[1]

output_file_names = plot_network("VisANTInput-5.txt","OBP2A",5)
output_file_names = plot_network("VisANTInput-10.txt","LCN2",10)
# 0.014264764778568

significant_threshold = rep(0,3)
n = 1
while (n < 4) 
{
  
}


#############
# topGO: gene ontology enrichment analysis
# use the paper script as reference, plot the top 10 enriched GO entities  
# Interfacing genes in a module with gene ontology enrichment analysis


# setwd("/home/lyang73/data")


GO_enrichment_scores <- function(input_file = "OBP2A_score.RData",
                                 interest_genes = "OBP2A",
                                 significant_threshold = 0.01)
{
  setwd("D:/major_project")
  library(WGCNA)
  library(data.table)
  library(grex);
  library(biomaRt);
  library(topGO);
  library(dplyr);
  library(ggplot2)
  # Read in the probe annotation
  # annot = read.csv(file = "00_TCGA breast Cancer_norm_all_ExpectedCnt.csv",header = TRUE)
  # annot = annot[,2:3]
  # Match probes in the data set to the probe IDs in the annotation file
  
  # load(file = "WGCNA_networkConstruction-auto_mydata08.RData")
  
  # table(moduleLabels)
  # # 
  # genes_module = names(moduleLabels)
  # module_labels = as.data.frame(moduleLabels)
  # module_labels$gene = rownames(module_labels)
  # genes_module = merge(ann33ot,module_labels,by = "gene")
  # write.csv(genes_module,file = "genes_module.csv", row.names = FALSE)
  # genes_module = read.csv("genes_module.csv")
  
  # # Annotate the Ensembl id of WGCNA result, with Entrez IDs, via function grex
  # annotComplete = grex(cleanid(genes_module$id))
  # genes_module$id = cleanid(genes_module$id)
  # annotComplete = merge(annotComplete,genes_module[,-1],by.x = "ensembl_id",by.y = "id")
  # 
  # write.csv(annotComplete, "topGO_moduleGeneAnnotated.csv", row.names = FALSE)
  # rm(annotComplete,interest_gene_score)
  annotComplete = read.csv("topGO_moduleGeneAnnotated.csv")
  
  
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
  
  # visualization of top 10 GO entities 
  # Generate a summary figure of topGO results via the ggplot2 R package
  #This sets the threshold p<0.05 and #.annotated.background.genes >= 30. 
  
  GO_bar = all_res_ks
  sapply(GO_bar,class);
  GO_bar = transform(GO_bar,KS = as.numeric(KS));
  sapply(GO_bar,class)
  
  GO_bar[is.na(GO_bar$KS),"KS"] =  1e-30
  GO_bar[is.nan(GO_bar$OR),"OR"] = 0
  GO_bar[GO_bar$OR == -Inf,"OR"] =  0
  GO_bar = GO_bar %>%
    arrange(KS) 
  GO_bar = GO_bar[1:10,]
  
  ggplot(GO_bar,
         aes(x = reorder(Term, -KS), y = Annotated, fill = KS))+
    geom_bar(stat = "identity", width = 0.9, color = "black") +
    coord_flip() +
    scale_fill_gradient(low="#feff2b",high="#fe0100")+
    labs(title = ~underline("Enriched GO Biological Processes, ks"),
         X = NULL,
         y = "number of annotated genes", fill = "p.value")+theme_bw()+
    theme(plot.title.position = "plot")+
    theme(plot.title = element_text(size = 12))+
    theme(axis.title.x = element_text(size = 10,face = "bold"),
          axis.title.y = element_text(size = 10,face = "bold"))+
    theme(axis.text.x = element_text(size = 10,face = "bold"),
          axis.text.y = element_text(size = 10,face = "bold"))+
    theme(legend.position = "right")
  ggsave(paste0("top 10 enriched GO for ",interest_genes," classic ks",".pdf"))
  output_file = c(output_file,paste0("top 10 enriched GO for ",interest_genes," classic ks",".pdf"))
  
  
  ## fisher results
  # sigTest = runTest(GOdata,
  #                   algorithm = "elim", statistic = "fisher")
  
  all_res_fisher = GenTable(GOdata,  Fis = resultFisher, orderBy = "Fis", 
                            topNodes = length(allGO))
  
  all_res_fisher$OR = log2((all_res_fisher$Significant/numSigGenes(GOdata))/
                             (all_res_fisher$Annotated/tot_background))
  
  GO_bar = all_res_fisher
  GO_bar = transform(GO_bar,Fis = as.numeric(Fis));
  
  GO_bar = GO_bar %>%
    arrange(Fis) 
  GO_bar = GO_bar[1:10,]
  
  ggplot(GO_bar,
         aes(x = reorder(Term, -Fis), y = OR, fill = Fis))+
    geom_bar(stat = "identity", width = 0.9, color = "black") +
    coord_flip() +
    scale_fill_gradient(low="#feff2b",high="#fe0100")+
    labs(title = ~underline("Enriched GO Biological Processes, fisher"),
         X = NULL,
         y = "Odds Ratio", fill = "p.value")+theme_bw()+
    theme(plot.title.position = "plot")+
    theme(plot.title = element_text(size = 12))+
    theme(axis.title.x = element_text(size = 10,face = "bold"),
          axis.title.y = element_text(size = 10,face = "bold"))+
    theme(axis.text.x = element_text(size = 10,face = "bold"),
          axis.text.y = element_text(size = 10,face = "bold"))+
    theme(legend.position = "right")
  ggsave(paste0("top 10 enriched GO for ",interest_genes," classic fisher",".pdf"))
  output_file = c(output_file,
                  paste0("top 10 enriched GO for ",interest_genes," classic fisher",".pdf"))
  return(output_file)
}

output_file_names = GO_enrichment_scores(input_file = "LCN2_score.RData",
                                         interest_genes = "LCN2",
                                         significant_threshold = 0.014)

output_file_names = GO_enrichment_scores(input_file = "OBP2A_score.RData",
                                         interest_genes = "OBP2A",
                                         significant_threshold = 0.01)

output_file_names = GO_enrichment_scores(input_file = "OBP2B_score.RData",
                                         interest_genes = "OBP2B",
                                         significant_threshold = 0.005)



GO_enrichment_module <- function(interest_module = 5)
{
  # Read in the probe annotation
  # annot = read.csv(file = "00_TCGA breast Cancer_norm_all_ExpectedCnt.csv",header = TRUE)
  # annot = annot[,2:3]
  # Match probes in the data set to the probe IDs in the annotation file
  
  # load(file = "WGCNA_networkConstruction-auto_mydata08.RData")
  
  # table(moduleLabels)
  # # 
  # genes_module = names(moduleLabels)
  # module_labels = as.data.frame(moduleLabels)
  # module_labels$gene = rownames(module_labels)
  # genes_module = merge(ann33ot,module_labels,by = "gene")
  # write.csv(genes_module,file = "genes_module.csv", row.names = FALSE)
  # genes_module = read.csv("genes_module.csv")
  
  # # Annotate the Ensembl id of WGCNA result, with Entrez IDs, via function grex
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
  
  
  # visualization of top 10 GO entities 
  # Generate a summary figure of topGO results via the ggplot2 R package
  #This sets the threshold p<0.05 and #.annotated.background.genes >= 30. 
  # GO_bar = all_res %>%
  #   filter(weightFisher < 0.05) %>% filter(Annotated >= 30);
  GO_bar = all_res 
  
  # GO_bar = GO_bar %>%
  #   select(Term, weightFisher, OR);
  
  GO_bar = transform(GO_bar,weightFisher = as.numeric(weightFisher));
  
  #This selects the top 10 most significant(i,e,,orcler by ascending p-values) GO terms for presentation. 
  GO_bar = GO_bar %>%
    arrange(weightFisher) %>% slice(1:10);
  
  
  ggplot(GO_bar,
         aes(x = reorder(Term, -weightFisher), y = OR, fill = weightFisher))+
    geom_bar(stat = "identity", width = 0.9, color = "black") +
    coord_flip() +
    scale_fill_gradient(low="#feff2b",high="#fe0100")+
    ylim(0,8)+
    labs(title = ~underline("Enriched GO Biological Processes, fisher"),
         X = NULL,
         y = "Odds Ratio", fill = "p.value")+theme_bw()+
    theme(plot.title.position = "plot")+
    theme(plot.title = element_text(size = 12))+
    theme(axis.title.x = element_text(size = 10,face = "bold"),
          axis.title.y = element_text(size = 10,face = "bold"))+
    theme(axis.text.x = element_text(size = 10,face = "bold"),
          axis.text.y = element_text(size = 10,face = "bold"))+
    theme(legend.position = "right");
  ggsave(paste0("top 10 enriched GO for module",interest_module,"classic fisher",".pdf"))
  
  
  
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
  
  #This selects the top 10 most significant(i,e,,orcler by ascending p-values) GO terms for presentation. 
  GO_bar = GO_bar %>%
    arrange(weightFisher) %>% slice(1:10);
  
  
  ggplot(GO_bar,
         aes(x = reorder(Term, -weightFisher), y = OR, fill = weightFisher))+
    geom_bar(stat = "identity", width = 0.9, color = "black") +
    coord_flip() +
    scale_fill_gradient(low="#feff2b",high="#fe0100")+
    ylim(0,8)+
    labs(title = ~underline("Enriched GO Biological Processes, fisher"),
         X = NULL,
         y = "Odds Ratio", fill = "p.value")+theme_bw()+
    theme(plot.title.position = "plot")+
    theme(plot.title = element_text(size = 12))+
    theme(axis.title.x = element_text(size = 10,face = "bold"),
          axis.title.y = element_text(size = 10,face = "bold"))+
    theme(axis.text.x = element_text(size = 10,face = "bold"),
          axis.text.y = element_text(size = 10,face = "bold"))+
    theme(legend.position = "right");
  ggsave(paste0("top 10 enriched GO for module",interest_module,"elim fisher.pdf"))
  
  output_file_names = c(paste0("GOdata_module_",interest_module,".RData"),
                        paste0("top 10 enriched GO for module",interest_module,"classic fisher",".pdf"),
                        paste0("top 10 enriched GO for module",interest_module,"elim fisher",".pdf"))
  return(output_file_names)
}
output_file_names = GO_enrichment_module(5)
output_file_names = GO_enrichment_module(10)
##########
# sanity check, pick the top 10 co-expressed genes with OBP2A, plot the pearson 
# correlation between these genes, check whether their correlation scores are high
library(dplyr)
library(limma) 
library(edgeR)
library(ggplot2)
library(corrplot)



setwd("D:/major_project")

load("OBP2A_score.RData")

# use topological overlap matrix as weights
edgeList_OBP = data.frame(Weight = interest_gene_score, genes = names(interest_gene_score))

candidates = edgeList_OBP[edgeList_OBP$genes %in% c("LCN2","OBP2A","OBP2B","OR51E1","OR2B6","OR51E2"),]
# edgeList_OBP = read.csv("OBP2A_top30_correlated_genes.csv")
edgeList_OBP = edgeList_OBP[sort(edgeList_OBP$Weight, decreasing = TRUE, index.return=TRUE)$ix[1:10],]
# 
# edgeList_OBP$TargetName[edgeList_OBP$SourceName != "OBP2A"] = edgeList_OBP[edgeList_OBP$SourceName != "OBP2A","SourceName"]
# edgeList_OBP$SourceName = "OBP2A"
interest_genes = edgeList_OBP$genes
write.csv(edgeList_OBP,file ="OBP2A_top10_correlated_genes.csv",row.names = FALSE)


norm_exprFinal = read.csv("00_TCGA breast Cancer_norm_all_ExpectedCnt.csv")

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

pdf(file = "sanity check of correlation in breast cancer.pdf")

corrplot(M,order = "AOE",addCoef.col = "white", mar=c(0,0,1,0),type = "lower",
         title = paste("expression correlation between genes in","breast cancer"))

dev.off()


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
         title = paste("correlation between co-expression scores of genes in breast cancer"))


dev.off()

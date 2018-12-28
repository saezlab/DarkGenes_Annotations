# ============================================== #
# TF and pathway enrichments of dark genes using
# DoRothEA (DRT) and PROGENy (PGN)
# Author: Panuwat Trairatphisan
# Last update: 03.09.2018
# ============================================== #

# Clean workspace
rm(list=ls()); cat("\014"); if(length(dev.list())>0){dev.off()}

# Create "results" directory to store results from the analysis
dir.create("results",showWarnings = FALSE)

# Load GO enrichment R-functions prepared by Janet Pinero
source("resources/GOenrichment_Janet.R")

# Load gene annotation file for Homo Sapien from Entrez 
# ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
library(readr)
background.file <- read_delim("resources/Homo_sapiens.gene_info", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
colnames(background.file)[1] <- "X.tax_id"
universe <- as.character(subset(background.file, X.tax_id == 9606)$GeneID)

# Analyse the proportion of genes that has GO/GO-BP annotations
print("universe just imported has length of:")
length(universe)
# [1] 60196
universeGO <- haveGO(universe)
print("universe haveGO has length of:")
length(universeGO)
# [1] 19303
universeGOBP <-  filterbyGO(universeGO, "BP")
length(universeGOBP)
# [1] 16986 (status 24.04.18) // 17380 (status 03.09.18)
# universeBP_Symbol <- NULL
# for (counter in 1:length(universeGOBP)) {
#   print(paste0("Mapping GeneSymbol ",counter,"/",length(universeGOBP)))
#   universeBP_Symbol[counter] <- background.file$Symbol[universeGOBP[counter]==background.file$GeneID]
# }

# Load gene set from TG-GATEs dataset with all annotated information
NetworkGeneList <- read.delim("resources/All_TGGATEs_genes.tsv", quote="", stringsAsFactors=FALSE)
ReactomeGeneList <- NetworkGeneList$Symbol[NetworkGeneList$reactome==1]
OmnipathGeneList <- NetworkGeneList$Symbol[NetworkGeneList$omnipath==1]
MSigDBGeneList <- NetworkGeneList$Symbol[NetworkGeneList$msigdb==1]
PWCommonGeneList <- NetworkGeneList$Symbol[NetworkGeneList$pathcom==1]

# =============================================================================== #
# === Check DoRothEA regulon with universeBP_Symbol AND interaction databases === #
# =============================================================================== #

# Load TF-regulons list from DoRothEA version 1
# https://github.com/saezlab/DoRothEA/releases/tag/version1
load('resources/CTFRs_v122016.rdata') 

# Extraction of DoRothEA version 1's regulon
Regulons <- NULL
for (counter in 1:length(CTFRs_genesets$GENES)) {
  if (length(CTFRs_genesets$GENES[[counter]])>0) {
    Regulons <- c(Regulons,CTFRs_genesets$GENES[[counter]])
  }
}
Regulons_Unique <- unique(names(Regulons))
Regulons_Unique <- matrix(Regulons_Unique,length(Regulons_Unique),1)
length(Regulons_Unique) 
# [1] 7996 -> Total number of DoRothEA regulons
colnames(Regulons_Unique) <- "hgnc_symbol"
# write.table(x = Regulons_Unique,file = "Regulons_from_DoRothEA.tsv",quote = F,row.names = F,col.names = F)

# Map DoRothEA regulon to the list of all genes with GO-BP annotation & not in the selected interaction databases (via GeneSymbol)
# [Note: The analysis on GO-BP was done in April 2018 so we use the results from that moment as the reference]
universeBP_Symbol <- read.table(file = "resources/universeBP_Symbol_April2018.tsv",header = F,col.names = F,stringsAsFactors = F,sep = ",")
universeBP_Symbol <- as.vector(t(universeBP_Symbol))
DoRothEA_Regulons_with_GOBP <- intersect(universeBP_Symbol,Regulons_Unique)
length(DoRothEA_Regulons_with_GOBP) 
# [1] 7416 -> Number of DoRothEA regulons with GO-BP term
IDX_DRTreg_wGOBP <- NULL
for (counter in 1:length(DoRothEA_Regulons_with_GOBP)) {
  # print(paste0("Now mapping: ",counter,"/ ",length(DoRothEA_Regulons_with_GOBP)))
  if (sum(DoRothEA_Regulons_with_GOBP[counter]==Regulons_Unique)>0) {
    IDX_DRTreg_wGOBP <- c(IDX_DRTreg_wGOBP,which(DoRothEA_Regulons_with_GOBP[counter]==Regulons_Unique))
  }
}
DoRothEA_Regulons_without_GOBP <- Regulons_Unique[-IDX_DRTreg_wGOBP]
length(DoRothEA_Regulons_without_GOBP)
# [1] 580 -> Number of DoRothEA regulons without GO-BP term

# 580/7996*100
# [1] 7.253627 -> 7% of DRT regulons do not have GO-BP term

DoRothEA_Regulons_without_GOBP_NotInNetwork1 <- setdiff(DoRothEA_Regulons_without_GOBP,ReactomeGeneList)
DoRothEA_Regulons_without_GOBP_NotInNetwork2 <- setdiff(DoRothEA_Regulons_without_GOBP_NotInNetwork1,OmnipathGeneList)
DoRothEA_Regulons_without_GOBP_NotInNetwork3 <- setdiff(DoRothEA_Regulons_without_GOBP_NotInNetwork2,MSigDBGeneList)
DoRothEA_Regulons_without_GOBP_NotInNetwork_All <- setdiff(DoRothEA_Regulons_without_GOBP_NotInNetwork3,PWCommonGeneList)
length(DoRothEA_Regulons_without_GOBP_NotInNetwork_All)
# [1] 485 -> Number of DoRothEA regulons without GO-BP term & not in interaction databases

# 485/7996*100
# [1] 6.065533 -> 6% of DRT regulons do not have GO-BP term AND are not in interaction databases

# ================================================================================= #
# === Check PROGENy signatures with universeBP_Symbol AND interaction databases === #
# ================================================================================= #

# load PROGENy signatures
# https://github.com/saezlab/progeny/tree/master/data
load('resources/model.rda')
ProgenyGenes <- matrix(rownames(model),nrow(model),1)
length(ProgenyGenes)
# [1] 1059 -> Total number of PROGENy signatures 
colnames(ProgenyGenes) <- "hgnc_symbol"
# write.table(x = ProgenyGenes,file = "ProgenyGenes.tsv",quote = F,row.names = F,col.names = F)

# Map PROGENy signatures to the list of all genes with GO-BP annotation & not in the selected interaction databases (via GeneSymbol)
# [Note: The analysis on GO-BP was done in April 2018 so we use the results from that moment as the reference]
universeBP_Symbol <- read.table(file = "resources/universeBP_Symbol_April2018.tsv",header = F,col.names = F,stringsAsFactors = F,sep = ",")
universeBP_Symbol <- as.vector(t(universeBP_Symbol))
PROGENy_Signatures_with_GOBP <- intersect(universeBP_Symbol,ProgenyGenes)
length(PROGENy_Signatures_with_GOBP)
# [1] 862 -> Number of PROGENy signatures with GO-BP term
IDX_PGNsig_wGOBP <- NULL
for (counter in 1:length(PROGENy_Signatures_with_GOBP)) {
  # print(paste0("Now mapping: ",counter,"/ ",length(PROGENy_Signatures_with_GOBP)))
  if (sum(PROGENy_Signatures_with_GOBP[counter]==ProgenyGenes)>0) {
    IDX_PGNsig_wGOBP <- c(IDX_PGNsig_wGOBP,which(PROGENy_Signatures_with_GOBP[counter]==ProgenyGenes))
  }
}
PROGENy_Signatures_without_GOBP <- ProgenyGenes[-IDX_PGNsig_wGOBP]
length(PROGENy_Signatures_without_GOBP)
# [1] 197 -> Number of PROGENy signature without GO-BP term

# 197/1059*100
# [1] 18.60246 -> 18% of PGN signatures do not have GO-BP term

PROGENy_Signatures_without_GOBP_NotInNetwork1 <- setdiff(PROGENy_Signatures_without_GOBP,ReactomeGeneList)
PROGENy_Signatures_without_GOBP_NotInNetwork2 <- setdiff(PROGENy_Signatures_without_GOBP_NotInNetwork1,OmnipathGeneList)
PROGENy_Signatures_without_GOBP_NotInNetwork3 <- setdiff(PROGENy_Signatures_without_GOBP_NotInNetwork2,MSigDBGeneList)
PROGENy_Signatures_without_GOBP_NotInNetwork_All <- setdiff(PROGENy_Signatures_without_GOBP_NotInNetwork3,PWCommonGeneList)
length(PROGENy_Signatures_without_GOBP_NotInNetwork_All)
# [1] 185 -> Number of PROGENy signature without GO-BP term AND are not in interaction databases

# 185/1059*100
# [1] 17.46931 -> 17% of PGN signatures do not have GO-BP term AND are not in interaction databases

# ================================================== #
# === Systematic annotation with GO-BP & DRT/PGN === #
# ================================================== #

# Load full information of gene set
NetworkGeneList <- read.delim("resources/All_TGGATEs_genes.tsv", quote="", stringsAsFactors=FALSE)
DarkGeneTableAll <- NetworkGeneList[which(NetworkGeneList$GOBP==0 & NetworkGeneList$pathcom==0 & NetworkGeneList$omnipath==0 & NetworkGeneList$msigdb==0 & NetworkGeneList$reactome==0),]
DarkGeneTable <- NetworkGeneList[which(NetworkGeneList$GOBP==0 & NetworkGeneList$pathcom==0 & NetworkGeneList$omnipath==0 & NetworkGeneList$msigdb==0 & NetworkGeneList$reactome==0 & NetworkGeneList$DEG==1),]
nrow(DarkGeneTable)
# [1] 916 -> There are 916 dark genes based on our definition

# Set logFC cutoff
FCcutoff <- 0.585 # LogFC of 1.5

# Toxic compounds:  
ToxCompound <- c("acetaminophen", "valproic_acid", "isoniazid", "diclofenac","nimesulide")

ToxGenesAll <- NULL
ToxGenesAnnoSummary <- matrix(NA,length(ToxCompound)+1,4)
colnames(ToxGenesAnnoSummary) <- c("Compound","NrGenes","UnmappedDRT(580)","UnmappedPGN(197)")
ToxGenesAnnoSummary[,1] <- c(ToxCompound,"AllToxCompounds")

AllToxSymbol <- NULL; AllToxDRT <- NULL; AllToxPGN <- NULL

for (counter in 1:length(ToxCompound)) {

  # Data from Janet's new DarkGene list 
  CompoundColID <- which(ToxCompound[counter]==colnames(DarkGeneTable))
  ToxGeneID <- DarkGeneTable[which(DarkGeneTable[,CompoundColID]==1),1]
  
  # Map to HGNC symbol
  ToxGeneSymbol <- NULL
  for (counter2 in 1:length(ToxGeneID)) {
    print(paste0("Now mapping: ",counter2," / ",length(ToxGeneID)))
    ToxGeneSymbol <- c(ToxGeneSymbol,background.file$Symbol[which(background.file$GeneID==ToxGeneID[counter2])])
  }
  AllToxSymbol <- c(AllToxSymbol,ToxGeneSymbol)
  
  # Map to DoRothEA and PROGENy
  Overlapped_Unmapped_ToxGene_DRT <- intersect(ToxGeneSymbol,DoRothEA_Regulons_without_GOBP)
  Overlapped_Unmapped_ToxGene_PGN <- intersect(ToxGeneSymbol,PROGENy_Signatures_without_GOBP)
  AllToxDRT <- c(AllToxDRT,Overlapped_Unmapped_ToxGene_DRT)
  AllToxPGN <- c(AllToxPGN,Overlapped_Unmapped_ToxGene_PGN)
  
  # Store results in the pre-defined table
  ToxGenesAnnoSummary[counter,2] <- length(ToxGeneSymbol)
  ToxGenesAnnoSummary[counter,3] <- length(Overlapped_Unmapped_ToxGene_DRT)
  ToxGenesAnnoSummary[counter,4] <- length(Overlapped_Unmapped_ToxGene_PGN)
  
  # Write results as files
  write.table(x = ToxGeneSymbol,file = paste0("results/",ToxCompound[counter],"_DEG_Symbol.tsv"),quote = F,sep = "\n",row.names = F,col.names = F)
  write.table(x = Overlapped_Unmapped_ToxGene_DRT,file = paste0("results/",ToxCompound[counter],"_NoGOBP_DRT.tsv"),quote = F,sep = "\n",row.names = F,col.names = F)
  write.table(x = Overlapped_Unmapped_ToxGene_PGN,file = paste0("results/",ToxCompound[counter],"_NoGOBP_PGN.tsv"),quote = F,sep = "\n",row.names = F,col.names = F)
  
}

ToxGenesAnnoSummary[length(ToxCompound)+1,2] <- length(unique(AllToxSymbol))
ToxGenesAnnoSummary[length(ToxCompound)+1,3] <- length(unique(AllToxDRT))
ToxGenesAnnoSummary[length(ToxCompound)+1,4] <- length(unique(AllToxPGN))

print(ToxGenesAnnoSummary)
write.table(x = ToxGenesAnnoSummary,file = paste0("results/Summary_ToxGene_Annotation.tsv"),quote = F,sep = "\t",row.names = F,col.names = T)

# Non-toxic compounds: 
NonToxCompound <- c("theophylline", "caffeine", "hydroxyzine", "chloramphenicol", "chlorpheniramine")

NonToxGenesAll <- NULL
NonToxGenesAnnoSummary <- matrix(NA,length(NonToxCompound)+1,4)
colnames(NonToxGenesAnnoSummary) <- c("Compound","NrGenes","UnmappedDRT(580)","UnmappedPGN(197)")
NonToxGenesAnnoSummary[,1] <- c(NonToxCompound,"AllNonToxCompounds")
AllNonToxSymbol <- NULL; AllNonToxDRT <- NULL; AllNonToxPGN <- NULL

for (counter in 1:length(NonToxCompound)) {
  
  # Data from Janet's new DarkGene list 
  NonCompoundColID <- which(NonToxCompound[counter]==colnames(DarkGeneTable))
  NonToxGeneID <- DarkGeneTable[which(DarkGeneTable[,NonCompoundColID]==1),1]
  
  # Map to HGNC symbol
  NonToxGeneSymbol <- NULL
  for (counter2 in 1:length(NonToxGeneID)) {
    print(paste0("Now mapping: ",counter2," / ",length(NonToxGeneID)))
    NonToxGeneSymbol <- c(NonToxGeneSymbol,background.file$Symbol[which(background.file$GeneID==NonToxGeneID[counter2])])
  }
  AllNonToxSymbol <- c(AllNonToxSymbol,NonToxGeneSymbol)
  
  # Map to DoRothEA and PROGENy
  Overlapped_Unmapped_NonToxGene_DRT <- intersect(NonToxGeneSymbol,DoRothEA_Regulons_without_GOBP)
  Overlapped_Unmapped_NonToxGene_PGN <- intersect(NonToxGeneSymbol,PROGENy_Signatures_without_GOBP)
  AllNonToxDRT <- c(AllNonToxDRT,Overlapped_Unmapped_NonToxGene_DRT)
  AllNonToxPGN <- c(AllNonToxPGN,Overlapped_Unmapped_NonToxGene_PGN)
  
  # Store results in the pre-defined table
  NonToxGenesAnnoSummary[counter,2] <- length(NonToxGeneSymbol)
  NonToxGenesAnnoSummary[counter,3] <- length(Overlapped_Unmapped_NonToxGene_DRT)
  NonToxGenesAnnoSummary[counter,4] <- length(Overlapped_Unmapped_NonToxGene_PGN)
  
  # Write results as files
  write.table(x = NonToxGeneSymbol,file = paste0("results/",NonToxCompound[counter],"_DEG_Symbol.tsv"),quote = F,sep = "\n",row.names = F,col.names = F)
  write.table(x = Overlapped_Unmapped_NonToxGene_DRT,file = paste0("results/",NonToxCompound[counter],"_NoGOBP_DRT.tsv"),quote = F,sep = "\n",row.names = F,col.names = F)
  write.table(x = Overlapped_Unmapped_NonToxGene_PGN,file = paste0("results/",NonToxCompound[counter],"_NoGOBP_PGN.tsv"),quote = F,sep = "\n",row.names = F,col.names = F)
  
}

NonToxGenesAnnoSummary[length(NonToxCompound)+1,2] <- length(unique(AllNonToxSymbol))
NonToxGenesAnnoSummary[length(NonToxCompound)+1,3] <- length(unique(AllNonToxDRT))
NonToxGenesAnnoSummary[length(NonToxCompound)+1,4] <- length(unique(AllNonToxPGN))

print(NonToxGenesAnnoSummary)
write.table(x = NonToxGenesAnnoSummary,file = paste0("results/Summary_NonToxGene_Annotation.tsv"),quote = F,sep = "\t",row.names = F,col.names = T)

# Merged results Tox/NonTox genes
write.table(x = rbind(ToxGenesAnnoSummary,NonToxGenesAnnoSummary),file = paste0("results/Summary_AllGene_Annotation.tsv"),quote = F,sep = "\t",row.names = F,col.names = T)


# ======================================================= #
# --- Now map back which TFs and pathways are involved -- #
# ======================================================= #

# Load results from Hepatotoxic (HPTx) compounds and annotate back

ToxCompound <- c("acetaminophen", "valproic_acid", "isoniazid", "diclofenac","nimesulide")
AllToxDRT <- NULL; AllToxPGN <- NULL

ToxGenesAnnoSummary <- matrix(NA,length(ToxCompound)+1,3)
colnames(ToxGenesAnnoSummary) <- c("Compound","Anno'ed-TF","Anno'ed-PW")
ToxGenesAnnoSummary[,1] <- c(ToxCompound,"AllToxCompounds")

for (counter_cpd in 1:length(ToxCompound)) {
  # DoRothEA
  DRT <- t(read.table(paste0("results/",ToxCompound[counter_cpd],"_NoGOBP_DRT.tsv"),sep = "\n"))
  DRT_list <- list()
  for (counter in 1:length(DRT)) {
    print(paste0("Processing ",counter,"/",length(DRT)))
    DRT_list[[counter]] <- list(name = DRT[counter], TF= NULL)
    
    for (counter2 in 1:length(CTFRs_genesets$GENES)) {
      if (sum(DRT[counter]==names(CTFRs_genesets$GENES[counter2][[1]]))>0) {
        DRT_list[[counter]]$TF <- c(DRT_list[[counter]]$TF,names(CTFRs_genesets$GENES[counter2]))
      }
    }
  }
  
  DRT_TF <- NULL
  for (counter3 in 1:length(DRT_list)) {
    DRT_TF <- c(DRT_TF,DRT_list[[counter3]]$TF)
  }
  DRT_TF_Table <- matrix(NA,2,length(unique(DRT_TF)))
  DRT_TF_Table[1,] <- names(table(DRT_TF))
  DRT_TF_Table[2,] <- table(DRT_TF)
  write.table(x = DRT_TF_Table,file = paste0("results/",ToxCompound[counter_cpd],"_TFlist_DRT.tsv"),quote = F,sep = "\t",row.names = F,col.names = F)
  AllToxDRT <- c(AllToxDRT,names(table(DRT_TF)))
  
  # PROGENy
  PGN <- t(read.table(paste0("results/",ToxCompound[counter_cpd],"_NoGOBP_PGN.tsv"),sep = "\n"))
  
  PGN_list <- list()
  for (counter in 1:length(PGN)) {
    print(paste0("Processing ",counter,"/",length(PGN)))
    PGN_list[[counter]] <- list(name = PGN[counter], Signature= NULL)
    
    # for (counter2 in 1:nrow(model)) {
    #   if (sum(PGN[counter]==rownames(model))>0) {
    PGN_list[[counter]]$Signature <- c(PGN_list[[counter]]$Signature,colnames(model)[which(model[which(PGN[counter]==rownames(model)),]!=0)])
    #   }
    # }
  }
  
  PGN_Sig <- NULL
  for (counter3 in 1:length(PGN_list)) {
    PGN_Sig <- c(PGN_Sig,PGN_list[[counter3]]$Signature)
  }
  PGN_Sig_Table <- matrix(NA,2,length(unique(PGN_Sig)))
  PGN_Sig_Table[1,] <- names(table(PGN_Sig))
  PGN_Sig_Table[2,] <- table(PGN_Sig)
  write.table(x = PGN_Sig_Table,file = paste0("results/",ToxCompound[counter_cpd],"_SigList_PGN.tsv"),quote = F,sep = "\t",row.names = F,col.names = F)
  AllToxPGN <- c(AllToxPGN,names(table(PGN_Sig)))
  
  ToxGenesAnnoSummary[counter_cpd,2] <- ncol(DRT_TF_Table)
  ToxGenesAnnoSummary[counter_cpd,3] <- ncol(PGN_Sig_Table)
  
}

ToxGenesAnnoSummary[length(ToxCompound)+1,2] <- length(sort(unique(AllToxDRT)))
ToxGenesAnnoSummary[length(ToxCompound)+1,3] <- length(sort(unique(AllToxPGN)))

DRT_TF_Table <- matrix(NA,2,length(unique(AllToxDRT)))
DRT_TF_Table[1,] <- names(table(AllToxDRT))
DRT_TF_Table[2,] <- table(AllToxDRT)
write.table(x = DRT_TF_Table,file = paste0("results/Pooled_ToxGene_TFlist_DRT.tsv"),quote = F,sep = "\t",row.names = F,col.names = F)

PGN_Sig_Table <- matrix(NA,2,length(unique(AllToxPGN)))
PGN_Sig_Table[1,] <- names(table(AllToxPGN))
PGN_Sig_Table[2,] <- table(AllToxPGN)
write.table(x = PGN_Sig_Table,file = paste0("results/Pooled_ToxGene_SigList_PGN.tsv"),quote = F,sep = "\t",row.names = F,col.names = F)

print(ToxGenesAnnoSummary)
write.table(x = ToxGenesAnnoSummary,file = paste0("results/Summary_ToxGene_TF_Sig.tsv"),quote = F,sep = "\t",row.names = F,col.names = T)

# Load results from non-HPTx compounds and annotate back
NonToxCompound <- c("theophylline", "caffeine", "hydroxyzine", "chloramphenicol", "chlorpheniramine")
AllNonToxDRT <- NULL; AllNonToxPGN <- NULL

NonToxGenesAnnoSummary <- matrix(NA,length(NonToxCompound)+1,3)
colnames(NonToxGenesAnnoSummary) <- c("Compound","Anno'ed-TF","Anno'ed-PW")
NonToxGenesAnnoSummary[,1] <- c(NonToxCompound,"AllNonToxCompounds")

for (counter_cpd in 1:length(NonToxCompound)) {
  
  # DoRothEA
  
  if (file.info(paste0("results/",NonToxCompound[counter_cpd],"_NoGOBP_DRT.tsv"))[1]!=0) {
    
    DRT <- t(read.table(paste0("results/",NonToxCompound[counter_cpd],"_NoGOBP_DRT.tsv"),sep = "\n"))
    DRT_list <- list()
    for (counter in 1:length(DRT)) {
      print(paste0("Processing ",counter,"/",length(DRT)))
      DRT_list[[counter]] <- list(name = DRT[counter], TF= NULL)
      
      for (counter2 in 1:length(CTFRs_genesets$GENES)) {
        if (sum(DRT[counter]==names(CTFRs_genesets$GENES[counter2][[1]]))>0) {
          DRT_list[[counter]]$TF <- c(DRT_list[[counter]]$TF,names(CTFRs_genesets$GENES[counter2]))
        }
      }
    }
    
    DRT_TF <- NULL
    for (counter3 in 1:length(DRT_list)) {
      DRT_TF <- c(DRT_TF,DRT_list[[counter3]]$TF)
    }
    DRT_TF_Table <- matrix(NA,2,length(unique(DRT_TF)))
    DRT_TF_Table[1,] <- names(table(DRT_TF))
    DRT_TF_Table[2,] <- table(DRT_TF)
    write.table(x = DRT_TF_Table,file = paste0("results/",NonToxCompound[counter_cpd],"_TFlist_DRT.tsv"),quote = F,sep = "\t",row.names = F,col.names = F)
    AllNonToxDRT <- c(AllNonToxDRT,names(table(DRT_TF)))
    
    NonToxGenesAnnoSummary[counter_cpd,2] <- ncol(DRT_TF_Table)
  } else {
    NonToxGenesAnnoSummary[counter_cpd,2] <- 0
    
  }
  
  # PROGENy
  
  if (file.info(paste0("results/",NonToxCompound[counter_cpd],"_NoGOBP_PGN.tsv"))[1]!=0) {
    
    PGN <- t(read.table(paste0("results/",NonToxCompound[counter_cpd],"_NoGOBP_PGN.tsv"),sep = "\n"))
    
    PGN_list <- list()
    for (counter in 1:length(PGN)) {
      print(paste0("Processing ",counter,"/",length(PGN)))
      PGN_list[[counter]] <- list(name = PGN[counter], Signature= NULL)
      
      # for (counter2 in 1:nrow(model)) {
      #   if (sum(PGN[counter]==rownames(model))>0) {
      PGN_list[[counter]]$Signature <- c(PGN_list[[counter]]$Signature,colnames(model)[which(model[which(PGN[counter]==rownames(model)),]!=0)])
      #   }
      # }
    }
    
    PGN_Sig <- NULL
    for (counter3 in 1:length(PGN_list)) {
      PGN_Sig <- c(PGN_Sig,PGN_list[[counter3]]$Signature)
    }
    PGN_Sig_Table <- matrix(NA,2,length(unique(PGN_Sig)))
    PGN_Sig_Table[1,] <- names(table(PGN_Sig))
    PGN_Sig_Table[2,] <- table(PGN_Sig)
    write.table(x = PGN_Sig_Table,file = paste0("results/",NonToxCompound[counter_cpd],"_SigList_PGN.tsv"),quote = F,sep = "\t",row.names = F,col.names = F)
    AllNonToxPGN <- c(AllNonToxPGN,names(table(PGN_Sig)))
    
    NonToxGenesAnnoSummary[counter_cpd,3] <- ncol(PGN_Sig_Table)
    
  } else {
    NonToxGenesAnnoSummary[counter_cpd,3] <- 0
  }
  
  
}

NonToxGenesAnnoSummary[length(NonToxCompound)+1,2] <- length(sort(unique(AllNonToxDRT)))
NonToxGenesAnnoSummary[length(NonToxCompound)+1,3] <- length(sort(unique(AllNonToxPGN)))

DRT_TF_Table <- matrix(NA,2,length(unique(AllNonToxDRT)))
DRT_TF_Table[1,] <- names(table(AllNonToxDRT))
DRT_TF_Table[2,] <- table(AllNonToxDRT)
write.table(x = DRT_TF_Table,file = paste0("results/Pooled_NonToxGene_TFlist_DRT.tsv"),quote = F,sep = "\t",row.names = F,col.names = F)

PGN_Sig_Table <- matrix(NA,2,length(unique(AllNonToxPGN)))
PGN_Sig_Table[1,] <- names(table(AllNonToxPGN))
PGN_Sig_Table[2,] <- table(AllNonToxPGN)
write.table(x = PGN_Sig_Table,file = paste0("results/Pooled_NonToxGene_SigList_PGN.tsv"),quote = F,sep = "\t",row.names = F,col.names = F)

write.table(x = NonToxGenesAnnoSummary,file = paste0("results/Summary_NonToxGene_TF_Sig.tsv"),quote = F,sep = "\t",row.names = F,col.names = T)

# Merged results Tox/NonTox genes
print(rbind(ToxGenesAnnoSummary,NonToxGenesAnnoSummary))
write.table(x = rbind(ToxGenesAnnoSummary,NonToxGenesAnnoSummary),file = paste0("results/Summary_AllGene_TF_Sig.tsv"),quote = F,sep = "\t",row.names = F,col.names = T)

# Additional annotation of TF and pathways on the dark gene table
DarkGeneTablePlus <- cbind(DarkGeneTable,rep(NA,nrow(DarkGeneTable)),rep(NA,nrow(DarkGeneTable)))
colnames(DarkGeneTablePlus)[(ncol(DarkGeneTablePlus)-1):ncol(DarkGeneTablePlus)] <- c("Regulating_TFs","Enriched_Pathways")

for (counter in 1:nrow(DarkGeneTable)) {
  
  print(paste0("Mapping DarkGene Nr: ",counter,"/",nrow(DarkGeneTablePlus)," - ",DarkGeneTable$Symbol[counter]))
  darkgene_symbol <- DarkGeneTable$Symbol[counter]
  
  # Map TF from DoRothEA
  
  TempTFList <- NULL
  for (counter2 in 1:length(CTFRs_genesets$GENES)) {
    if (sum(grepl(DarkGeneTable$Symbol[counter],names(CTFRs_genesets$GENES[counter2][[1]]),fixed = T))>0) {
      TempTFList <- c(TempTFList,names(CTFRs_genesets$GENES[counter2]))
    }
  }
  # print(TempTFList)
  if (!is.null(TempTFList)) {
    TempTFList <- paste(TempTFList, collapse="|")
    DarkGeneTablePlus$Regulating_TFs[counter] <- TempTFList
  } else {
    DarkGeneTablePlus$Regulating_TFs[counter] <- "--"
  }
  
  # Map Pathway from PROGENy
  
  TempPWList <- NULL
  for (counter3 in 1:ncol(model)) {
    if (sum(grepl(DarkGeneTable$Symbol[counter],rownames(model)[which(model[,counter3]!=0)],fixed = T))>0) {
      TempPWList <- c(TempPWList,colnames(model)[counter3])
    }
  }
  # print(TempPWList)
  if (!is.null(TempPWList)) {
    TempPWList <- paste(TempPWList, collapse="|")
    DarkGeneTablePlus$Enriched_Pathways[counter] <- TempPWList
  } else {
    DarkGeneTablePlus$Enriched_Pathways[counter] <- "--"
  }
  
}

print(DarkGeneTablePlus)
write.table(x = DarkGeneTablePlus,file = paste0("results/All_Annotated_DarkGenes.tsv"),quote = F,sep = "\t",row.names = F,col.names = T)

# --- End of the script --- #


#=== libraries ====================
require(dplyr)
require(igraph)
require(GOSemSim)
require(org.Hs.eg.db)
require(pbapply)

wangMethod <- function(t1, t2, ont) {
  matrix( mapply( wangMethod_internal,
                  rep( t1, length(t2) ),
                  rep( t2, each=length(t1) ),
                  MoreArgs = list( ont = ont ) ),
          dimnames = list( t1, t2 ), ncol=length(t2) )
}

wangMethod_internal <- function(ID1, ID2, ont) {
  if (ID1 == ID2)
    return (sim=1)
  if (ont == "DO") {
    .DOSEEnv <- get(".DOSEEnv", envir=.GlobalEnv)
    rel_df <- get("dotbl", envir=.DOSEEnv)
  } else if (ont %in% c("BP", "CC", "MF")) {
    if (!exists(".GOSemSimEnv")) .initial()
    .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)
    rel_df <- get("gotbl", envir=.GOSemSimEnv)
  } else {
    .meshesEnv <- get(".meshesEnv", envir=.GlobalEnv)
    rel_df <- get("meshtbl", envir=.meshesEnv)
  }


  sv.a <- getSV(ID1, ont, rel_df)
  sv.b <- getSV(ID2, ont, rel_df)

  if(all(is.na(sv.a)) || all(is.na(sv.b)))
    return (NA)

  idx         <- intersect(names(sv.a), names(sv.b))
  inter.sva   <- sv.a[idx]
  inter.svb   <- sv.b[idx]
  if (is.null(inter.sva) ||
      is.null(inter.svb) ||
      length(inter.sva) == 0 ||
      length(inter.svb) ==0) {
    return (NA)
  }

  sim <- sum(inter.sva,inter.svb) / sum(sv.a, sv.b)
  return(sim)
}
getSV <- function(ID, ont, rel_df, weight=NULL) {
  if (!exists(".SemSimCache")) .initial()
  .SemSimCache <- get(".SemSimCache", envir=.GlobalEnv)

  if( exists(ID, envir=.SemSimCache) ) {
    sv <- get(ID, envir=.SemSimCache)
    return(sv)
  }

  if (ont == "DO") {
    topNode <- "DOID:4"
  } else {
    topNode <- "all"
  }

  if (ID == topNode) {
    sv <- 1
    names(sv) <- topNode
    return (sv)
  }

  if (is.null(weight)) {
    weight <- c(0.8, 0.6, 0.7)
    names(weight) <- c("is_a", "part_of", "other")
  }

  rel_df <- rel_df[rel_df$Ontology == ont,]
  if (! 'relationship' %in% colnames(rel_df))
    rel_df$relationship <- "other"

  rel_df$relationship[!rel_df$relationship %in% c("is_a", "part_of")] <- "other"


  sv <- 1
  names(sv) <- ID
  allid <- ID

  idx <- which(rel_df[,1] %in% ID)
  while (length(idx) != 0) {
    p <- rel_df[idx,]
    pid <- p$parent
    allid <- c(allid, pid)

    sv <- c(sv, weight[p$relationship]*sv[p[,1]])
    names(sv) <- allid
    idx <- which(rel_df[,1] %in% pid)
  }

  sv <- sv[!is.na(names(sv))]
  sv <- sv[!duplicated(names(sv))]

  if(ont != "DO")
    sv[topNode] <- 0

  if( ! exists(ID, envir=.SemSimCache) ) {
    assign(ID,
           sv,
           envir=.SemSimCache)
  }

  return(sv)
}
.initial <- function() {
  pos <- 1
  envir <- as.environment(pos)
  assign(".GOSemSimEnv", new.env(), envir = envir)
  assign(".SemSimCache", new.env(), envir = envir)
  .GOSemSimEnv <- get(".GOSemSimEnv", envir=.GlobalEnv)

  tryCatch(utils::data(list="gotbl",
                       package="GOSemSim"))
  gotbl <- get("gotbl")
  assign("gotbl", gotbl, envir = .GOSemSimEnv)
  rm(gotbl, envir = .GlobalEnv)
}

gene2GO <- function(gene, godata, dropCodes) {
  goAnno <- godata@geneAnno
  if (! "EVIDENCE" %in% colnames(goAnno)) {
    warning("Evidence codes not found, 'drop' parameter will be ignored...")
  } else {
    goAnno <- goAnno[!goAnno$EVIDENCE %in% dropCodes,]
  }
  go <- as.character(unique(goAnno[goAnno[,1] == gene, "GO"]))
  go[!is.na(go)]
}


Degrees <- data.frame(Genes = names(V(net)),
                      degree = degree(net, v = V(net)),
                      stringsAsFactors = F) %>%
  mutate(Bin = case_when( .$degree <= 100 ~ 'G1',
                          .$degree %in% 101:200 ~ 'G2',
                          .$degree %in% 201:300 ~ 'G3',
                          .$degree %in% 301:400 ~ 'G4',
                          .$degree %in% 401:500 ~ 'G5',
                          .$degree %in% 501:600 ~ 'G6',
                          .$degree %in% 601:700 ~ 'G7',
                          .$degree %in% 701:800 ~ 'G8',
                          .$degree %in% 801:900 ~ 'G9',
                          .$degree %in% 901:1000 ~ 'G10',
                          .$degree %in% 1001:1100 ~ 'G11',
                          .$degree %in% 1101:1200 ~ 'G12',
                          .$degree %in% 1201:1300 ~ 'G13',
                          .$degree %in% 1301:1400 ~ 'G14',
                          .$degree %in% 1401:1500 ~ 'G15',
                          .$degree %in% 1501:1600 ~ 'G16',
                          .$degree %in% 1601:1700 ~ 'G17',
                          .$degree %in% 1701:1800 ~ 'G18',
                          .$degree > 1801~ 'G19',
                          TRUE ~ 'Other'))


# data prepare ------------------------------------------------------------

# Get the network nodes (as list)
NetworkNodes <- as.list(V(net)$name)

### Get the background Go BP Data for ENTREZID
BP_DATA <- godata(OrgDb = org.Hs.eg.db, ont = "BP", keytype = "ENTREZID")

MF_DATA <- godata(OrgDb = org.Hs.eg.db, ont = "MF", keytype = "ENTREZID")

CC_DATA <- godata(OrgDb = org.Hs.eg.db, ont = "CC", keytype = "ENTREZID")

saveRDS(BP_DATA,MF_DATA,CC_DATA, file=paste0("./GO_sim/backround_godata",".rds"))

### Get the pairwise similarity scores between the GO terms connected to the network genes
# Run the function for all genes in the network
Node2GO_BP <- pblapply(NetworkNodes, function(x) gene2GO(x, BP_DATA, "IEA"))
names(Node2GO_BP) <- NetworkNodes

Node2GO_MF <- pblapply(NetworkNodes, function(x) gene2GO(x, MF_DATA, "IEA"))
names(Node2GO_MF) <- NetworkNodes

Node2GO_CC <- pblapply(NetworkNodes, function(x) gene2GO(x, CC_DATA, "IEA"))
names(Node2GO_CC) <- NetworkNodes

save(Node2GO_BP,Node2GO_MF,Node2GO_CC, file=paste0("./GO_sim/Node2GO",".Rdata"))


GO_all_BP <- unique(unlist(Node2GO_BP))

GO_all_MF <- unique(unlist(Node2GO_MF))

GO_all_CC <- unique(unlist(Node2GO_CC))


PDL1_driver <- import('PDL1_driver_gene.txt')

DiseaseModules <- list(PDL1_driver$V2)

# Get the Go terms associated to the modules (list)
GO_Modules_BP <- lapply(DiseaseModules, function(x)
  unique(unlist(Node2GO_BP[names(Node2GO_BP) %in% x])))


GO_Modules_MF <- lapply(DiseaseModules, function(x)
  unique(unlist(Node2GO_MF[names(Node2GO_MF) %in% x])))


GO_Modules_CC <- lapply(DiseaseModules, function(x)
  unique(unlist(Node2GO_CC[names(Node2GO_CC) %in% x])))


save(GO_Modules_BP,GO_Modules_MF,GO_Modules_CC, file=paste0("./GO_sim/PDL1_GO_Modules",".Rdata"))

# calculate the pairwise similarity between the Module GO terms and the network GO terms

GO_pairwise_BP <- lapply(GO_Modules_BP, function(x) wangMethod(x, GO_all_BP, "BP"))

GO_pairwise_MF <- lapply(GO_Modules_MF, function(x) wangMethod(x, GO_all_MF, "MF"))

GO_pairwise_CC <- lapply(GO_Modules_CC, function(x) wangMethod(x, GO_all_CC, "CC"))

save(GO_pairwise_BP,GO_pairwise_MF,GO_pairwise_CC, file=paste0("./GO_sim/PDL1_GO_pairwise",".Rdata"))

load(file=paste0("./GO_sim/PDL1_GO_pairwise",".Rdata"))

GO_pairwise_BP <- GO_pairwise_BP[[1]]
GO_pairwise_MF <- GO_pairwise_MF[[1]]
GO_pairwise_CC <- GO_pairwise_CC[[1]]

DTI <- import('DTI201811.txt')

DTI_targets <- with(DTI, split(V2, V1))


ResGOsim <- vector("list", length(DTI_targets))

for (i in 1:length(DTI_targets)) {

  print(i)

  drug <- DTI_targets[[i]]
  name <- names(DTI_targets[i])

  GOterms_BP <- unlist(unique(Node2GO_BP[drug]))
  GOterms_MF <- unlist(unique(Node2GO_MF[drug]))
  GOterms_CC <- unlist(unique(Node2GO_CC[drug]))

  if(!length(GOterms_BP)== 0&!length(GOterms_MF)== 0&!length(GOterms_CC)== 0){

    GOSIM_BP <- combineScores(GO_pairwise_BP[, colnames(GO_pairwise_BP) %in% GOterms_BP], combine = "BMA")
    GOSIM_MF <- combineScores(GO_pairwise_MF[, colnames(GO_pairwise_MF) %in% GOterms_MF], combine = "BMA")
    GOSIM_CC <- combineScores(GO_pairwise_CC[, colnames(GO_pairwise_CC) %in% GOterms_CC], combine = "BMA")

    GOSIM <- (GOSIM_BP+GOSIM_MF+GOSIM_CC)/3
    # Get the degree bins of the target of the real drug modules
    TargetBin <- filter(Degrees, Genes %in% drug)
    # for each real drug module: select those genes in the degree bin of the target (removing the real target)
    GoodIntervalGenes <- filter(Degrees, Bin == TargetBin$Bin)$Genes %>% setdiff(TargetBin$Genes)
    #GoodIntervalGenes <- filter(Degrees, Bin %in% TargetBin$Bin)$Genes %>% setdiff(TargetBin$Genes)

    RandGOsim <- vector("list", length = 1000)
    for (j in 1:1000) {

      set.seed(j)
      # Select 100 random genes from the same bin in which the real drug target lies as random Targets
      RandTargets <- sample(GoodIntervalGenes, length(drug), replace = F)

      GOterms_BP <- unlist(unique(Node2GO_BP[RandTargets]))
      GOterms_MF <- unlist(unique(Node2GO_MF[RandTargets]))
      GOterms_CC <- unlist(unique(Node2GO_CC[RandTargets]))

      GOSIM_BP <- combineScores(GO_pairwise_BP[, colnames(GO_pairwise_BP) %in% GOterms_BP], combine = "BMA")
      GOSIM_MF <- combineScores(GO_pairwise_MF[, colnames(GO_pairwise_MF) %in% GOterms_MF], combine = "BMA")
      GOSIM_CC <- combineScores(GO_pairwise_CC[, colnames(GO_pairwise_CC) %in% GOterms_CC], combine = "BMA")

      RandGOsim[j] <- (GOSIM_BP+GOSIM_MF+GOSIM_CC)/3

    }
    permuteScore <- do.call(cbind,RandGOsim)

    permuteScore[is.na(permuteScore)] <- 0

    ## Compute the p-value based on bootstrap method
    pValue <- rowSums(abs(permuteScore) >= abs(GOSIM)) / 1000
    ## Compute the adjusted p-value. The adjusting method can be reseted
    ## (Refer to p.adjust()).
    pAdjust <- stats::p.adjust(pValue, "fdr")

    ResGOsim[[i]] <- list(GOSIM,pValue,pAdjust,name)
  }else{
    ResGOsim[[i]] <- list(0,1,1,name)

  }

}


RandGO.permuteScore <- data.frame(do.call(rbind,ResGOsim),check.rows = F,check.names = F)

names(RandGO.permuteScore) <- c('GO_score','p','fdr','name')

head(RandGO.permuteScore$GO_score)

RandGO.permuteScore$GO_score = sapply(RandGO.permuteScore$GO_score, function(x) x[1])
RandGO.permuteScore$p = sapply(RandGO.permuteScore$p, function(x) x[1])
RandGO.permuteScore$fdr = sapply(RandGO.permuteScore$fdr, function(x) x[1])
RandGO.permuteScore$name = sapply(RandGO.permuteScore$name, function(x) x[1])


sig_go <- subset(RandGO.permuteScore,GO_score>0.6&p<0.05)


write.csv(RandGO.permuteScore,file = 'RandGO.permuteScore_GO_drug.csv',quote = F,row.names = F)

save(RandGO.permuteScore,sig_go,file = './GO_sim/PDL1_all_fda_GO.Rdata')

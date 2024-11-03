
library(igraph)
library(rio)
library(tidyverse)
library(dnet)
library(xlsx)
library(clusterProfiler)
library(org.Hs.eg.db)

# PPI networks
net <- HumanNet_PI3_LCC
set <- V(net)$name
# RWR for PDL1_driver genes ====

PDL1_driver <- import('PDL1_driver_gene.txt')


driver_PDL1_gene <- PDL1_driver$V2

# define set of seeds
seeds_PDL1 <- rep(1, length(driver_PDL1_gene))
seeds_PDL1 <- data.frame(seeds_PDL1)
rownames(seeds_PDL1) <- driver_PDL1_gene

n1 <- length(intersect(rownames(seeds_PDL1), set))

# drug ---------------------------------------------------------

DTI <- import('DTI201811.txt')

length(unique(DTI$V1))

#drug overlap with PPI ---------------------------------------------------

DTI_final <- DTI[DTI$V2%in%set,]

drug_num <- length(unique(DTI_final$V1))


PTmatrix_disease <- suppressWarnings(suppressMessages(dRWR(g = net, normalise = "row",
                                                           setSeeds = seeds_PDL1, restart = 0.75,
                                                           parallel = FALSE)))

load(file = 'disease_random_rwr.Rdata')

z_score <- NULL


library(furrr)
plan(multisession, workers = availableCores()-30)

options(future.globals.maxSize= 3221222500)

system.time(
  result <- furrr::future_map(1:drug_num, function(i) {

    CID = unique(DTI_final$V1)[i]

    CID_target= unique(DTI_final[DTI_final$V1%in%CID,"V2"])

    seeds_target <- rep(1, length(CID_target))

    seeds_target <- data.frame(seeds_target)

    rownames(seeds_target) <-  CID_target

    n2 <- length(intersect(rownames(seeds_target), set))

    PTmatrix_target <- suppressWarnings(suppressMessages(dRWR(g = net, normalise = "row",
                                                              setSeeds = seeds_target, restart = 0.75,
                                                              parallel = TRUE)))

    pcc=cor(as.numeric(PTmatrix_disease), as.numeric(PTmatrix_target), method = "pearson")

      rv <- purrr::map(1:1000, function(j) {
      seeds_target_random_net <- rep(1, n2)
      seeds_target_random_net <- data.frame(seeds_target_random_net)
      row.names(seeds_target_random_net) <- set[sample(1:length(set), n2)]

      PTmatrix_target_random_net <- suppressWarnings(suppressMessages(dRWR(g = net, normalise = "row",
                                                                           setSeeds = seeds_target_random_net, restart = 0.75,
                                                                           parallel = TRUE)))

      pcc_seed_random <- cor(as.numeric(disease_random[j][[1]]), as.numeric(PTmatrix_target_random_net), method = "pearson")
      pcc_seed_random
    })


    drug_random <- do.call(cbind,rv)

    p = sum(pcc<drug_random)/1000

    mean_pcc_seed_random <- mean(drug_random)

    std_pcc_seed_random <- sd(drug_random)

    z_score <- (pcc - mean_pcc_seed_random)/std_pcc_seed_random

    res=list(CID =CID,cor=pcc,pvalue=p,Z=z_score)
    res
  }, .progress = TRUE)
)

save(result,file = 'PDL1_drug_rwr.Rdata')

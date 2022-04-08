library(tidyverse)
library(qtl)

setwd("/data/LM_qtl_manuscript/qtlMapping")

lm_cross <- read.cross("csvs", "/data/LM_qtl_manuscript", "mapConstruction/lmMap_2020-03-09.csv", 
                       "genotype_values/lm_AR1_BLUEs.csv", 
                       genotypes=c("A", "X", "B"), alleles=c("A", "B"), error.prob = .0001,
                       estimate.map = FALSE, crosstype = "riself")

lm_cross <- jittermap(lm_cross)
lm_cross <- calc.genoprob(lm_cross, step = 1, error.prob = .01)


phenoNames <- names(lm_cross$pheno)[-1]

#Get tibble of names, run control QTL mapping using Id number to generate frame...
allLods <- cim(lm_cross, pheno.col = "Id") %>% as_tibble() %>% rename("Id" = "lod")

#Get point estimates of heading date QTL with CIM, then refine with multi-qtl mapping

lodThresholds <- c()
for (pheno in phenoNames) {
  lodVals <- cim(lm_cross, pheno.col = pheno) %>%
    as_tibble()
  
  lodSig <- cim(lm_cross, pheno.col = pheno, n.perm = 1000)

  allLods <- full_join(allLods, lodVals, by = c("chr", "pos"))
  lodThresholds <- c(lodThresholds, summary(lodSig))
}

names(allLods)[-c(1:3)] <- phenoNames

#Get a QTL shortlist for each phenotype.

phenoTypes <- str_replace(phenoNames, "[Ral, Kin, Pla]+\\d+_", "") %>% unique()
allQTLPos <- NULL

for (pheno in phenoTypes) {
  lodSumByLocus <- select(allLods, chr, pos, contains(pheno)) %>% 
          pivot_longer(contains(pheno)) %>% 
          group_by(chr, pos) %>% 
          summarize(combLod = sum(value)) %>%
          filter(combLod > (2.2 * (ncol(select(allLods, contains(pheno)))))) %>% #Some evidence for QTL presence
          ungroup() %>%
          mutate(chr = as.character(chr))
  
  #Get groups
  groupVec <- c()
  curGroup <- 1
  for (row in c(2:nrow(lodSumByLocus))) {
    if (lodSumByLocus[row,1] != lodSumByLocus[row-1,1]) {
      curGroup <- curGroup + 1
    }
    
    if ((lodSumByLocus[row,2] - lodSumByLocus[row-1, 2]) > 20) { #If nearest sig. marker is greater than 10 cm away
      curGroup <- curGroup + 1 
    }
    groupVec <- c(groupVec, curGroup)
  }
  
  groupVec <- c(1, groupVec) #For lag reasons first row needs to be added
  lodSumByLocus <- add_column(lodSumByLocus, Group = groupVec)
  
  #Now get best marker for each group
  bestQTLPos <- group_by(lodSumByLocus, Group) %>%
    arrange(desc(combLod)) %>%
    slice(1)
  
  bestQTLPos <- cbind(phenoType = rep(pheno, nrow(bestQTLPos)), bestQTLPos)
  
  allQTLPos <- rbind(allQTLPos, bestQTLPos)
}


for (pheno in phenoTypes) {
  phenoQTL <- filter(allQTLPos, phenoType == pheno)
  phenoMod <- makeqtl(testSim, chr = phenoQTL$chr, pos = phenoQTL$pos, qtl.name = paste0(phenoQTL$chr, as.character(phenoQTL$Group)))
  
  qtlNum <- c(1:nrow(filter(allQTLPos, phenoType == pheno))) %>% as.character()
  qtlNum <- paste0("Q", qtlNum) %>% paste(collapse = " + ")
  qtlNum <- paste0("y ~ ", qtlNum)
  
  refQTL <- 
  
  print(phenoMod)
  print(qtlNum)
  
  qtlMod <- fitqtl(testSim, pheno.col = "Pla19_HD", qtl = phenoMod, formula = qtlNum)
  print(summary(qtlMod))
}



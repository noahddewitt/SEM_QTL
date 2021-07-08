library(qtl)

#Need to at least read in pheno data for filtering of genotype
rilVals <- read.csv("lm_AR1_BLUEs.csv", row.names = 1)

#####Get genome-wide markers for testing
lm_cross <- read.cross("csvs", ".", "lmMap_2020-03-09.csv", 
                       "lm_AR1_BLUEs.csv", 
                       genotypes=c("A", "X", "B"), alleles=c("A", "B"), error.prob = .0001,
                       estimate.map = FALSE, crosstype = "riself")

lm_cross <- jittermap(lm_cross)
lm_cross <- calc.genoprob(lm_cross, step = 1, error.prob = .01)

lm_crossS <- sim.geno(lm_cross, step = 2, n.draws = 500, err = .001)
rm(lm_cross)
#Bind GBS-derived genome-wide genotype probabilities to data frame:
subVals_cross <- subset(lm_crossS, ind = rownames(rilVals))
rm(lm_crossS)
sub_genoProbs <- pull.genoprob(subVals_cross)
rm(subVals_cross)

sub_genoProbs <- sub_genoProbs[,grep("BB", colnames(sub_genoProbs))]

#We also want to replace combined markers and the :BB as lavaan doesn't like colons
colnames(sub_genoProbs) <- gsub(":.*", "", colnames(sub_genoProbs))
sub_genoProbs <- scale(sub_genoProbs)

write.csv(sub_genoProbs, "scaled_genoProbs.csv")

#####Get phenotype data for this environment and scale/filter
rilVals <- rilVals[rownames(sub_genoProbs), ] #Order and toss parents
rilVals <- rilVals[, grep("Kin19", colnames(rilVals))] #Just get our environment
colnames(rilVals) <- gsub("Kin19_", "", colnames(rilVals))
rilVals <- scale(rilVals)

write.csv(rilVals, "Kin19_rilPhenoVals.csv")

#Many plant growth effects here also have direct YC effects... this should come out organically in the model. 

qNames <- c(rep(NA, 6))
combN <- c("5", "1", rep(NA, 4))
spsStr <- c("SPS", "HD", "PH", rep(NA, 3))
kpsStr <- c("KPS", "PH", "rustScore", rep(NA, 3))
kwStr <- c("KWidth", "HD", "PH", rep(NA, 3))
klStr <- c("KLength", "rustScore", rep(NA, 4))
tkwStr <- c("TKW", "KArea", "KWidth", "rustScore", "HD", "PH")

baseGraph <- rbind(qNames, combN, spsStr, kpsStr, kwStr, klStr, tkwStr)

write.csv(baseGraph, "Kin19_baseGraph.csv")





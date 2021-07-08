library(parallel)
library(lavaan)

source("graphConstructionFunctions.R")

args = commandArgs(trailingOnly=TRUE)

runStr <- as.character(args[1])
curCores <- as.numeric(args[2])
baseStr <- "Kin19_"
baseGraph <- read.csv(paste0("allGraphs/", baseStr, runStr, ".csv"), colClasses =  "character")

graphVec <- c()

#Turn graph df into vector of strings for yield components
for (i in c(3:nrow(baseGraph))) { #First and last row contain different info
  rhsStr <- paste(baseGraph[i,-c(1:2)], collapse = "+")
  rhsStr <- gsub("\\+NA", "", rhsStr)
  nodeStr <- paste0(" \n",baseGraph[[i,2]], "~", rhsStr)
  graphVec <- c(graphVec, nodeStr)
}

#Modify TKW because it depends on KArea
graphVec[5] <- paste0(" \nKArea~KLength+KWidth", graphVec[5])

#Only grab up to ten QTL per phenotype for iterations...
iterGraphVec <- c()
#Turn graph df into vector of strings for yield components
for (i in c(3:nrow(baseGraph))) { 
  rowLength <- sum(!is.na(baseGraph[i,]))
  nonQTL <- sum(baseGraph[i,-c(1:2)] %in% c("KArea", "KWidth", "PH", "HD", "rustScore")) 
  if (rowLength > (7 + nonQTL)) {
    rowLength <- (7 + nonQTL)
  }
  rhsStr <- paste(baseGraph[i,c(3:rowLength)], collapse = "+")
  rhsStr <- gsub("\\+NA", "", rhsStr)
  nodeStr <- paste0(" \n",baseGraph[[i,2]], "~", rhsStr)
  iterGraphVec <- c(iterGraphVec, nodeStr)
}

#Modify TW because it depends on KArea
iterGraphVec[5] <- paste0(" \nKArea~KLength+KWidth", iterGraphVec[5])

#Read in phenotype and genotype data
phenoDF <- read.csv("Kin19_rilPhenoVals.csv", row.names = 1)
genoDF <- read.csv("scaled_genoProbs.csv", row.names = 1)

maxTime <- 100 #hrs
startTime <- Sys.time()
combnVec <- as.numeric(unlist(baseGraph[2,c(2,3)]))

baseGraphVec <- graphVec

while(as.numeric(difftime(Sys.time(), startTime, units='hours')) < maxTime) {
	#Creating this iteration's current testing arcs
	arcsVec <- combn(c(1:5), combnVec[1])[,combnVec[2]]
	subGraphVec <- graphVec[arcsVec]
	print("Starting iterations")
	sigThreshold <- setLODThreshold(iterGraphVec[arcsVec], phenoDF, genoDF, iter = 1000)
	print("Threshold")
	print(sigThreshold)

	newSubGraphVec <- graphIteration(subGraphVec, genoDF, phenoDF, sigThreshold)

	#Update full graph with recovered subgraph
 	for (i in c(1:length(newSubGraphVec))) {
		graphVec[arcsVec[i]] <- newSubGraphVec[i]
	}

	#Update combnVec
	maxComb <- 120 / (factorial(combnVec[1]) * factorial(5-combnVec[1]))
	if (combnVec[2] < maxComb) {
		combnVec[2] <- combnVec[2] + 1
	} else {
		combnVec <- c(combnVec[1] - 1, 1)
		if (combnVec[1] == 0) {
			break #You're done!
		}
	}

	#Write out updated graph -- name of file = NEXT combn it will do
	writeGraph(baseGraph, graphVec, combnVec, baseStr) 	

}

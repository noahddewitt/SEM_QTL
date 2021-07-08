testMarkerPerm <- function(index, jointDF, varIndices, altModStr, nulModStr) {
  jointDF <- jointDF[,unique(c(varIndices, index), fromLast = T)]
  colnames(jointDF)[ncol(jointDF)] <- "M"
  return(as.numeric(lavaan::logLik(lavaan::sem(altModStr, data = jointDF, fixed.x = F))) - as.numeric(lavaan::logLik(lavaan::sem(nulModStr, data = jointDF, fixed.x = F))))
}

testMarkerAltPerm <- function(index, jointDF, varIndices, altModStr) {
  jointDF <- jointDF[,unique(c(varIndices, index), fromLast = T)]
  colnames(jointDF)[ncol(jointDF)] <- "M"
  return(as.numeric(lavaan::logLik(lavaan::sem(altModStr, data = jointDF, fixed.x = F))))
}

testMarker <- function(index, jointDF, varIndices, altModStr, nulModStr) {
  jointDF <- jointDF[,unique(c(varIndices, index), fromLast = T)] #This is a little bit more work but it's really not going to be the limiting factor once we hvae a lot of QTL in there. 
  #FromLast is important to make sure if already in tehre the correct column renamed.

  #This should only be an issue for this specific data set.
  if (index %in% varIndices) {
    return(0)
  } else {
    colnames(jointDF)[ncol(jointDF)] <- "M"
    return(as.numeric(try(lavaan::logLik(lavaan::sem(altModStr, data = jointDF, fixed.x = F)), silent = T)) - as.numeric(try(lavaan::logLik(lavaan::sem(nulModStr, data = jointDF, fixed.x = F))), silent = T))
  }
}

setLODThreshold <- function(subGraph, phenoDF, genoDF, curCores = 72, iter = 1000, perNum = 5, propScan = .1) {
	qNames <- c()
	for (trait in subGraph) {
	  trait <- regmatches(trait, gregexpr(" \\n[^ ]+$", trait))[[1]] #Only gets last  trait in multi-trait strings
	  traitQs <- regmatches(trait,gregexpr("[\\+|~][A-Za-z0-9\\.\\_]+",trait))[[1]]
	  traitQs <- gsub("^[\\+|~]", "", traitQs)
	  qNames <- c(qNames, traitQs)
	}
	qNames <- unique(qNames)
	qNames <- qNames[!qNames %in% c("rustScore", "HD", "PH", "KArea", "KWidth")]

	#Here, we want known factors to not vary, and unknown factors to vary
	phenoDF <- cbind(phenoDF, genoDF[,qNames])
	genoDF <- genoDF[,!colnames(genoDF) %in% qNames]

	varIndices <- c(1:ncol(phenoDF))

	nulModStr <- paste0(paste(c(subGraph, ""), collapse = " "), "\nM~1")
	altModStr <- paste0(paste(c(subGraph, ""), collapse = "+M "), "\nM~1")
	print(altModStr)

	storeLODs <- c()

	while(length(storeLODs) <= iter) {
		ranGenoDF <- genoDF[sample(c(1:nrow(genoDF)),  nrow(genoDF), replace  = F),]
		jointDF <- cbind(phenoDF, ranGenoDF)

		if ((length(storeLODs) %% 100) == 0) {
			print(paste0("Starting Iteration ", as.character(length(storeLODs))))
		}


                if ((length(storeLODs) %% 10) ==  0)  {
                        write.csv(storeLODs, "allGraphs/curIterLODs.csv")
                }

		#Sample every Nth index
		sparseIndices <- c(ncol(phenoDF) + (1:(round(ncol(ranGenoDF)/perNum)))*perNum)

		sparseLods <-  mclapply(sparseIndices, testMarkerAltPerm, jointDF, varIndices, altModStr, mc.cores = (curCores - 1))
		sparseLods <- unlist(sparseLods)

		#Divide into bins of 100 TOTAL markers (incl. hidden, so if fifths then 20 per bin)
		binNum <- (length(sparseLods) / (100/perNum))

		maxBinLods <- c()
		for (y in c(1:binNum)) {
		  binPos <- y * (100 / perNum)
		  binRange <- c((binPos - (50/perNum)):(binPos + (50/perNum)))
		  maxBinLod <- which(sparseLods[binRange] == max(sparseLods[binRange]))
		  maxBinLods <- c(maxBinLods, binRange[maxBinLod]) #maxBinLod is the index relative to binRange vector, need to get actual index
		}

		topSparseLods <- maxBinLods[order(sparseLods[maxBinLods], decreasing = T)][1:round(length(maxBinLods) / 10)]
		topSparseLods <- topSparseLods * perNum #Return to full scale

		maxLOD <- 0
		allPosSweeps <- c()
		for (pos in topSparseLods) {
		  sweepPos <- c((pos-4):(pos+4))
		  allPosSweeps <- c(allPosSweeps, sweepPos)
		}
		
		#Start in the randomized marker portion of the DF
		allPosSweeps <- allPosSweeps[order(allPosSweeps)] + ncol(phenoDF) 

		allSweepLODs <- mclapply(allPosSweeps, testMarkerPerm, jointDF, varIndices, altModStr, nulModStr, mc.cores = (curCores - 1))
		storeLODs <- c(storeLODs, max(unlist(allSweepLODs)))
	}

	storeLODs <- storeLODs[order(storeLODs)]
	sigCutOff <- storeLODs[round(length(storeLODs) * .95)]
	return(sigCutOff)
}

graphIteration <- function(subGraph, genoDF, phenoDF, subGraphCutoff, curCores = 72) {
	#We're going to convert the strings in the subGraph vector into vectors of "QTL" (which include phenotypes) stored in list
	qtlList <- list()
	for (trait in subGraph) {
	  traitName <- regmatches(trait,gregexpr("^ \\n.+~",trait))[[1]]

	  #Removes \n and ~
	  traitName <- substr(traitName, 3, nchar(traitName)-1) 

	  #Only gets last  trait in multi-trait strings
	  trait <- regmatches(trait, gregexpr(" \\n[^ ]+$", trait))[[1]] 	  traitQs <- regmatches(trait,gregexpr("[\\+|~][A-Za-z0-9\\.\\_]+",trait))[[1]]
	  traitQs <- gsub("^[\\+|~]", "", traitQs)
	  qtlList[[traitName]] <- traitQs
	}
	qNames <- unique(unlist(qtlList))
	qNames <- qNames[!qNames %in% c("rustScore", "HD", "PH", "KArea", "KLength", "KWidth")]

	jointDF <- cbind(phenoDF, genoDF)

	noMoreSigQTL = FALSE

	while(noMoreSigQTL == F){
		#In giant df, which are the ones we want to always keep?
		varIndices <- c(1:ncol(jointDF))[colnames(jointDF) %in% c(colnames(phenoDF), qNames)] 

		nulModStr <- "M~1"
		altModStr <- "M~1"
		for (trait in rev(names(qtlList))) {
		  nulModTraitStr <- paste0(trait, "~", paste0(qtlList[[trait]], collapse = "+"))
		  nulModStr <- paste0(nulModTraitStr, " \n", nulModStr)
		  
		  altModTraitStr <- paste0(trait, "~", paste0(c(qtlList[[trait]], "M"), collapse = "+"))
		  altModStr <- paste0(altModTraitStr, " \n", altModStr)
		}
		
		
		testMarkers <- c((ncol(phenoDF) + 1):ncol(jointDF))
		testMarkersExclude <- which(grepl("6B", colnames(jointDF)))
		testMarkers <- testMarkers[!testMarkers %in% testMarkersExclude]

		allLODs <- mclapply(testMarkers, testMarker, jointDF, varIndices, altModStr, nulModStr, mc.cores = (curCores-1))
		allLODs <- unlist(allLODs)

		#This happens sometimes with markers linked to QTL
		allLODs[is.na(allLODs)] <- 0 		names(allLODs) <- colnames(jointDF)[testMarkers]

		print("Iteration Done") 
		print(max(allLODs))

		if (typeof(max(allLODs)) == "character") {
		  write.csv(allLODs, "Kin19_graph2_brokenRun.csv")
		  print("broken run")
		  noMoreSigQTL <- T
		} else if(max(allLODs) > subGraphCutoff) {
		  newQ <- names(allLODs[order(allLODs, decreasing = T)])[1]
		  print(newQ)

		  #Get correlations of known markers with new QTL
		  corVec <- cor(genoDF[,c(qNames, newQ)])[,newQ] 
		  corVec <- corVec[-length(corVec)] #Corr'ns just b/w new QTL & others
		  maxCor <- max(corVec)

		  if(maxCor > .95) {
			newQIndexJoint <- which(colnames(jointDF) == newQ)
			newQIndexGeno <- which(colnames(genoDF) == newQ)

			jointDF <- jointDF[,-newQIndexJoint] 
			genoDF <- genoDF[,-newQIndexGeno]
			next
		  } else {
			qNames <- c(qNames, newQ)
		  }

		  for (trait in names(qtlList)) {
			qtlList[[trait]] <- c(qtlList[[trait]], newQ)
		  }

		  #Here's where we remove individual arcs that aren't significant
		  allGoodArcs <- FALSE
		  while (allGoodArcs == FALSE) {
			curModStr <- ""
			for (trait in names(qtlList)) {
			  curModTraitStr <- paste0(trait, "~", paste0(qtlList[[trait]], collapse = "+"))
			  curModStr <- paste0(curModStr, " \n", curModTraitStr) 
			}

			print(curModStr)
			curMod <- tryCatch(
      					sem(curModStr, data = jointDF, fixed.x = F),
      					error=function(e) e
  				  )
			
  			if(inherits(curMod, "error")){  
  			        newQIndexJoint <- which(colnames(jointDF) == newQ)
                        	newQIndexGeno <- which(colnames(genoDF) == newQ)
                	        jointDF <- jointDF[,-newQIndexJoint] 
        	                genoDF <- genoDF[,-newQIndexGeno]
	                        next

			}
		
			curModSum <- summary(curMod)

			qInfo <- curModSum$PE[(curModSum$PE$rhs == newQ) & (curModSum$PE$lhs %in% c("SPS", "KPS", "KWidth", "KLength", "TKW")),]

			#get arc w/ highest p value
			qWorstArc <- qInfo[order(qInfo$pvalue, decreasing = T)[1],] 			if (qWorstArc[[8]] > .1) {
			  
			  ycName <- qWorstArc[[1]]
			  if (ycName == "TKW") {ycName <- "KArea~KLength+KWidth \nTKW"}
			  qtlList[[ycName]] <- qtlList[[ycName]][qtlList[[ycName]] != newQ] 

			} else {
			  allGoodArcs <- TRUE
			  print("Done removing arcs!")
			}
		  }
		  allQInfo <- curModSum$PE[curModSum$PE$lhs %in% c("SPS", "KPS", "KWidth", "KLength", "TKW") & curModSum$PE$op == "~",]
		  allWorstArc <- allQInfo[order(allQInfo$pvalue, decreasing = T)[1],]
		  
		  write.table(t(c(newQ, as.character(Sys.time()), paste0(subGraph, collapse = "_"), gsub("\\n", "", nulModStr), gsub("\\n", "", altModStr), max(allLODs), subGraphCutoff, allLODs)), 
				"Kin19_graph2_iterLODs.csv", sep = ",", append = T, col.names = F)

		  if (allWorstArc[[8]] > .1) {
			ycName <- allWorstArc[[1]]
			if (ycName == "TKW") {ycName <- "KArea~KLength+KWidth \nTKW"}
			
			badQ <- allWorstArc[[3]]
			
			print("No longer significant QTL or phenotype")
			print(ycName)
			print(badQ) 
			#Remove qtl from vector of QTL associated with given trait
			qtlList[[ycName]] <- qtlList[[ycName]][qtlList[[ycName]] != badQ] 
		  }
		} else {
			noMoreSigQTL <- TRUE
			print("No more sig. qtl Final LOD:")
			print(max(allLODs))
		}
	}

	#Convert back from list object to vector of string
	newSubGraph <- c()	for (trait in names(qtlList)) {
		curModTraitStr <- paste0(" \n",trait, "~", paste0(qtlList[[trait]], collapse = "+"))
		newSubGraph <- c(newSubGraph, curModTraitStr)
	}
	return(newSubGraph) 
}

writeGraph <- function(baseGraph, graphVec, combnVec, baseStr) {
	qtlList <- list()
	for (trait in graphVec) {
	  trait <- regmatches(trait, gregexpr(" \\n[^ ]+$", trait))[[1]] 
	  traitName <- regmatches(trait,gregexpr("^ \\n.+~",trait))[[1]]
	  traitName <- substr(traitName, 3, nchar(traitName)-1)

	  traitQs <- regmatches(trait,gregexpr("[\\+|~][A-Za-z0-9\\.\\_]+",trait))[[1]]
	  traitQs <- gsub("^[\\+|~]", "", traitQs)
	  qtlList[[traitName]] <- traitQs
	}
	qNames <- unique(unlist(qtlList))
	qNames <- qNames[!qNames %in% c("rustScore", "HD", "PH", "KArea", "KWidth")]

	longestVecLength <- max(c(length(qNames), length(combnVec), lengths(qtlList)))
	length(qNames) <- longestVecLength
	length(combnVec) <- longestVecLength

	newGraph <- rbind(qNames, combnVec)
	for (trait in names(qtlList)) {
		newRow <- c(trait, qtlList[[trait]])
		length(newRow) <- longestVecLength
		newGraph <- rbind(newGraph, newRow)
	}

	newGraph <- cbind(baseGraph[,1], newGraph)

	graphStr <- paste0("allGraphs/", baseStr, paste0(as.character(combnVec[1:2]), collapse = ""), ".csv")
	write.csv(newGraph, graphStr, row.names = F)
}

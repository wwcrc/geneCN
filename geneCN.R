#!/usr/bin/env Rscript

# Copy number analysis script
# 11/11/2015
# Copyright (C) 2015 14MG, 2017-2018 University of Glasgow
# Author: Susie Cooke
# Version 2.0.2
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

########################## Get Command Line Options ##################################

args <- commandArgs(trailingOnly = TRUE) # Read in command line parameters
data <- read.csv(args[1], header = FALSE) # Output of perl script
static <- read.table(args[2], sep = "\t") # Static file of annotated bin regions
SampleID <- args[3] # Sample name, will appear in output file names and headers
GenomeBuild <- args[4] # Reference genome version
genders <- c('female', 'male', 'unknown')

if (length(args) < 4 || length(args) > 6) {
  stop("The following arguments are needed: data file, regions file, sample name, genome build, (thresholds file), (sample gender)", call.=FALSE)
} else if (length(args) == 4) {
  thresholds <- 'empty'
  gender <- 'unknown'
} else if (length(args) == 5) {
  thresholds <- read.table(args[5], sep = "\t", header = TRUE) # Static file of thresholds for calling
  gender <- 'unknown'
} else if (length(args) == 6) {
  thresholds <- read.table(args[5], sep = "\t", header = TRUE) # Static file of thresholds for calling
  gender <- args[6]
}

if (!gender %in% c(genders)) {
  stop("Specified gender must be female/male", call.=FALSE)
}
rm(genders)

########################## Set Global Parameters #####################################

# Default y axis height for whole genome plots
ymax <- 4
ymin <- -4
myDataMax <- 8
myDataMin <- -6

mySpecialCases <- c('background', 'other')

# Colours for whole genome plots
colour1 <- 'deepskyblue'
colour2 <- 'darkslateblue'

# Chromosomes to be plotted
Chrs <- unique(static$V1)

# Phasing of colours for data points
mycol1Chrs <- Chrs[c(TRUE, FALSE)]
mycol2Chrs <- Chrs[c(FALSE, TRUE)]

###################################### Functions ###################################

# Function to add a bit of noise to zero values to prevent errors
increaseZero <- function(x) {
  minValue <- 100
  for (value in x) {
    if (value < minValue && value > 0) {
      minValue <- value
    }
  }
  rm(value)
  myNewValues <- c()
  for (value in x) {
    if (value == 0) {
      myNewValues <- c(myNewValues, minValue/10)
    } else {myNewValues <- c(myNewValues, value) }
  }
  rm(value, minValue)
  return(myNewValues)
}

# Function to cap data prior to plotting
capped <- function(x, Min, Max) {
  CappedData <- c() 
  for (value in x) {
    if (is.na(value)) {
      CappedData <- c(CappedData, value)
      next
    } else if (value < Min) {
     CappedData <- c(CappedData, Min) #Need data to be within plot window
   } else if (value > Max) {
      CappedData <- c(CappedData, Max) #Need data to be within plot window
   } else {
     CappedData <- c(CappedData, value)
   }
  }
  return(CappedData)
}

# Function to calculate y-axis limits
limit <- function(x, Min, Max) {
  if (max(x, na.rm = TRUE) > Max) { 
    mySampleMax <- max(x, na.rm = TRUE) #Extend axis to include all data
  } else {
    mySampleMax <- Max #Need axes to be sensible heights
  }
  if (min(x, na.rm = TRUE) < Min) {
    mySampleMin <- min(x, na.rm = TRUE) #Extend axis to include all data
  } else {
    mySampleMin <- Min #Need axes to be sensible heights
  }
  return(c(mySampleMin, mySampleMax))
}

################################### Normalisation #####################################

normed_data <- data[,4]/static[,6]
data[,4] <- normed_data

# Remove zero values to prevent -Inf after logging
newData <- increaseZero(data[,4])
data[,4] <- newData

# Loess normalise
GCloess <- loess(data[,4] ~ static[,4])
myPredictions <- predict(GCloess, data.frame(static[,4]))

# Internal normalisation
med <- median(data[,4]/myPredictions)
dataNorm <- log2(data[,4]/myPredictions/med)
data[,4] <- dataNorm

# Tidy up
rm(med, GCloess, myPredictions, dataNorm, newData)

############################ Generate whole genome plots ###########################################

# Calculate chromosome offsets within the genome
chrLengths <- c()
myOffsets <- c(0)
increment <- 0

for(chromosome in Chrs) {
  numberBins <- length(static[static[,1] == chromosome, 1])
  chrLengths <- c(chrLengths, numberBins) 
}
myDf <- data.frame(Chr = Chrs, Length = chrLengths)

for (bins in myDf$Length) {
  myOffsets <- c(myOffsets, bins+increment)
  increment <- bins+increment
}

myDf$Offsets <- myOffsets[1:nrow(myDf)]
rm(bins, increment, chrLengths, myOffsets, numberBins, chromosome)

# Cap data prior to plotting
myCappedData <- capped(data[,4], myDataMin, myDataMax) 

# Calculate y-axis limits
myYlimits <- limit(myCappedData, ymin, ymax)

# Generate plots as .pdf
pdf(paste(SampleID, '_CNplot.pdf', sep = ""), width = 30)

# Set up plot area but don't add any data
plot(data[,4], ylab = 'Log2 Ratio of Normalised Depths', ylim = myYlimits, 
     xlab = 'Bin Index', cex = 0.5, pch = 19, type = "n", axes = TRUE, main = SampleID)
# Draw different chromosome data points in different colours
for(chr in Chrs) {
    myStart <- myDf[myDf$Chr == chr, "Offsets"] + 1
    myEnd <- myDf[myDf$Chr == chr, "Offsets"] + length(data[data[,1] == chr, 1])
    x <- c(myStart : myEnd)
    y <- myCappedData[data[,1] == chr]
    if (is.element(chr, mycol1Chrs)) {
      myColour <- colour1
    } else if (is.element(chr, mycol2Chrs)) {
      myColour <- colour2
    }
    points(x, y, col = myColour, cex = 0.5, pch = 19)
  }
rm(x, y, myStart, myEnd, chr)

dev.off()

################################## Chromosome by Chromosome Plots ##################################

# Colours for gene features
featurecolour1 <- 'red'
featurecolour2 <- 'deepskyblue'

for (chr in Chrs) {
  pdf(paste(SampleID, chr, 'plot.pdf', sep = "_"), width = 20)
  myFeatureList <- unique(static$V5[static$V1 == chr])
  myFeatureList <- myFeatureList[!myFeatureList %in% mySpecialCases]
  col1features <- myFeatureList[c(TRUE, FALSE)]
  col2features <- myFeatureList[c(FALSE, TRUE)]
  plot(data[,4][data$V1 == chr], pch = 19, cex = 0.5, ylim = myYlimits, main = paste(SampleID, '_', chr), ylab = 'Log2 Ratio of Normalised Depths', xlab = 'Bin Index', type = "n", axes = TRUE)
  for (feature in unique(static$V5[static$V1 == chr])) {
    if (feature %in% mySpecialCases) {
      points(which(static$V5[static$V1 == chr] == feature), data[,4][static$V1 == chr][which(static$V5[static$V1 == chr] == feature)], pch = 19)
    } else if (is.element(feature, col1features)){
      points(which(static$V5[static$V1 == chr] == feature), data[,4][static$V1 == chr][which(static$V5[static$V1 == chr] == feature)], pch = 19, col = featurecolour1)
    } else if (is.element(feature, col2features)){
      points(which(static$V5[static$V1 == chr] == feature), data[,4][static$V1 == chr][which(static$V5[static$V1 == chr] == feature)], pch = 19, col = featurecolour2)
    }
  }
  dev.off()
}
rm(featurecolour1, featurecolour2, myFeatureList, col1features, col2features, feature, chr)

################################### Gene CN state Calling ######################################

# Pull out chr, lowest coord and highest coord for each feature
myFeatureChr <- list()
myFeatureL1 <- list()
myFeatureH1 <- list()

for(feature in levels(static[,5])) {
  myFeatureChr[[feature]] <- droplevels(unique(static[static[,5] == feature , 1]))
  myFeatureL1[[feature]] <- min(static[static[,5] == feature, 2])
  myFeatureH1[[feature]] <- max(static[static[,5] == feature, 3])
}

# Calculate QC metric
mySDs <- c()
for(feature in levels(static[,5])) {
  mySDs <- c(mySDs, sd(data[static[,5] == feature, 4], na.rm = TRUE))
}
myQCval <- median(mySDs)
rm(mySDs, feature)

# Thresholds
myPvalCutoff <- 0.01

# Define reference distribution
myBackground <- data[static[,5] == 'background', 4]

# Initiate output file
myOutput <- paste(SampleID, 'CNcalls.txt', sep = "_") 
sink(myOutput) #Start writing to output file

#Print header lines
cat('##fileDate=',format(Sys.Date()),"\n", 
    '##CollaboratorSampleID=',SampleID,"\n",
    '##GenomeBuild=',GenomeBuild,"\n",
    '##QCvalue=', myQCval,"\n",
    '##Chr="The chromosome of the feature"',"\n",
    '##NA="Column not in use"',"\n",
    '##L0="zero-based coordinate for the start of the feature"',"\n",
    '##L1="one-based coordinate for the start of the feature"',"\n",
    '##Feature="The name of the feature"',"\n",
    '##State="The state of the feature if thresholds were provided else the pvalue"',"\n",
    '##H0="zero-based coordinate for the end of the feature"',"\n",
    '##H1="one-based coordinate for the end of the feature"',"\n",
    '##Score="Score"',"\n",
    sep="")
cat('#Chr','NA','L0','L1','Feature','State', 'Chr','NA','H0','H1','Feature','State','Score', sep="\t")
cat("\n")
# Test features and write results to file
for(feature in levels(static[,5])) {
  if (!feature %in% mySpecialCases) {
    myGene <- data[static[,5] == feature , 4]
    myResult <- t.test(myBackground, myGene)
    myMeanDiff <- mean(myGene, na.rm = TRUE) - mean(myBackground, na.rm = TRUE)
    if (thresholds == 'empty' || gender == 'unknown') {
      cat(levels(myFeatureChr[[feature]]),'.',myFeatureL1[[feature]] -1,myFeatureL1[[feature]],feature,myResult$p.value
          ,levels(myFeatureChr[[feature]]),'.',myFeatureH1[[feature]] -1,myFeatureH1[[feature]],feature,myResult$p.value,myMeanDiff, sep = "\t")
      cat("\n")
    }
    if (gender == 'female') {
      if ((myMeanDiff > thresholds$upper_limit_female[thresholds$Feature == feature]) && (myResult$p.value < myPvalCutoff)) {
        myState <- 'gain'
        cat(levels(myFeatureChr[[feature]]),'.',myFeatureL1[[feature]] -1,myFeatureL1[[feature]],feature,myState
            ,levels(myFeatureChr[[feature]]),'.',myFeatureH1[[feature]] -1,myFeatureH1[[feature]],feature,myState,myMeanDiff, sep = "\t")
        cat("\n")
      }
      else if ((myMeanDiff < thresholds$lower_limit_female[thresholds$Feature == feature]) && (myResult$p.value < myPvalCutoff)) {
        myState <- 'loss'
        cat(levels(myFeatureChr[[feature]]),'.',myFeatureL1[[feature]] -1,myFeatureL1[[feature]],feature,myState
            ,levels(myFeatureChr[[feature]]),'.',myFeatureH1[[feature]] -1,myFeatureH1[[feature]],feature,myState,myMeanDiff, sep = "\t")
        cat("\n")
      }
    } else if (gender == 'male') {
      if ((myMeanDiff > thresholds$upper_limit_male[thresholds$Feature == feature]) && (myResult$p.value < myPvalCutoff)) {
        myState <- 'gain'
        cat(levels(myFeatureChr[[feature]]),'.',myFeatureL1[[feature]] -1,myFeatureL1[[feature]],feature,myState
            ,levels(myFeatureChr[[feature]]),'.',myFeatureH1[[feature]] -1,myFeatureH1[[feature]],feature,myState,myMeanDiff, sep = "\t")
        cat("\n")
      }
      else if ((myMeanDiff < thresholds$lower_limit_male[thresholds$Feature == feature]) && (myResult$p.value < myPvalCutoff)) {
        myState <- 'loss'
        cat(levels(myFeatureChr[[feature]]),'.',myFeatureL1[[feature]] -1,myFeatureL1[[feature]],feature,myState
          ,levels(myFeatureChr[[feature]]),'.',myFeatureH1[[feature]] -1,myFeatureH1[[feature]],feature,myState,myMeanDiff, sep = "\t")
        cat("\n")
      }
    } 
  }
}
sink()

########################################## END ##############################################

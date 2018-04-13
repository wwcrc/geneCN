#!/usr/bin/env Rscript

# Copy number analysis script
# 11/11/2015
# Copyright (C) 2015 14MG
# Copyright (C) 2017-2018 University of Glasgow
# Author: Susie Cooke
# Version 1.0
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

################################### SETUP ############################################
########################## Get Command Line Options ##################################

args <- commandArgs(trailingOnly = TRUE) # Read in command line parameters
data <- read.csv(args[1], header = FALSE) # Output of perl script
static <- read.table(args[2], sep = "\t") # Static file of annotated bin regions
SampleID <- args[3]

########################## Set Global Parameters #####################################

# Default y axis height for whole genome plots
ymax <- 4
ymin <- -4
myDataMax <- 8
myDataMin <- -6

# Colours for whole genome plots
colour1 <- 'deepskyblue'
colour2 <- 'darkslateblue'

# Chromosomes to be plotted
Chrs <- c(levels(static$V1))

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
    } else if (is.element(chr, mycol3Chrs)) {
      myColour <- colour3
    }
    points(x, y, col = myColour, cex = 0.5, pch = 19)
  }
  rm(x, y, myStart, myEnd, chr)

dev.off()

############################ Output some values about the data ##################################

myOutput <- paste(SampleID, '_values.txt', sep = "")

sink(myOutput) #Start writing to output file
cat(SampleID, "\t", mean(data$V4[grep("[^xy]", tolower(static$V1))]), "\t", sd(data$V4[grep("[^xy]", tolower(static$V1))]))
cat("\t", max(data$V4[grep("[^xy]", tolower(static$V1))]), "\t", min(data$V4[grep("[^xy]", tolower(static$V1))]),"\n")
sink()

########################################## END ##############################################

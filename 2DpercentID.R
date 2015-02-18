#!/usr/bin/env Rscript

# Figure 2

# load ggplot2 for creating a nice graph
library(ggplot2)

# The arguments from the calling script are passed in using the commandArgs function
# these arguments include the path to use as well as the name of the strain
cmdArgs <- commandArgs(trailingOnly=TRUE)

# get the path
path <- getwd()
setwd(path)

# load the file of interest
fileofinterest <- cmdArgs[1]

# Load the data
graph <- read.table(fileofinterest, header=TRUE)

# Initialise the required variables
siftedResults <- vector()
fCOut <- vector()
rLOut <- vector()
percentIDOut <- vector()
SDIDOut <- vector()
foldCovOut <- vector()
foldCovSDOut <- vector()
results <- vector()

# Populate the results vector with the appropriate data
# 'AspA' is occasionally used as a negative control, and should not be included in 'results'
results <- graph[graph$target != 'aspA', ]

# find the unique fold coverages present in the study
uniqueFoldCoverages <- unique(results$foldCoverage)

# readlength could be taken from the dataset, but, instead, it's set manually
# previously, when different read lengths were being compared, there was a 
# uniqueReadLengths variable used
readlength <- 21

for (foldcoverage in uniqueFoldCoverages) {
  # for each iteration, sifted results only contains data of the current fold coverage
  siftedResults <- results[results$readLength == readlength & results$foldCoverage == foldcoverage, ]
  # populate the vectors
  rLOut <- append(rLOut, readlength, after = length(rLOut))
  fCOut <- append(fCOut, foldcoverage, after = length(fCOut))
  percentIDOut <- append(percentIDOut, mean(siftedResults$MedianPercentID), after = length(percentIDOut))
  SDIDOut <- append(SDIDOut, sd(siftedResults$MedianPercentID), after = length(SDIDOut))
  foldCovOut <- append(foldCovOut, mean(siftedResults$MedianFoldCoverage), after = length(foldCovOut))
  foldCovSDOut <- append(foldCovSDOut, mean(siftedResults$FoldCoverageSD), after = length(foldCovSDOut))
}

# populate the data frame to be used in generating the figure
output <- data.frame(readLength = rLOut, foldCoverage = fCOut, meanpercentID = percentIDOut, SDPercentID = SDIDOut, meanFoldCoverage = foldCovOut, foldCovSD = foldCovSDOut)

# output is printed to console - allows visual inspection
print(output)

# plot a nice output graph using ggplot2
p <- ggplot(output, aes(x=foldCoverage,y=percentIDOut)) 
# add scatterplot and a simple theme
p + geom_point() + theme_bw() + theme(axis.text=element_text(size=10)) +
  # error bars
  geom_errorbar(aes(ymin=percentIDOut-SDPercentID, ymax=percentIDOut+SDPercentID)) +
  # axis labels
  xlab("Fold-Coverage") + ylab("Percent Identity") + ylim(-5, 115)

filename <- "Figure2_multipanel_7by19.pdf"

# saves the chart as 'filename', in with dimensions provided below
ggsave(filename, width=19, height=7, units=c("cm"))

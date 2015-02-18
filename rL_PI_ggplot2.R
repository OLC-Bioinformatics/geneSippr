#!/usr/bin/env Rscript

# Figure 1

# load ggplot2 for creating a nice graph
library(ggplot2)

# load the file of interest
# hard-coded for now
path <- "/home/blais/Desktop/"
setwd(path)
# Files are hard-coded, too
fileofinterest <- "SipprModelling_1418048737.49.csv"

# load the data
graph <- read.table(fileofinterest, header=TRUE)
# Initialise the required variables
siftedResults <- vector()
rLOut <- vector()
readsOut <- vector()
fcOut <- vector()
percentIDOut <- vector()
SDIDOut <- vector()
foldCovOut <- vector()
foldCovSDOut <- vector()
results <- vector()

# Populate the results vector with the appropriate data
# 'AspA' is occasionally used as a negative control, and should not be included in 'results'
results <- graph[graph$target != 'aspA', ]
# find the unique kmers present in the study
uniqueReadLengths <- unique(results$readLength)
uniqueFoldCoverages <- unique(results$foldCoverage)

for (readlength in uniqueReadLengths) {
  # for each iteration, sifted results only contains data of the current read length
  for (foldcoverage in uniqueFoldCoverages) {
    # can manually set the acceptable values for fold coverage to be used in the analysis
    # for instance: & results$foldCoverage >= 10.0
    siftedResults <- results[results$readLength == readlength & results$foldCoverage == foldcoverage, ]
    # populate the vectors
    rLOut <- append(rLOut, readlength, after = length(rLOut))
    # Note that this is only for this script - the fold-coverage values are estimates, while the
    # number of reads is exact.  
    readLengthTitle <- foldcoverage * 100000
    readsOut <- append(readsOut, readLengthTitle, after = length(readsOut))
    fcOut <- append(fcOut, foldcoverage, after = length(fcOut))
    percentIDOut <- append(percentIDOut, mean(siftedResults$MedianPercentID), after = length(percentIDOut))
    SDIDOut <- append(SDIDOut, sd(siftedResults$MedianPercentID), after = length(SDIDOut))
    foldCovOut <- append(foldCovOut, mean(siftedResults$MedianFoldCoverage), after = length(foldCovOut))
    foldCovSDOut <- append(foldCovSDOut, mean(siftedResults$FoldCoverageSD), after = length(foldCovSDOut))
  }
}

# populate the data frame to be used in generating the figure
output <- data.frame(readLength = rLOut, numReads = readsOut, foldCoverage = fcOut, meanpercentID = percentIDOut, SDPercentID = SDIDOut, meanFoldCoverage = foldCovOut, foldCovSD = foldCovSDOut)

# output is printed to console - allows visual inspection
outputOrdered <- output[order(output$foldCoverage, output$readLength),]
print(outputOrdered)

# plot a nice output graph using ggplot2
ggplot(outputOrdered, aes(readLength,meanpercentID)) +
  geom_point() + theme_bw() + theme(axis.text=element_text(size=10)) +
  geom_errorbar(aes(ymin=meanpercentID-SDPercentID, ymax=meanpercentID+SDPercentID)) +
  xlab("Read Length") + ylab("Percent Identity") +
  scale_y_continuous(breaks=c(50,60, 70, 80, 90, 100)) +
  # Facet wrap allows for the tiled multi-figure graphs created by this script
  facet_wrap( ~ numReads)

filename <- "Figure1_multipanel_15by19.pdf"

# saves the chart as 'filename', in with dimensions provided below
ggsave(filename, width=19, height=15, units=c("cm"))

# Note that the titles still need to be edited in order to make them "pretty"

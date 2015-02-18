#!/usr/bin/env Rscript

# Figure 3

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
graph <- read.csv(fileofinterest, header=TRUE)

# Initialise the required variables
siftedResults <- vector()
kmerOut <- vector()
FCOut <- vector()
percentIDOut <- vector()
SDIDOut <- vector()
foldCovOut <- vector()
foldCovSDOut <- vector()
results <- vector()

# Populate the results vector with the appropriate data
# 'AspA' is occasionally used as a negative control, and should not be included in 'results'
# the foldCoverage can also be set to '>=' a value
#  for instabce: & graph$foldCoverage == 2.0
results <- graph[graph$target != 'aspA', ]
# find the unique kmers present in the study
uniqueKmers <- unique(results$kmerLength)
uniquefoldcoverages <- unique(results$foldCoverage)
# perform the analysis for each kmer length
for (kmerlength in uniqueKmers) {
  # for each iteration, sifted results only contains data of the current kmer
  for (foldcoverage in uniquefoldcoverages) {
    siftedResults <- results[results$kmerLength == kmerlength & results$foldCoverage == foldcoverage, ]
    # populate the vectors
    # Since the titles are going to be added to the top of each plot, they must be formatted appropriately
    if (foldcoverage %% 1 == 0) {
      foldcoverageTitle <- paste(foldcoverage, ".0 Fold-Coverage", sep="")
    } else {
      foldcoverageTitle <- paste(foldcoverage, " Fold-Coverage", sep="")
    }
    FCOut <- append(FCOut, foldcoverageTitle, after = length(FCOut))
    kmerOut <- append(kmerOut, kmerlength, after = length(kmerOut))
    percentIDOut <- append(percentIDOut, mean(siftedResults$MedianPercentID), after = length(percentIDOut))
    SDIDOut <- append(SDIDOut, sd(siftedResults$MedianPercentID), after = length(SDIDOut))
    # fold coverage isn't present in the figure, but if there's odd results, these data might be 
    # useful in finding out where the issue is
    foldCovOut <- append(foldCovOut, mean(siftedResults$MedianFoldCoverage), after = length(foldCovOut))
    foldCovSDOut <- append(foldCovSDOut, mean(siftedResults$FoldCoverageSD), after = length(foldCovSDOut))
  }
}

# populate the data frame to be used in generating the figure
output <- data.frame(kmerLength = kmerOut, FC = FCOut, meanpercentID = percentIDOut, SDPercentID = SDIDOut, meanFoldCoverage = foldCovOut, foldCovSD = foldCovSDOut)

# output is printed to console - allows visual inspection
print(output)

# plot a nice output graph using ggplot2
sp <- ggplot(output, aes(x=kmerLength, y=meanpercentID))
# add scatterplot and a simple theme
sp + geom_point() + theme_bw() + theme(axis.text=element_text(size=10)) +
  # error bars
  geom_errorbar(aes(ymin=meanpercentID-SDPercentID, ymax=meanpercentID+SDPercentID), width=.2) +
  # axis labels - ylim could probably be determined automatically, but it's not. Needs to be set manually
  xlab("Kmer Length (bp)") + ylab("Percent Identity") + ylim(35, 120) +
  scale_y_continuous(breaks=c(30, 40, 50,60, 70, 80, 90, 100)) +
  facet_wrap( ~ FC)

filename <- "Figure3_multipanel_15by19.pdf"

# saves the chart as 'filename', in with dimensions provided below
ggsave(filename, width=19, height=15, units=c("cm"))

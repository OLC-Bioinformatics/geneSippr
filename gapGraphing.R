#!/usr/bin/env Rscript

# Figure 5

library("ggplot2", lib.loc="/usr/local/lib/R/site-library")

# If this script is called from an external script, then comment and uncomment the appropriate
# cmdArgs, name, outPath, and csvFile lines as necessary

# The arguments from the calling script are passed in using the commandArgs function
# these arguments include the path to use as well as the name of the strain
cmdArgs <- commandArgs(trailingOnly=TRUE)

# The name will be used often enough that it needed its own variable
name <- cmdArgs[2]

# get the path
path <- cmdArgs[1]
setwd(path)

# Create a string of the name and the necessary additions for the full name of the csv file
csvFile <- paste(path, name, "_gaps.csv", sep='')

# Read the tables of insert sizes
a = read.csv(csvFile, header=TRUE)

# Create a plot of GapSize vs GapFrequency
p <- ggplot(a, aes(x=GapSize,y=GapFrequency)) 

# Histogram
p + geom_histogram(stat="identity") + xlab("Gap Size (bp)") + ylab("Gap Frequency") + 
  # This will have to be changed if the scale changes 
  # Uncomment if you want a scale with ticks every five bp instead of one bp
  #   + scale_x_continuous(breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60))
  # I'm thinking something like finding the maximum gap size, and then either making a variable something like range(0,maxGap,5) with 5 being the step
  # size - either that or calculating the step size to be 1/10 
  theme_bw() +
  theme(axis.text.x = element_text(size=10))

# Set the filename
filename <- paste(name, "_insert_sizes.pdf", sep="")

# saves the chart as 'filename', in with dimensions provided below
ggsave(filename, width=19, height=7)
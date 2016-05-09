#!/usr/bin/env Rscript
#./diveRsity.R input=nameOfGenePopFile output=nameOfROutputFile
library(diveRsity)
options(warn=-1)
for(i in commandArgs()){
	tmp = strsplit(i, "=")
	if(tmp[[1]][1] == "input"){input = tmp[[1]][2]}
	if(tmp[[1]][1] == "output"){output = tmp[[1]][2]}
}
a=diffCalc(input, fst=T, pairwise=F, outfile=output)


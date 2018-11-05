#!/usr/bin/Rscript
library(cobs)
#read data of each galaxy's redshift and overdensity
df <- read.table('~/ZF_environ/nearest_neighbor/newstat/overdensity_redshift.txt',header=FALSE)
z <- df[,1] 
log1pd <- log10(df[,2])
cobsxy05 <- cobs(z,log1pd, ic='BIC', tau=0.05) 
cobsxy25 <- cobs(z,log1pd, ic='BIC', tau=0.25) #25%-tile 
cobsxy50 <- cobs(z,log1pd, ic='BIC', tau=0.5)  #median 
cobsxy75 <- cobs(z,log1pd, ic='BIC', tau=0.75) #75%-tile
#write down the redshift, 25%-tile, median, and 75%-tile of each galaxy to a file
write.table(data.frame(cobsxy25$x,cobsxy25$fitted,cobsxy50$fitted,cobsxy75$fitted,df[,3]),'~/ZF_environ/nearest_neighbor/newstat/cobs_output.txt',sep="\t", quote = FALSE)



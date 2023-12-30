# load the climatol package
library(climatol) 

# set working directory as needed in every case
setwd('C:/Users/Administrator/Documents/CLINO')

# read the stations file
sta=read.table('SIstation.txt',sep='\t',header=TRUE)

# add extension txt to the file names
sta$CODE=sprintf('%s.txt',sta$CODE)

# write new stations file
write.table(sta,'stations.txt',sep='\t',row.names=FALSE)
rclimdex2climatol('stations.txt',stcol=c(1,3,2,4,5),kvar=1,chrcod=c(1,7))
rclimdex2climatol('stations.txt',stcol=c(1,3,2,4,5),kvar=2,chrcod=c(1,7))
rclimdex2climatol('stations.txt',stcol=c(1,3,2,4,5),kvar=3,chrcod=c(1,7))

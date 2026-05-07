#!/user/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(plyr)

d=read.table(args[1])
options(scipen=999)
data=NULL
bin=args[2]

for(i in 1:nrow(d))
{
  #i=grep(paste0("chr",chrom,"$"),d$V1)
  chr=as.character(d[i,1])
  print(chr)
  print(d[i,2])
  print(round_any(d[i,2],as.numeric(bin),f=ceiling))
  to=round_any(d[i,2],as.numeric(bin),f=ceiling)
  start=seq(as.numeric(0),to,by=as.numeric(bin))
  end=start+as.numeric(bin)
  #print(chr)
  #print(start)
  pos=paste(chr,start,end,sep="_")
  dd=data.frame(chr,start,end,pos)
  print(nrow(dd))
  data=rbind(data,dd)
}
#Going to generate 530992 bins
fileName=args[3]
write.table(data,file=fileName,col.names = F,row.names = F,quote = F,sep="\t")


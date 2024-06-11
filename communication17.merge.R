#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")


library(limma)
library(tidyverse)
outFile="./merge.txt"     
  

geneList=list()
data = read_tsv("./GSE176171 2.txt", col_names = TRUE)
rt = read.table("./GSE176171 2.txt", header=T, sep="\t", check.names=F)
header = unlist(strsplit(inputFile, "\\.|\\-|\\_"))
if("TP53" %in% row.names(rt) | "GAPDH" %in% row.names(rt) | "CCL4" %in% row.names(rt)){
  geneList[[header[1]]]=row.names(rt)
}else{
  geneList[[header[1]]]=as.vector(rt[,1])
}

interGenes = Reduce(intersect, geneList)


allTab=data.frame()
for(i in 1:length(files)){
    inputFile=files[i]
    if(inputFile==outFile){next}
    header=unlist(strsplit(inputFile, "\\.|\\-|\\_"))
   
    rt=read.table(inputFile, header=T, sep="\t", check.names=F)
    if("TP53" %in% row.names(rt) | "GAPDH" %in% row.names(rt) | "CCL4" %in% row.names(rt)){
    	rt=avereps(rt)
    }else{
	    rt=as.matrix(rt)
	    rownames(rt)=rt[,1]
	    exp=rt[,2:ncol(rt)]
	    dimnames=list(rownames(exp),colnames(exp))
	    data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	    rt=avereps(data)
	}
    colnames(rt)=paste0(header[1], "_", colnames(rt))

    if(ncol(allTab)==0){
    	allTab=rt[interGenes,]
    }else{
    	allTab=cbind(allTab, rt[interGenes,])
    }
}

outTab=rbind(geneNames=colnames(allTab), allTab)
write.table(outTab, file=outFile, sep="\t", quote=F, col.names=F)


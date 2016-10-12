#Create a file to build the network
setwd("media/tandrean/Elements/PhD/Natalia_Gut_Development/motif_enrichment/paired_end/motif")

#import the set of motif and proteins bounded by the splicing factor predicted by http://rbpmap.technion.ac.il/
a <- read.table("motif.txt",stringsAsFactors=F)
b <- read.table("RNA_binding_proteins_defenitive.txt",stringsAsFactors=F)
d <- matrix(0, nrow=94,ncol=134)
a <- a$V1
s <- b$V1
s <- t(s)
rownames(d) <- s
colnames(d) <- a
length(a)
length(s)
d[is.na(d)] <- 0
head(d)

#Import the set of proteins with relative motifs predicted
motif_scanned <- scan("reference_motif.txt",what="chr")
head(motif_scanned,40)
str(motif_scanned)
idx <- grep("chr",motif_scanned)
for (i in 1:length(idx)){
  if(i<length(idx))
    d[motif_scanned[(idx[i]+1):(idx[i+1]-1)],motif_scanned[idx[i]]] = 1
  else
    d[motif_scanned[(idx[i]+1):length(motif_scanned)],motif_scanned[idx[i]]] = 1
}
head(d)
head(d)
write.table(d,"network_table_definitive.tsv",row.names=T,col.names=T,quote=F,sep="\t")

#Map network in protein space
network <- d%*%t(d)
head(network[1:10,1:10])

#Map network in motif space
nework_motif <- t(d)%*%(d)
head(nework_motif)

#Save the files
save(d,network,nework_motif, file = "splice_network.RData")

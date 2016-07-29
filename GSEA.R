a <- rep(c(0),each=8939)
gene <- read.table(file="GeneName.txt",header=TRUE)
names(a) <- gene$x

down <- read.table(file="DownKOover2xGeneCountover6Noprespore.txt",sep="\t",header=TRUE)

b <- rep(c(1),each=668)
names(b) <- down$V1

c <- a[-which(names(a) %in% names(b))]
d <- a[which(names(a) %in% names(b))]

genes <- c(c,b)
genes1 <- genes[sort(names(genes))]

library(goseq)
load("CountBias.RData")
countbias1 <- countbias[sort(names(countbias))]
pdf("PlotFit-GSEA.pdf")
pwf <- nullp(genes1,bias.data=countbias1,plot.fit=TRUE)
dev.off()
load("GoSets.RData")
pvals <- goseq(pwf,gene2cat=go.sets.dd)
library(GO.db)

 # Significant cateories
        sign <- (pvals[pvals$over_represented_pvalue < 0.05, ])
        write.table(sign,file="SignCatPvalueF.txt",sep="\t")

        a <- read.table(file="OntologyUsed",header=TRUE)

	xx <- as.data.frame(GOTERM)
	total <- xx[which(xx$go_id %in% sign$category),]
	total1 <- total[!duplicated(total[,c('go_id','go_id')]),]	
	total2 <- total1[,c(1,3,4,5,6)]
	
	Final <- cbind(sign,total2[match(sign$category,total2$go_id), ])
	write.table(Final,file="Combine-.txt",sep="\t")	

      #  onlySig <- a[which(a$GO %in% sign$category),]
       # for(go in unique(onlySig$GO)){
       # sink("SigncatPvalueTerms.txt",append=TRUE)
       # print(GOTERM[[go]])
       # cat("--------------------------------------\n")
       # sink()
       # }

###### Multiple Hypothesis correction @ GSEA

        enriched.Go <- pvals$category[p.adjust(pvals$over_represented_pvalue,method="BH") < .05]

	En <- Final[which(Final$category %in% enriched.Go),]
	write.table(En,file="EnrichedGo-DownKOover2xGeneCountover6Noprespore",sep="\t")



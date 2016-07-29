	files <- list.files(pattern = ".count$")
	data = lapply(files, function(x)read.table(x))	
	count = as.data.frame(data)
	CountDF <- count[,-c(3,5,7,9,11,13,15)]
	rownames(CountDF) <- CountDF[,1]
	CountDF <- CountDF[,-1]
	colnames(CountDF) <- c("DGC-KO1-A","DGC-KO1-B","DGC-KO2-A","DGC-KO2-B","DGC-RI1-A","DGC-RI1-B","DGC-RI2-A","DGC-RI2-B")	
	CountDF <- CountDF[-c(13863:13867),]	
	save(CountDF,file="DictyEnsembleCount.RData")	

	library(edgeR)
	targets <- read.table(file="targets",header=TRUE)
	#load("DictyEnsembleCount.RData")

	### Filtering of read on the basis of rowsum

	keep <- rowSums(cpm(CountDF)>1) >=4
	y <- CountDF[keep,]
	nrow(CountDF)
	nrow(y)
	
	#### Correlation Plot between Samples
	library(ape)
        t <- cor(y,method="spearman")
        hc <- hclust(dist(1-t))
        pdf("CorrelationSamples.pdf")
        plot.phylo(as.phylo(hc),type="p",edge.col="4",edge.width="2",show.node.label=TRUE,no.margin=TRUE)
        dev.off()
	
	### DEG
	
	group <- factor(targets$Condition,levels=c("Control","KO"))
	design <- model.matrix(~group)	

	y <- DGEList(counts=y,group=group)
	y <- calcNormFactors(y)
	write.table(y$counts,file="FilterGenesRawCount.txt",sep="\t")
        #### just to check the group and coulmn name symmetry
        #### data.frame(samples = colnames(y),group)

	pdf("PlotMDS.pdf")
	plotMDS(y)
	dev.off()

	rownames(design) <- colnames(y)	
	y <- estimateGLMCommonDisp(y,design,verbose=TRUE)
	y <- estimateGLMTrendedDisp(y,design)
	y <- estimateGLMTagwiseDisp(y,design)
	fit <- glmFit(y,design)
	diff <- glmLRT(fit)
	save(diff,file="DiffExpGenes.RData")

	##### DEG at FDR < 0.05 using Benjamin Hochberg method. BH is thw default method in decideTestsDGE
	summary(de <- decideTestsDGE(diff,p.value=0.05,adjust.method="BH"))
	detags <- rownames(y)[as.logical(de)]
	pdf("LogFCvsAvelogCPM.pdf")
	plotSmear(diff, de.tags=detags)
	abline(h=c(-1, 1), col="blue")
	dev.off()
	DiffPvalue <- diff$table[which(rownames(diff$table) %in% detags),]
	write.table(DiffPvalue,file="DEGPvalue.txt",sep="\t")

	##### Extracting non-DEGs
	NotDiffPvalue <- diff$table[-which(rownames(diff$table) %in% detags),]
	write.table(NotDiffPvalue,file="Not-DifferentialExpressed/NotDiffPvalue.txt",sep="\t",quote=FALSE)
	
	###### Heatmap of significamtaly enriches Down/Up regulated genes

        down <-  DiffPvalue[which( DiffPvalue$logFC < 0), ]
        y1 <- cpm(y,prior.count=2, log=TRUE)
        countD <- y1[rownames(down), ]
        #y1 <- cpm(countD[1:40,], prior.count=2, log=TRUE)
        pdf("DownHeatmap40.pdf")
        heatmap(countD[1:40,],margins = c(10, 10))
        dev.off()

        up <-  DiffPvalue[which( DiffPvalue$logFC > 0), ]
        countU <- y1[rownames(up), ]
        #y1 <- cpm(countU[1:40,], prior.count=2, log=TRUE)
        pdf("UpHeatmap40.pdf")
        heatmap(countU[1:40,],margins = c(10, 10))
        dev.off()

        counts.per.m <- cpm(y, normalized.lib.sizes=TRUE)
        downNorC <- counts.per.m[rownames(down), ]
        write.table(downNorC,file="DownGenesNormalizedCount.txt",sep="\t")
        ### or counts.per.m <- cpm(y) gives the same results as the above

        upNorC <- counts.per.m[rownames(up), ]
        write.table(upNorC,file="UpGenesNormalizedCount.txt",sep="\t")

	unchanged <- counts.per.m[rownames(NotDiffPvalue),]
	write.table(unchanged,file="Not-DifferentialExpressed/NotDiff-NormCount.txt",sep="\t",quote=FALSE)

        write.table(counts.per.m,file="FilterGenesNormalizedCount.txt",sep="\t")

	

	######## Gene Set Enrichment Analysis using goseq

	genes=as.integer(p.adjust(diff$table$PValue[diff$table$logFC!=0],method="BH")<.05)
        names(genes)=row.names(diff$table[diff$table$logFC!=0,])
        library(goseq)
	countbias <- rowSums(y$counts)[rowSums(y$counts) != 0]
        pdf("PlotFit-GSEA.pdf")
        pwf <- nullp(genes,bias.data=countbias,plot.fit=TRUE)
        dev.off()
        load("GoSets.RData")
	pvals <- goseq(pwf,gene2cat=go.sets.dd) 
	library(GO.db)
	for(go in pvals$category){
	sink("AllEnrichedGoCategory.txt",append=TRUE)
	print(GOTERM[[go]])	
	cat("--------------------------------------\n")
	sink()
	}
        write.table(pvals,file="EnrichedGeneCategory.txt",sep="\t")

	# Significant cateories
	sign <- (pvals[pvals$over_represented_pvalue < 0.05, ])
	write.table(sign,file="SignCatPvalueF.txt",sep="\t")

	a <- read.table(file="OntologyUsed",header=TRUE)	

	onlySig <- a[which(a$GO %in% sign$category),]
	for(go in unique(onlySig$GO)){
	sink("SigncatPvalueTerms.txt",append=TRUE)
	print(GOTERM[[go]])
        cat("--------------------------------------\n")
	sink()
	}
	
	lfc <- DiffPvalue[which(rownames(DiffPvalue) %in% onlySig$ID),]
	write.table(lfc,file="DEG-OnlySignGoCat",sep="\t")	
	
	###### Multiple Hypothesis correction @ GSEA

	enriched.Go <- pvals$category[p.adjust(pvals$over_represented_pvalue,method="BH") < .05]
	onlyEnr <- a[which(a$GO %in% enriched.Go),]
	lfc1 <- DiffPvalue[which(rownames(DiffPvalue) %in% onlyEnr$ID),]
	write.table(lfc1,file="SigniGenesTwoSignCategory.txt",sep="\t")

